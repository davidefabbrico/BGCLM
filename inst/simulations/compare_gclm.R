# ============================================================
# BGCLM (spike-and-slab bayesiano) vs gclm (LASSO frequentista, Varando et al.)
#
# Uso:  Rscript inst/simulations/compare_gclm.R
#
# Entrambi i metodi risolvono B Sigma + Sigma B^T + C = 0 con la stessa
# convenzione (gclm::clyap e BGCLM::lyapunov concordano), quindi le matrici B
# sono confrontabili voce per voce senza trasposizioni.
# ============================================================

suppressMessages({
  library(BGCLM)
  library(gclm)
})

# ---- configurazione ----------------------------------------------------
N_REP     <- 10
N_VALS    <- c(100, 500, 2000)
P         <- 5
N_EDGE    <- 5
N_ITER    <- 8000
BURNIN    <- 4000
THIN      <- 4
LAMBDAS   <- exp(seq(log(0.6), log(0.004), length.out = 30))
SCENARIOS <- c("gershgorin", "generic", "weak", "hetero_C", "marginal")
OUT_RDS   <- "inst/simulations/compare_gclm_results.rds"

# ========================================================================
# Generatori
#
# Ogni scenario restituisce il vero B, la vera Sigma, il pattern di archi e la
# vera C. "dom" segnala se B e' diagonalmente dominante, cioe' se il vero cade
# dentro la classe che BGCLM sa rappresentare.
# ========================================================================

is_dominant <- function(B) all(-diag(B) > rowSums(abs(B)) - abs(diag(B)) + 1e-10)

rand_off <- function(p, n_edge, lo, hi) {
  Boff <- matrix(0, p, p)
  Boff[rbind(c(2, 1), c(3, 2), c(1, 3))] <- c(1, -1, 1) * runif(3, lo, hi)
  for (k in seq_len(n_edge - 3)) {
    free <- which(diag(p) == 0 & Boff == 0 & t(Boff) == 0, arr.ind = TRUE)
    if (!nrow(free)) break
    Boff[free[sample(nrow(free), 1), , drop = FALSE]] <-
      sample(c(-1, 1), 1) * runif(1, lo, hi)
  }
  Boff
}

# Diagonale di Gershgorin: il vero e' dentro la classe del modello bayesiano.
gen_gershgorin <- function(p, n_edge, lo = 0.4, hi = 0.8, a = 0.5) {
  Boff <- rand_off(p, n_edge, lo, hi)
  B <- build_B(Boff, rep(a, p))
  list(B = B, C = diag(p), delta = (Boff != 0) * 1L, Sigma = lyapunov(B, diag(p)))
}

# B stabile generico: la diagonale e' uno shift spettrale, non la somma delle
# righe, quindi in generale NON e' diagonalmente dominante. E' il test che
# conta: il vero cade fuori dalla parametrizzazione di BGCLM.
gen_generic <- function(p, n_edge, lo = 0.4, hi = 0.8, a = 0.5) {
  Boff <- rand_off(p, n_edge, lo, hi)
  B <- Boff - (max(Re(eigen(Boff, only.values = TRUE)$values)) + a) * diag(p)
  list(B = B, C = diag(p), delta = (Boff != 0) * 1L, Sigma = lyapunov(B, diag(p)))
}

gen_weak <- function(p, n_edge) gen_gershgorin(p, n_edge, lo = 0.10, hi = 0.25)

# C vera diagonale ma non identita': entrambi i metodi la fissano a I.
gen_hetero_C <- function(p, n_edge, a = 0.5) {
  Boff <- rand_off(p, n_edge, 0.4, 0.8)
  B <- build_B(Boff, rep(a, p))
  cd <- exp(runif(p, log(0.35), log(2.8)))
  cd <- cd / prod(cd)^(1 / p)            # tr geometrica 1: toglie la deriva di scala
  C <- diag(cd)
  list(B = B, C = C, delta = (Boff != 0) * 1L, Sigma = lyapunov(B, C))
}

# Modello marginalizzato: si genera su p + 3 nodi e se ne osservano p. Le
# latenti inducono una C non diagonale ("blunt edges"), quindi C = I e' falsa
# per costruzione. Il B bersaglio e' il sotto-blocco, un target solo
# approssimato: le metriche sugli archi qui vanno lette con cautela.
gen_marginal <- function(p, n_edge, n_lat = 3, a = 0.5) {
  pf <- p + n_lat
  Boff <- rand_off(pf, n_edge + 2 * n_lat, 0.4, 0.8)
  Bf <- build_B(Boff, rep(a, pf))
  Sf <- lyapunov(Bf, diag(pf))
  o <- seq_len(p)
  S <- Sf[o, o]; B <- Bf[o, o]
  list(B = B, C = -(B %*% S + S %*% t(B)), delta = (Boff[o, o] != 0) * 1L,
       Sigma = S, marginal = TRUE)
}

GEN <- list(gershgorin = gen_gershgorin, generic = gen_generic, weak = gen_weak,
            hetero_C = gen_hetero_C, marginal = gen_marginal)

# ========================================================================
# Stimatori
# ========================================================================

# Il campionatore C++ usa Y'Y/n come statistica sufficiente (media nota = 0),
# non cov(Y): gclm deve ricevere esattamente la stessa matrice.
shat <- function(Y) crossprod(Y) / nrow(Y)

# lungo il path il LASSO puo' produrre un B non stabile: lyapunov() fa stop().
safe_lyap <- function(B, C) tryCatch(lyapunov(B, C), error = function(e) NULL)

loglik <- function(S, Sh, n) {
  if (is.null(S)) return(-Inf)
  R <- tryCatch(chol(S), error = function(e) NULL)
  if (is.null(R)) return(-Inf)
  Ri <- backsolve(R, diag(ncol(S)))
  -0.5 * n * (ncol(S) * log(2 * pi) + 2 * sum(log(diag(R))) + sum((Ri %*% t(Ri)) * Sh))
}

fit_bgclm <- function(Y) {
  t0 <- proc.time()[3]
  f <- gclm_mcmc(Y, n_iter = N_ITER, burnin = BURNIN, thin = THIN, verbose = FALSE)
  f$time <- proc.time()[3] - t0
  f
}

# free_C = FALSE riproduce l'assunzione di BGCLM (C = I fissa): lambdac = -1
# tiene C bloccata. free_C = TRUE stima C diagonale con ridge verso C0 = 1, che
# e' anche cio' che ancora il raggio di scala Sigma(aB, aC) = Sigma(B, C).
fit_gclm <- function(Y, free_C = FALSE) {
  p <- ncol(Y); n <- nrow(Y); Sh <- shat(Y)
  t0 <- proc.time()[3]
  path <- tryCatch(
    gclm.path(Sh, lambdas = LAMBDAS, C = rep(1, p),
              lambdac = if (free_C) 0.01 else -1,
              eps = 1e-8, maxIter = 2000),
    error = function(e) NULL)
  if (is.null(path)) return(NULL)

  off <- which(diag(p) == 0)
  ok  <- vapply(path, function(z) !is.null(z$B) && all(is.finite(z$B)), logical(1))
  path <- path[ok]
  if (!length(path)) return(NULL)

  # punteggio per la ROC: il lambda piu' grande a cui l'arco e' ancora attivo.
  # E' l'ordine d'ingresso lungo il path, il duale frequentista del PIP.
  score <- rep(0, length(off))
  for (z in path) {
    act <- abs(z$B[off]) > 1e-8
    score[act] <- pmax(score[act], z$lambda)
  }

  bic <- vapply(path, function(z) {
    S <- safe_lyap(z$B, diag(z$C, nrow = p))
    df <- sum(abs(z$B[off]) > 1e-8) + p + if (free_C) p - 1 else 0
    -2 * loglik(S, Sh, n) + log(n) * df
  }, numeric(1))
  bic[!is.finite(bic)] <- Inf
  if (all(is.infinite(bic))) return(NULL)
  best <- path[[which.min(bic)]]

  list(B = best$B, Sigma = safe_lyap(best$B, diag(best$C, nrow = p)),
       score = score, lambda = best$lambda, time = proc.time()[3] - t0)
}

# ========================================================================
# Metriche
#
# cov95 e' separata fra archi veri e non archi: mediata su tutti gli
# off-diagonali sarebbe dominata dagli zeri, il cui intervallo marginale e'
# degenere in [0, 0] e contiene lo zero vero per costruzione.
# ========================================================================

metrics <- function(B_hat, Sigma_hat, score, truth, qs = NULL) {
  p <- ncol(truth$B); off <- which(diag(p) == 0)
  lab <- truth$delta[off]
  if (is.null(Sigma_hat)) Sigma_hat <- matrix(NA_real_, p, p)
  sel <- if (is.null(qs)) as.integer(abs(B_hat[off]) > 1e-8) else as.integer(score > 0.5)
  out <- data.frame(
    auc     = if (all(is.na(score))) NA_real_ else auc(score, lab),
    tpr     = if (any(lab == 1)) mean(sel[lab == 1] == 1) else NA_real_,
    fpr     = if (any(lab == 0)) mean(sel[lab == 0] == 1) else NA_real_,
    rmse_B  = sqrt(mean((B_hat[off] - truth$B[off])^2)),
    err_Sig = norm(Sigma_hat - truth$Sigma, "F") / norm(truth$Sigma, "F"),
    cov_e   = NA_real_, cov_0 = NA_real_,
    n_edge  = sum(sel))
  if (!is.null(qs)) {
    inside <- truth$B[off] >= qs[1, , ][off] & truth$B[off] <= qs[2, , ][off]
    out$cov_e <- if (any(lab == 1)) mean(inside[lab == 1]) else NA_real_
    out$cov_0 <- if (any(lab == 0)) mean(inside[lab == 0]) else NA_real_
  }
  out
}

# ========================================================================
# Loop principale
# ========================================================================

res <- list()
grid <- expand.grid(rep = seq_len(N_REP), n = N_VALS, scen = SCENARIOS,
                    stringsAsFactors = FALSE)
cat(sprintf("%d configurazioni (%d scenari x %d n x %d repliche)\n\n",
            nrow(grid), length(SCENARIOS), length(N_VALS), N_REP))

for (g in seq_len(nrow(grid))) {
  scen <- grid$scen[g]; n <- grid$n[g]; rep_i <- grid$rep[g]
  set.seed(1000 * match(scen, SCENARIOS) + 10 * rep_i)

  truth <- GEN[[scen]](P, N_EDGE)
  Y <- rmvn0(n, truth$Sigma)

  fb <- fit_bgclm(Y)
  qs <- apply(fb$B, c(1, 2), stats::quantile, probs = c(0.025, 0.975))
  rows <- cbind(method = "BGCLM",
                metrics(fb$B_mean, fb$Sigma_mean, fb$PIP[which(diag(P) == 0)],
                        truth, qs = qs),
                time = fb$time)

  for (fc in c(FALSE, TRUE)) {
    fg <- fit_gclm(Y, free_C = fc)
    nm <- if (fc) "gclm (C diag. stimata)" else "gclm (C = I)"
    rows <- rbind(rows, if (is.null(fg))
      cbind(method = nm, metrics(matrix(NA, P, P), matrix(NA, P, P),
                                 rep(NA_real_, P * (P - 1)), truth), time = NA)
      else cbind(method = nm, metrics(fg$B, fg$Sigma, fg$score, truth),
                 time = fg$time))
  }

  res[[g]] <- cbind(scen = scen, n = n, rep = rep_i,
                    dom = is_dominant(truth$B),
                    offC = max(abs(truth$C[diag(P) == 0])),
                    ratC = max(diag(truth$C)) / min(diag(truth$C)), rows)

  if (g %% 10 == 0 || g == nrow(grid))
    cat(sprintf("  %3d / %d\n", g, nrow(grid)))
}

res <- do.call(rbind, res)
dir.create(dirname(OUT_RDS), showWarnings = FALSE, recursive = TRUE)
saveRDS(res, OUT_RDS)

# ========================================================================
# Riepilogo
# ========================================================================

num <- c("auc", "tpr", "fpr", "rmse_B", "err_Sig", "cov_e", "cov_0", "n_edge", "time")
agg <- aggregate(res[num], by = res[c("scen", "n", "method")],
                 FUN = function(x) mean(x, na.rm = TRUE))
agg <- agg[order(match(agg$scen, SCENARIOS), agg$n, agg$method), ]

info <- aggregate(res[c("dom", "offC", "ratC")], by = res["scen"],
                  FUN = function(x) mean(as.numeric(x)))
cat("\n=== scenari ===\n")
cat("dom  = frazione di B veri diagonalmente dominanti (1 = dentro la classe BGCLM)\n")
cat("offC = massima voce fuori diagonale della C vera (0 = C diagonale)\n")
cat("ratC = rapporto max/min della diagonale di C (1 = C proporzionale a I)\n")
cat("entrambi i metodi assumono C = I, quindi offC > 0 o ratC > 1 = misspecificato\n\n")
print(info, row.names = FALSE, digits = 3)

cat("\n=== medie su", N_REP, "repliche ===\n")
cat("cov_e / cov_0 = copertura al 95% su archi veri / non archi (solo BGCLM)\n")
cat("n_edge = archi selezionati (veri:", N_EDGE, "su", P * (P - 1), ")\n\n")
for (s in SCENARIOS) {
  cat("---", s, "---\n")
  print(agg[agg$scen == s, setdiff(names(agg), "scen")], row.names = FALSE, digits = 3)
  cat("\n")
}
cat("risultati grezzi salvati in", OUT_RDS, "\n")
