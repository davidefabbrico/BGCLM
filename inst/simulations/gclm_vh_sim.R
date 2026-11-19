# ============================================================
# gclm con la procedura di selezione di Varando & Hansen (2020, Sez. 5) negli
# scenari simulati, come terza variante accanto alle due con BIC di
# compare_gclm.R. Stessi dati (stessi semi), quindi BGCLM non va rifatto.
#
# Procedura: C diagonale (lambda_C = 0.01), path di 100 lambda fra 6e-4 e 6 sulle
# correlazioni del train; 200 divisioni casuali a meta'; su ogni divisione si
# ristima ogni supporto del path (lambda = 0) e si sceglie quello con la
# verosimiglianza migliore sul test; si tengono gli archi scelti in almeno
# l'85% delle divisioni. Sigma si ottiene ristimando B (supporto fisso, job = 10)
# e C diagonale sui dati completi. Punteggio per l'AUC: frequenza di selezione.
#
# Uso:  Rscript inst/simulations/gclm_vh_sim.R <scenario>
# Scrive inst/simulations/compare_gclm_results_vh_<scenario>.rds
# ============================================================
Sys.setenv(OMP_NUM_THREADS = "1")
suppressMessages({ library(BGCLM); library(gclm) })
src <- "inst/simulations/compare_gclm.R"
lines <- readLines(src)
eval(parse(text = paste(lines[1:(grep("^GEN <- list", lines) - 1)], collapse = "\n")))
for (e in parse(src)) if (is.call(e) && identical(e[[1]], as.name("<-")) && is.call(e[[3]]) &&
                          identical(e[[3]][[1]], as.name("function"))) eval(e)   # metrics, shat, ...
SC <- commandArgs(TRUE)[1]
VH_LAMBDAS <- 6 * 10^seq(-4, 0, length.out = 100); N_SPLIT <- as.integer(Sys.getenv("N_SPLIT", "200")); SOGLIA <- 0.85

vh_select <- function(Y) {
  p <- ncol(Y); freq <- matrix(0, p, p)
  for (r in seq_len(N_SPLIT)) {
    ix <- sample(nrow(Y), floor(nrow(Y) / 2))
    Str <- cov(Y[-ix, ]); dd <- diag(1 / sqrt(diag(Str))); Ste <- dd %*% cov(Y[ix, ]) %*% dd
    path <- tryCatch(gclm.path(cov2cor(Str), lambdas = VH_LAMBDAS, B = -0.5 * solve(cov2cor(Str)), C = rep(1, p),
                               lambdac = 0.01, eps = 1e-6, maxIter = 1000), error = function(e) NULL)
    if (is.null(path)) next
    path <- Filter(function(z) !is.null(z$B) && all(is.finite(z$B)), path)
    refit <- lapply(path, function(z) gclm(cov2cor(Str), B = z$B, C = z$C, lambda = 0, lambdac = -1, eps = 1e-10, job = 10))
    mll <- sapply(refit, function(z) { P <- tryCatch(solve(z$Sigma), error = function(e) NULL)
      if (is.null(P) || !all(is.finite(P))) Inf else -determinant(P, logarithm = TRUE)$modulus + sum(Ste * P) })
    s <- (abs(refit[[which.min(mll)]]$B) > 1e-8) * 1; diag(s) <- 0; freq <- freq + s
  }
  freq / N_SPLIT
}

res <- list()
grid <- expand.grid(rep = seq_len(N_REP), n = N_VALS)
for (g in seq_len(nrow(grid))) {
  n <- grid$n[g]; rep_i <- grid$rep[g]
  set.seed(1000 * match(SC, SCENARIOS) + 10 * rep_i)             # stessi dati di compare_gclm.R
  truth <- get(paste0("gen_", SC))(P, N_EDGE); Y <- rmvn0(n, truth$Sigma)
  set.seed(4000 + 1000 * match(SC, SCENARIOS) + 10 * rep_i + n)
  t0 <- proc.time()[3]
  freq <- vh_select(Y); sel <- (freq >= SOGLIA) * 1L
  Sh <- shat(Y)                                                   # ristima sui dati completi, supporto fisso
  fit <- tryCatch(gclm(Sh, B = -0.5 * diag(P) + 0.1 * sel, C = rep(1, P), lambda = 0, lambdac = 0.01,
                       eps = 1e-10, maxIter = 5000, job = 10), error = function(e) NULL)
  B_hat <- if (is.null(fit)) matrix(NA, P, P) else fit$B
  S_hat <- if (is.null(fit)) NULL else safe_lyap(fit$B, diag(fit$C, nrow = P))
  m <- metrics(B_hat, S_hat, freq[which(diag(P) == 0)], truth, sel = sel)
  res[[g]] <- cbind(scen = SC, n = n, rep = rep_i, dom = is_dominant(truth$B),
                    offC = max(abs(truth$C[diag(P) == 0])), ratC = max(diag(truth$C)) / min(diag(truth$C)),
                    method = "gclm (200 split)", m, time = proc.time()[3] - t0)
  if (g %% 10 == 0) message(sprintf("%s: %d / %d", SC, g, nrow(grid)))
}
res <- do.call(rbind, res)
saveRDS(res, sprintf("inst/simulations/compare_gclm_results_vh_%s.rds", SC))
print(aggregate(res[c("auc", "tpr", "fpr", "tpr_sk", "fpr_sk", "err_Sig", "n_edge", "time")], res["n"], mean), digits = 3)
