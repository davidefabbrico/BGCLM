# ============================================================
# Diagonale uguale (B_ii = -d) contro diagonale di Gershgorin libera, modello A1
# con parallel tempering in entrambi i casi.
#
#  corretto: il vero ha diagonale uguale, d = max_i sum_j |B_ij| + 0.5 (stessi
#            archi e stessi semi di gershgorin / random); si stimano entrambi i
#            modelli.
#  misspec:  il vero e' quello degli scenari originali (a_i = 0.5, diagonale
#            disuguale); si stima solo il modello a diagonale uguale. Il modello
#            a diagonale libera sugli stessi dati e' in results/data/simple_vs_free.csv.
#
# Uso:  Rscript inst/simulations/equal_diag.R <corretto|misspec> <scenario>
# Scrive inst/simulations/results/data/equal_diag_<modo>_<scenario>.csv
# ============================================================
suppressMessages(library(BGCLM))
src <- readLines("inst/simulations/compare_gclm.R")
eval(parse(text = paste(src[1:(grep("^GEN <- list", src) - 1)], collapse = "\n")))
a <- commandArgs(TRUE); MODE <- a[1]; SC <- a[2]
NV <- c(500, 2000, 10000)
D_SHAPE <- 2; D_RATE <- 1                                  # d ~ Gamma(2, 1)
OUT <- sprintf("inst/simulations/results/data/equal_diag_%s_%s.csv", MODE, SC)
off <- which(diag(P) == 0)

ll_true <- function(Y, S) { Sh <- crossprod(Y) / nrow(Y); R <- chol(S)
  -0.5 * nrow(Y) * (P * log(2 * pi) + 2 * sum(log(diag(R))) + sum(chol2inv(R) * Sh)) }
metrics <- function(f, tr, Y) {
  if (is.null(f))                                           # nessun draw senza coppie reciproche
    return(data.frame(auc = NA, tpr = NA, fpr = NA, grafo_esatto = NA, orient_ok = NA, P_grafo_vero = NA,
                      d_media = NA, d_vero = NA, frac_draw = 0, ll_gap = NA, sec = NA))
  lab <- tr$delta[off]; sel <- select_graph(f)
  data.frame(auc = auc(f$PIP[off], lab), tpr = mean(sel[off][lab == 1]), fpr = mean(sel[off][lab == 0]),
             grafo_esatto = all(sel == tr$delta),
             orient_ok = mean(sel[tr$delta == 1] == 1),
             P_grafo_vero = mean(apply(f$delta, 3, function(d) all(d == tr$delta))),
             d_media = if (isTRUE(f$equal_diag)) mean(f$d) else NA_real_,
             d_vero = if (diff(range(diag(tr$B))) < 1e-12) -tr$B[1, 1] else NA_real_,
             frac_draw = f$frac_kept, ll_gap = mean(f$loglik) - ll_true(Y, tr$Sigma), sec = f$sec)
}
fit <- function(Y, eq) {
  t0 <- proc.time()[3]
  f <- gclm_mcmc(Y, n_iter = N_ITER_PT, burnin = BURNIN_PT, thin = THIN_PT, simple = TRUE,
                 eps_bridge = EPS_BRIDGE, pt_betas = PT_BETAS, equal_diag = eq,
                 a_shape = if (eq) D_SHAPE else 2, a_rate = if (eq) D_RATE else 2, verbose = FALSE)
  f$sec <- proc.time()[3] - t0; f
}

res <- list()
for (n in NV) for (r in 1:10) {
  set.seed(1000 * match(SC, SCENARIOS) + 10 * r)
  tr <- get(paste0("gen_", SC))(P, N_EDGE)
  if (MODE == "corretto") {                                 # stessi archi, diagonale uguale
    Bo <- tr$B; diag(Bo) <- 0
    tr$B <- Bo - (max(rowSums(abs(Bo))) + 0.5) * diag(P); tr$Sigma <- lyapunov(tr$B, diag(P))
  }
  Y <- rmvn0(n, tr$Sigma)
  for (eq in if (MODE == "corretto") c(FALSE, TRUE) else TRUE) {
    set.seed(77 + r)
    f <- tryCatch(suppressWarnings(fit(Y, eq)), error = function(e) NULL)
    res[[length(res) + 1]] <- cbind(modo = MODE, scenario = SC, n = n, rep = r,
                                    modello = if (eq) "diagonale uguale" else "diagonale libera", metrics(f, tr, Y))
  }
}
res <- do.call(rbind, res)
write.csv(res, OUT, row.names = FALSE)
num <- c("P_grafo_vero", "grafo_esatto", "orient_ok", "tpr", "fpr", "auc", "d_media", "d_vero", "frac_draw", "ll_gap", "sec")
ag <- aggregate(res[num], by = res[c("n", "modello")], FUN = function(z) mean(z, na.rm = TRUE))
print(ag[order(ag$n, ag$modello), ], row.names = FALSE, digits = 3)
