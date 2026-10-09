# Verifica del parallel tempering su tutte le repliche: per ogni replica 4
# catene (2 dal grafo vuoto, 2 dal grafo vero) con 10 temperature.
# Uso: Rscript inst/simulations/sampler_checks/pt_multichain.R <A1|libero> <n>
suppressMessages(library(BGCLM))
args <- commandArgs(TRUE); MOD <- args[1]; N <- as.integer(args[2])
src <- readLines("inst/simulations/compare_gclm.R")
eval(parse(text = paste(src[1:(grep("^GEN <- list", src) - 1)], collapse = "\n")))
off <- which(diag(P) == 0)
rhat <- function(ch) { nn <- min(lengths(ch)); x <- sapply(ch, function(z) tail(z, nn))
  W <- mean(apply(x, 2, var)); Bv <- nn * var(colMeans(x)); sqrt(((nn - 1) / nn * W + Bv / nn) / W) }
res <- list()
for (r in 1:10) {
  set.seed(1000 + 10 * r); tr <- gen_gershgorin(P, N_EDGE); Y <- rmvn0(N, tr$Sigma)
  fits <- lapply(1:4, function(k) { set.seed(500 + 10 * r + k)
    suppressWarnings(gclm_mcmc(Y, n_iter = N_ITER_PT, burnin = BURNIN_PT, thin = THIN_PT,
      simple = (MOD == "A1"), eps_bridge = EPS_BRIDGE, pt_betas = PT_BETAS,
      B_init = if (k > 2) tr$B else NULL, verbose = FALSE)) })
  pip <- sapply(fits, function(f) f$PIP[off])
  pv <- sapply(fits, function(f) mean(apply(f$delta, 3, function(d) all(d == tr$delta))))
  res[[r]] <- data.frame(modello = MOD, n = N, rep = r,
    diff_vuoto_vs_vero = max(abs(rowMeans(pip[, 1:2]) - rowMeans(pip[, 3:4]))),
    diff_fra_4 = max(apply(pip, 1, function(z) diff(range(z)))),
    P_vero_vuoto = mean(pv[1:2]), P_vero_dal_vero = mean(pv[3:4]),
    rhat_E = rhat(lapply(fits, function(f) apply(f$delta, 3, sum))),
    rhat_loglik = rhat(lapply(fits, `[[`, "loglik")),
    scambio_min = min(sapply(fits, function(f) min(f$swap))))
}
res <- do.call(rbind, res)
write.csv(res, sprintf("inst/simulations/sampler_checks/pt_multichain_%s_%d.csv", MOD, N), row.names = FALSE)
print(res, row.names = FALSE, digits = 3)
