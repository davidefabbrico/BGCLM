# Taratura della scala delle temperature: tasso di scambio fra livelli adiacenti
# per K = 6, 8, 10 temperature geometriche fra 1 e 0.01, a n = 2000 e 10000.
suppressMessages(library(BGCLM))
src <- readLines("inst/simulations/compare_gclm.R")
eval(parse(text = paste(src[1:(grep("^GEN <- list", src) - 1)], collapse = "\n")))
for (n in c(2000, 10000)) {
  set.seed(1000 + 10 * 1); tr <- gen_gershgorin(P, N_EDGE); Y <- rmvn0(n, tr$Sigma)
  for (K in c(6, 8, 10)) {
    b <- 0.01^(seq(0, 1, length.out = K)); set.seed(9); t0 <- proc.time()[3]
    f <- gclm_mcmc(Y, n_iter = 10000, burnin = 2000, thin = 2, simple = TRUE, eps_bridge = 0.1, pt_betas = b, verbose = FALSE)
    cat(sprintf("n=%5d K=%2d  scambi: %s   (min %.2f)  [%.0fs]\n", n, K, paste(sprintf("%.2f", f$swap), collapse = " "), min(f$swap), proc.time()[3] - t0))
  }
}
