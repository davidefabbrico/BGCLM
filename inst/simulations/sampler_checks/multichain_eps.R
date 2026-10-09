# Catene multiple per il campionatore con augmentazione: 2 catene dal grafo
# vuoto e 2 dal grafo vero (semi diversi). Se mescolano, arrivano alla stessa
# posteriori: MPPI, probabilita' del grafo vero e log-verosimiglianza devono
# coincidere entro l'errore Monte Carlo, e R-hat deve essere vicino a 1.
# Uso: Rscript inst/simulations/sampler_checks/multichain_eps.R <eps> [n_iter] [thin]
suppressMessages(library(BGCLM))
args <- commandArgs(TRUE)
EPS <- as.numeric(args[1])
src <- readLines("inst/simulations/compare_gclm.R")
eval(parse(text = paste(src[1:(grep("^GEN <- list", src) - 1)], collapse = "\n")))
off <- which(diag(P) == 0)
NIT <- if (length(args) >= 2) as.integer(args[2]) else N_ITER
THN <- if (length(args) >= 3) as.integer(args[3]) else THIN

rhat <- function(chains) {                     # Gelman-Rubin su catene di lunghezza diversa
  m <- length(chains); nn <- min(lengths(chains)); x <- sapply(chains, function(z) tail(z, nn))
  W <- mean(apply(x, 2, var)); Bv <- nn * var(colMeans(x))
  sqrt(((nn - 1) / nn * W + Bv / nn) / W)
}
res <- list()
for (n in c(2000, 10000)) for (r in 1:5) {
  set.seed(1000 + 10 * r); tr <- gen_gershgorin(P, N_EDGE); Y <- rmvn0(n, tr$Sigma)
  fits <- list()
  for (k in 1:4) {
    set.seed(500 + 10 * r + k)
    fits[[k]] <- suppressWarnings(gclm_mcmc(Y, n_iter = NIT, burnin = BURNIN, thin = THN,
      simple = TRUE, eps_bridge = EPS, B_init = if (k > 2) tr$B else NULL, verbose = FALSE))
  }
  pip <- lapply(fits, function(f) f$PIP[off])
  pv  <- sapply(fits, function(f) mean(apply(f$delta, 3, function(d) all(d == tr$delta))))
  ll  <- sapply(fits, function(f) mean(f$loglik))
  pip_vuoto <- (pip[[1]] + pip[[2]]) / 2; pip_vero <- (pip[[3]] + pip[[4]]) / 2
  res[[length(res) + 1]] <- data.frame(n = n, rep = r,
    max_diff_MPPI_vuoto_vs_vero = max(abs(pip_vuoto - pip_vero)),
    max_diff_MPPI_fra_4 = max(sapply(1:20, function(e) diff(range(sapply(pip, `[`, e))))),
    P_vero_vuoto = mean(pv[1:2]), P_vero_dal_vero = mean(pv[3:4]),
    ll_vuoto = mean(ll[1:2]), ll_dal_vero = mean(ll[3:4]),
    rhat_loglik = rhat(lapply(fits, `[[`, "loglik")),
    rhat_E = rhat(lapply(fits, function(f) apply(f$delta, 3, sum))),
    frac_tenuti = mean(sapply(fits, `[[`, "frac_kept")))
}
res <- do.call(rbind, res)
write.csv(res, sprintf("inst/simulations/sampler_checks/multichain_eps_%s_iter%d.csv", EPS, NIT), row.names = FALSE)
options(width = 200); print(res, row.names = FALSE, digits = 3)
