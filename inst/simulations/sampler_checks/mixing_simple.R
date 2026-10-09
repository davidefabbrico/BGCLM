# Mixing del campionatore per grafi semplici a n crescente, contro la posteriori
# libera filtrata sui grafi semplici. ll_gap < 0 = la catena resta su grafi che
# si adattano peggio del vero: segnale che non esplora.
suppressMessages(library(BGCLM))
src <- readLines("inst/simulations/compare_gclm.R")
eval(parse(text = paste(src[1:(grep("^GEN <- list", src) - 1)], collapse = "\n")))
off <- which(diag(P) == 0)
llg <- function(S, Sh, n) { R <- chol(S); Ri <- backsolve(R, diag(P))
  -0.5 * n * (P * log(2 * pi) + 2 * sum(log(diag(R))) + sum((Ri %*% t(Ri)) * Sh)) }
met <- function(f, tr, Y) {
  sel <- (f$PIP > 0.5) * 1L; sk <- function(d) (d | t(d)) & upper.tri(d)
  pr <- which(sk(tr$delta), arr.ind = TRUE)
  data.frame(auc = auc(f$PIP[off], tr$delta[off]), grafo_esatto = all(sel == tr$delta),
    orient = mean(apply(pr, 1, function(ij) sel[ij[2], ij[1]] == tr$delta[ij[2], ij[1]] && sel[ij[1], ij[2]] == tr$delta[ij[1], ij[2]])),
    P_vero = mean(apply(f$delta, 3, function(d) all(d == tr$delta))),
    ll_gap = mean(f$loglik) - llg(tr$Sigma, crossprod(Y) / nrow(Y), nrow(Y)))
}
res <- list()
for (n in c(500, 2000, 10000)) for (r in 1:10) {
  set.seed(1000 + 10 * r); tr <- gen_gershgorin(P, N_EDGE); Y <- rmvn0(n, tr$Sigma)
  set.seed(77 + r); fn <- gclm_mcmc(Y, n_iter = N_ITER, burnin = BURNIN, thin = THIN, simple = TRUE, verbose = FALSE)
  set.seed(77 + r); fl <- gclm_mcmc(Y, n_iter = N_ITER, burnin = BURNIN, thin = THIN, verbose = FALSE)
  res[[length(res) + 1]] <- cbind(n = n, rep = r, versione = "nuovo (prior coppie, nascite vietate)", met(fn, tr, Y))
  res[[length(res) + 1]] <- cbind(n = n, rep = r, versione = "libera + filtro", met(restrict_simple(fl, 50), tr, Y))
}
res <- do.call(rbind, res); saveRDS(res, "inst/simulations/sampler_checks/mixing_simple.rds")
print(aggregate(cbind(auc, grafo_esatto, orient, P_vero, ll_gap) ~ n + versione, res, mean), digits = 3)
cat("\nrepliche con ll_gap < -2 (catena peggiore del vero di piu' di 2 nat):\n")
print(aggregate(ll_gap ~ n + versione, res, function(z) sum(z < -2)))
