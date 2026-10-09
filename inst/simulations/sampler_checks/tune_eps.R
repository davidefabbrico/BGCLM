# Taratura di eps (augmentazione con stato "entrambi gli archi").
# Per ogni eps: frazione di draw tenuti (senza 2-cicli), ESS della
# log-verosimiglianza e del numero di archi sui soli draw tenuti, recupero del
# grafo e ll_gap = media della log-verosimiglianza (draw tenuti) - quella del
# vero. ll_gap molto negativo = catena bloccata su grafi peggiori del vero.
suppressMessages(library(BGCLM))
src <- readLines("inst/simulations/compare_gclm.R")
eval(parse(text = paste(src[1:(grep("^GEN <- list", src) - 1)], collapse = "\n")))
off <- which(diag(P) == 0)
llg <- function(S, Sh, n) { R <- chol(S); Ri <- backsolve(R, diag(P))
  -0.5 * n * (P * log(2 * pi) + 2 * sum(log(diag(R))) + sum((Ri %*% t(Ri)) * Sh)) }
met <- function(f, tr, Y) {
  sel <- (f$PIP > 0.5) * 1L
  data.frame(frac_tenuti = if (is.null(f$frac_kept)) 1 else f$frac_kept,
    draw_tenuti = length(f$loglik), ess_loglik = ess(f$loglik), ess_E = ess(apply(f$delta, 3, sum)),
    auc = auc(f$PIP[off], tr$delta[off]), grafo_esatto = all(sel == tr$delta),
    P_vero = mean(apply(f$delta, 3, function(d) all(d == tr$delta))),
    ll_gap = mean(f$loglik) - llg(tr$Sigma, crossprod(Y) / nrow(Y), nrow(Y)))
}
res <- list()
for (n in c(500, 2000, 10000)) for (r in 1:10) {
  set.seed(1000 + 10 * r); tr <- gen_gershgorin(P, N_EDGE); Y <- rmvn0(n, tr$Sigma)
  for (e in c(0, 0.01, 0.05, 0.1)) {
    set.seed(77 + r)
    f <- tryCatch(suppressWarnings(gclm_mcmc(Y, n_iter = N_ITER, burnin = BURNIN, thin = THIN,
                  simple = TRUE, eps_bridge = e, verbose = FALSE)), error = function(z) NULL)
    res[[length(res) + 1]] <- if (is.null(f))
      data.frame(n = n, rep = r, eps = e, frac_tenuti = 0, draw_tenuti = 0, ess_loglik = NA, ess_E = NA,
                 auc = NA, grafo_esatto = NA, P_vero = NA, ll_gap = NA)
    else cbind(n = n, rep = r, eps = e, met(f, tr, Y))
  }
  set.seed(77 + r)
  f <- gclm_mcmc(Y, n_iter = N_ITER, burnin = BURNIN, thin = THIN, verbose = FALSE)
  res[[length(res) + 1]] <- cbind(n = n, rep = r, eps = NA, met(f, tr, Y))
}
res <- do.call(rbind, res); res$versione <- ifelse(is.na(res$eps), "libero (simple=FALSE)", paste0("eps = ", res$eps))
saveRDS(res, "inst/simulations/sampler_checks/tune_eps.rds")
ag <- aggregate(res[c("frac_tenuti", "draw_tenuti", "ess_loglik", "ess_E", "auc", "grafo_esatto", "P_vero", "ll_gap")],
                by = res[c("n", "versione")], FUN = function(z) mean(z, na.rm = TRUE))
ag$repliche_senza_draw <- aggregate(res["frac_tenuti"], by = res[c("n", "versione")], FUN = function(z) sum(z == 0))$frac_tenuti
options(width = 200); print(ag[order(ag$n, ag$versione), ], row.names = FALSE, digits = 3)
cat("\nll_gap minimo per cella (peggior replica):\n")
print(aggregate(res["ll_gap"], by = res[c("n", "versione")], FUN = function(z) min(z, na.rm = TRUE)), row.names = FALSE, digits = 3)
