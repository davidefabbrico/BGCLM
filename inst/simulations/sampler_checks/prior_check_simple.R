# ============================================================
# Test del campionatore per grafi semplici (gclm_mcmc(..., simple = TRUE)).
# (a) verosimiglianza spenta: la catena deve restituire la prior sugli stati
#     delle coppie e sul numero di archi E;
# (b) nessun draw con coppie reciproche;
# (c) ogni B campionata e' stabile.
# Uso:  Rscript inst/simulations/sampler_checks/prior_check_simple.R
# ============================================================
suppressMessages(library(BGCLM))
p <- 5; P <- p * (p - 1) / 2
pairs <- which(upper.tri(diag(p)), arr.ind = TRUE)

# distribuzione esatta: E | rho ~ Binom(P, rho), rho ~ Beta(a, b) => beta-binomiale;
# stato di una coppia: P(nessuno) = b/(a+b), P(i->j) = P(j->i) = a/(2(a+b))
exact_E <- function(a, b) sapply(0:P, function(k) choose(P, k) * beta(a + k, b + P - k) / beta(a, b))

for (hp in list(c(1, 1), c(2, 5))) {
  set.seed(11)
  f <- gclm_mcmc(matrix(rnorm(50 * p), 50), n_iter = 220000, burnin = 20000, thin = 10,
                 a_rho = hp[1], b_rho = hp[2], simple = TRUE, prior_only = TRUE, verbose = FALSE)
  K <- dim(f$delta)[3]
  E <- apply(f$delta, 3, sum)
  st <- t(apply(f$delta, 3, function(d) {
    s <- apply(pairs, 1, function(ij) if (d[ij[2], ij[1]] == 1) 1 else if (d[ij[1], ij[2]] == 1) 2 else 0)
    tabulate(s + 1, 3) / nrow(pairs) }))
  essE <- ess(E)
  cat(sprintf("\n==== a_rho = %g, b_rho = %g   (%d draw, ESS di E = %.0f) ====\n", hp[1], hp[2], K, essE))
  ex <- exact_E(hp[1], hp[2]); em <- tabulate(E + 1, P + 1) / K
  print(round(rbind(E = 0:P, empirica = em, esatta = ex), 3))
  cat(sprintf("media di E: %.3f (esatta %.3f, errore MC ~ %.3f)\n", mean(E), sum((0:P) * ex), sd(E) / sqrt(essE)))
  cat(sprintf("stati delle coppie  nessuno / i->j / j->i:  empirici %.3f / %.3f / %.3f   esatti %.3f / %.3f / %.3f\n",
      mean(st[, 1]), mean(st[, 2]), mean(st[, 3]),
      hp[2] / sum(hp), hp[1] / (2 * sum(hp)), hp[1] / (2 * sum(hp))))
  cat(sprintf("rho medio %.3f (esatto %.3f);  margini a_i: media %.3f (1), varianza %.3f (0.5)\n",
      mean(f$rho), hp[1] / sum(hp), mean(f$a), var(as.vector(f$a))))
  rs <- t(apply(f$B, 3, function(b) rowSums(abs(b)) - abs(diag(b))))
  cat(sprintf("max |cor(a_i, somma della riga)| = %.3f\n", max(abs(diag(cor(f$a, rs))))))
  cat(sprintf("(b) draw con coppie reciproche: %d\n", sum(apply(f$delta, 3, function(d) any(d == 1 & t(d) == 1)))))
  cat(sprintf("(c) draw instabili: %d\n", sum(apply(f$B, 3, function(B) max(Re(eigen(B, only.values = TRUE)$values)) >= 0))))
  cat(sprintf("accettazione inversione %.3f, nascite rifiutate per 2-ciclo %d\n", f$acc$rev, f$acc$blocked))
}
