# Identificabilita' debole. Per ogni replica in cui il grafo selezionato (regola
# per coppia, MPPI salvate da simple_vs_free.R) e' diverso dal vero, quanto si
# avvicina un drift della classe Gershgorin con quel supporto alla Sigma vera?
# KL = min_B KL(N(0, Sigma_vero) || N(0, Sigma(B))): n * KL e' il vantaggio atteso
# del vero in log-verosimiglianza. Se il grafo selezionato contiene il vero il
# minimo e' 0 (basta porre a zero gli archi in piu') e non si ottimizza.
# Uso: Rscript inst/simulations/identifiability/05_near_equivalence.R
source("inst/simulations/identifiability/common_ident.R")
pv <- readRDS("inst/simulations/results/data/simple_vs_free_pip.rds")

kl <- function(S0, S) 0.5 * (sum(diag(solve(S, S0))) - p +
                             determinant(S)$modulus - determinant(S0)$modulus)
best_kl <- function(S0, supp, starts = 20) {
  idx <- which(supp == 1)
  f <- function(th) {
    Bo <- matrix(0, p, p); Bo[idx] <- th[seq_along(idx)]
    S <- tryCatch(lyapunov(build_B(Bo, exp(th[-seq_along(idx)])), diag(p)), error = function(e) NULL)
    if (is.null(S)) return(1e6)
    v <- tryCatch(as.numeric(kl(S0, S)), error = function(e) 1e6)
    if (is.finite(v)) v else 1e6
  }
  min(sapply(seq_len(starts), function(s) { set.seed(s)
    optim(c(rnorm(length(idx), 0, 0.5), log(runif(p, 0.2, 1))), f, method = "BFGS",
          control = list(maxit = 2000, reltol = 1e-12))$value }))
}

res <- list()
for (z in Filter(function(z) z$n >= 2000, pv)) {
  tr <- true_model(z$scen, z$rep)
  s <- select_graph(list(PIP = z$A1, simple = TRUE))
  if (all(s == tr$delta)) next
  sopra <- all(s[tr$delta == 1] == 1)                     # contiene il vero
  k <- if (sopra) 0 else best_kl(tr$Sigma, s)
  res[[length(res) + 1]] <- data.frame(scenario = z$scen, n = z$n, rep = z$rep,
    grafo_sel = edges_lab(s), archi_sel = sum(s), contiene_vero = sopra,
    kl = k, n_kl = z$n * k, n_per_10_nat = if (k > 0) 10 / k else Inf)
}
res <- do.call(rbind, res)
write.csv(res, "inst/simulations/identifiability/output/05_near_equivalence.csv", row.names = FALSE)
print(res, row.names = FALSE, digits = 3)
