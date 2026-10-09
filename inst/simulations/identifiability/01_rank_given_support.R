# Fase 2.1 - A1: identificabilita' di B dato il supporto.
# Sul supporto vero, B' = B + K Sigma^-1 ha gli stessi zeri fuori diagonale se e
# solo se (K Sigma^-1)_ji = 0 su tutte le posizioni zero del grafo vero. Se la
# matrice di questo sistema ha rango p(p-1)/2, l'unica soluzione e' K = 0.
# Uso: Rscript inst/simulations/identifiability/01_rank_given_support.R
source("inst/simulations/identifiability/common_ident.R")
res <- list()
for (scen in names(SEEDS)) for (r in 1:10) {
  tr <- true_model(scen, r)
  M <- fiber_matrix(tr$Sigma); Z <- which(tr$delta[OFF] == 0)
  sv <- svd(M[Z, ])$d
  res[[length(res) + 1]] <- data.frame(scenario = scen, rep = r, grafo = edges_lab(tr$delta),
    zeri = length(Z), incognite = ncol(M), rango = qr(M[Z, ])$rank,
    sv_min = min(sv), cond = max(sv) / min(sv))
}
res <- do.call(rbind, res)
write.csv(res, "inst/simulations/identifiability/output/01_rank_given_support.csv", row.names = FALSE)
print(res, row.names = FALSE, digits = 3)
cat(sprintf("\nrango = p(p-1)/2 = %d in %d casi su %d\n", nrow(PAIRS), sum(res$rango == nrow(PAIRS)), nrow(res)))
