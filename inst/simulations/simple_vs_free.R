# ============================================================
# Recupero del grafo: modello libero (2-cicli ammessi) contro modello A1
# (grafi semplici, prior sulle coppie, augmentazione eps), entrambi con
# parallel tempering.
#
# Uso:  Rscript inst/simulations/simple_vs_free.R [scenario n]
# Scrive inst/simulations/results/data/simple_vs_free[_scenario_n].csv
# ============================================================

suppressMessages(library(BGCLM))
src <- readLines("inst/simulations/compare_gclm.R")
eval(parse(text = paste(src[1:(grep("^GEN <- list", src) - 1)], collapse = "\n")))
OUT <- "inst/simulations/results/data/simple_vs_free.csv"
a <- commandArgs(TRUE)
FIT_FREE <- !("--senza-libero" %in% a); a <- setdiff(a, "--senza-libero")
SC <- c("gershgorin", "weak"); NV <- c(500, 2000, 10000)
if (length(a) == 2) { SC <- a[1]; NV <- as.integer(a[2]); OUT <- sub(".csv", sprintf("_%s_%s.csv", a[1], a[2]), OUT, fixed = TRUE) }
off <- which(diag(P) == 0)

graph_metrics <- function(f, tr) {
  lab <- tr$delta[off]; sel <- select_graph(f)   # coppia per A1, freccia per il libero
  sk <- function(d) (d | t(d)) & upper.tri(d)
  pairs <- which(sk(tr$delta), arr.ind = TRUE)
  data.frame(
    auc          = auc(f$PIP[off], lab),
    tpr          = mean(sel[off][lab == 1]),
    fpr          = mean(sel[off][lab == 0]),
    grafo_esatto = all(sel == tr$delta),
    orient_ok    = mean(apply(pairs, 1, function(ij)
                     sel[ij[2], ij[1]] == tr$delta[ij[2], ij[1]] &&
                     sel[ij[1], ij[2]] == tr$delta[ij[1], ij[2]])),
    P_grafo_vero = mean(apply(f$delta, 3, function(d) all(d == tr$delta))),
    frac_draw    = if (isTRUE(f$simple)) f$frac_kept else 1,
    ess_loglik   = ess(f$loglik))
}

res <- list(); pips <- list()
for (scen in SC) for (n in NV) for (r in 1:10) {
  set.seed(1000 * match(scen, SCENARIOS) + 10 * r)
  tr <- get(paste0("gen_", scen))(P, N_EDGE); Y <- rmvn0(n, tr$Sigma)
  set.seed(77 + r)
  f <- if (FIT_FREE) gclm_mcmc(Y, n_iter = N_ITER_PT, burnin = BURNIN_PT, thin = THIN_PT,
                 pt_betas = PT_BETAS, verbose = FALSE)
  set.seed(77 + r)
  s <- gclm_mcmc(Y, n_iter = N_ITER_PT, burnin = BURNIN_PT, thin = THIN_PT, simple = TRUE,
                 eps_bridge = EPS_BRIDGE, pt_betas = PT_BETAS, verbose = FALSE)
  if (FIT_FREE) res[[length(res) + 1]] <- cbind(scenario = scen, n = n, rep = r, posterior = "libera", graph_metrics(f, tr))
  res[[length(res) + 1]] <- cbind(scenario = scen, n = n, rep = r, posterior = "A1", graph_metrics(s, tr))
  pips[[length(pips) + 1]] <- list(scen = scen, n = n, rep = r, delta = tr$delta, A1 = s$PIP, free = f$PIP)
}
res <- do.call(rbind, res)
dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE)
write.csv(res, OUT, row.names = FALSE)
saveRDS(pips, sub(".csv", "_pip.rds", OUT, fixed = TRUE))
num <- c("frac_draw", "ess_loglik", "auc", "tpr", "fpr", "grafo_esatto", "orient_ok", "P_grafo_vero")
ag <- aggregate(res[num], by = res[c("scenario", "n", "posterior")], FUN = mean)
print(ag[order(ag$scenario, ag$n, ag$posterior), ], row.names = FALSE, digits = 2)
