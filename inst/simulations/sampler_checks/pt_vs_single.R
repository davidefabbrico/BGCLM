# Serve il parallel tempering? Modello A1 con catena singola a parita' di costo
# (100.000 iterazioni = 10 livelli x 10.000) contro il tempering. Per ogni
# replica 4 catene, 2 dal grafo vuoto e 2 dal grafo vero: se la catena mescola,
# arrivano alla stessa posteriori. Metriche del grafo sulla catena 1 (regola per
# coppia).
# Uso: Rscript inst/simulations/sampler_checks/pt_vs_single.R <pt|singola> <scenario> <n>
suppressMessages(library(BGCLM))
args <- commandArgs(TRUE); MODE <- args[1]; SC <- args[2]; N <- as.integer(args[3])
src <- readLines("inst/simulations/compare_gclm.R")
eval(parse(text = paste(src[1:(grep("^GEN <- list", src) - 1)], collapse = "\n")))
off <- which(diag(P) == 0)
rhat <- function(ch) { nn <- min(lengths(ch)); x <- sapply(ch, function(z) tail(z, nn))
  W <- mean(apply(x, 2, var)); Bv <- nn * var(colMeans(x)); sqrt(((nn - 1) / nn * W + Bv / nn) / W) }

fit <- function(Y, B0) {
  t0 <- proc.time()[3]
  f <- if (MODE == "pt")
    gclm_mcmc(Y, n_iter = N_ITER_PT, burnin = BURNIN_PT, thin = THIN_PT, simple = TRUE,
              eps_bridge = EPS_BRIDGE, pt_betas = PT_BETAS, B_init = B0, verbose = FALSE)
  else
    gclm_mcmc(Y, n_iter = 100000, burnin = 20000, thin = 20, simple = TRUE,
              eps_bridge = EPS_BRIDGE, B_init = B0, verbose = FALSE)
  f$time <- proc.time()[3] - t0; f
}

res <- list()
for (r in 1:10) {
  set.seed(1000 * match(SC, SCENARIOS) + 10 * r)
  tr <- get(paste0("gen_", SC))(P, N_EDGE); Y <- rmvn0(N, tr$Sigma)
  fits <- lapply(1:4, function(k) { set.seed(500 + 10 * r + k)
    suppressWarnings(fit(Y, if (k > 2) tr$B else NULL)) })
  pip <- sapply(fits, function(f) f$PIP[off])
  pv  <- sapply(fits, function(f) mean(apply(f$delta, 3, function(d) all(d == tr$delta))))
  f1 <- fits[[1]]; sel <- select_graph(f1); lab <- tr$delta[off]
  res[[r]] <- data.frame(modo = MODE, scenario = SC, n = N, rep = r,
    diff_vuoto_vs_vero = max(abs(rowMeans(pip[, 1:2]) - rowMeans(pip[, 3:4]))),
    diff_fra_4 = max(apply(pip, 1, function(z) diff(range(z)))),
    P_vero_vuoto = mean(pv[1:2]), P_vero_dal_vero = mean(pv[3:4]),
    rhat_E = rhat(lapply(fits, function(f) apply(f$delta, 3, sum))),
    rhat_loglik = rhat(lapply(fits, `[[`, "loglik")),
    auc = auc(f1$PIP[off], lab), tpr = mean(sel[off][lab == 1]), fpr = mean(sel[off][lab == 0]),
    grafo_esatto = all(sel == tr$delta), P_grafo_vero = pv[1],
    frac_draw = mean(sapply(fits, `[[`, "frac_kept")),
    sec = mean(sapply(fits, `[[`, "time")))
}
res <- do.call(rbind, res)
write.csv(res, sprintf("inst/simulations/sampler_checks/pt_vs_single_%s_%s_%d.csv", MODE, SC, N),
          row.names = FALSE)
print(res, row.names = FALSE, digits = 3)
