# ============================================================
# Figure per il confronto BGCLM / gclm.
#
# Uso:  Rscript inst/simulations/plot_results.R
#
# Produce in inst/simulations/figures/:
#   graph_correct_spec.png grafo a posteriori del modello A1, gershgorin e weak (n = 500)
#   graphs_by_n.png        modello A1, gershgorin a n crescente
#   auc_coverage.png       AUC e copertura vs n, dai risultati gia' salvati
# ============================================================

suppressMessages({
  library(BGCLM)
  library(gclm)
})

SRC <- "inst/simulations/compare_gclm.R"
RDS <- "inst/simulations/compare_gclm_results.rds"
FIG <- "inst/simulations/figures"
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

# riusa generatori e configurazione senza rieseguire le simulazioni
src <- readLines(SRC)
eval(parse(text = paste(src[1:(grep("^GEN <- list", src) - 1)], collapse = "\n")))
GEN <- list(gershgorin = gen_gershgorin, generic = gen_generic, weak = gen_weak,
            hetero_C = gen_hetero_C, marginal = gen_marginal)

fit_one <- function(scen, n, rep_i = 1) {
  set.seed(1000 * match(scen, SCENARIOS) + 10 * rep_i)
  truth <- GEN[[scen]](P, N_EDGE)
  Y <- rmvn0(n, truth$Sigma)
  s <- gclm_mcmc(Y, n_iter = N_ITER_PT, burnin = BURNIN_PT, thin = THIN_PT, simple = TRUE,
                 eps_bridge = EPS_BRIDGE, pt_betas = PT_BETAS, verbose = FALSE)
  list(fit = s, truth = truth)
}

# ---- 1. un grafo per scenario, a n = 500 ------------------------------
cat("grafi per scenario...\n")
graph_legend <- function() {
  plot(NA, xlim = 0:1, ylim = 0:1, axes = FALSE, xlab = "", ylab = "")
  legend("center", bty = "n", cex = 0.95, seg.len = 2.4,
         legend = c("true edge (halo)", "high MPPI", "low MPPI", "negative coefficient"),
         col = c("#E39B2E", gray(0), gray(0.65), gray(0)),
         lwd = c(5, 3.5, 1.2, 2), lty = c(1, 1, 1, 2))
}
fits500 <- lapply(setNames(c("gershgorin", "weak"), c("gershgorin", "weak")), function(s) fit_one(s, 500))
png(file.path(FIG, "graph_correct_spec.png"), width = 2100, height = 800, res = 210)
par(mfrow = c(1, 3))
for (s in names(fits500))
  plot_gclm_graph(fits500[[s]]$fit, truth = fits500[[s]]$truth, main = sprintf("%s  (n = 500)", s))
graph_legend()
dev.off()

# ---- 2. lo stesso scenario a n crescente ------------------------------
cat("grafi a n crescente...\n")
fitsn <- lapply(N_VALS, function(n) fit_one("gershgorin", n))
png(file.path(FIG, "graphs_by_n.png"), width = 2100, height = 800, res = 210)
par(mfrow = c(1, 3))
for (k in seq_along(N_VALS))
  plot_gclm_graph(fitsn[[k]]$fit, truth = fitsn[[k]]$truth, main = sprintf("gershgorin, n = %d", N_VALS[k]))
dev.off()

# ---- 3. AUC e copertura vs n -----------------------------------------
if (file.exists(RDS)) {
  cat("curve riassuntive...\n")
  r <- readRDS(RDS)
  ag <- aggregate(r[c("auc", "cov_e")], by = r[c("scen", "n", "method")],
                  FUN = function(z) mean(z, na.rm = TRUE))
  meths <- c("BGCLM (A1)", "BGCLM", "gclm (C = I)", "gclm (C diag. stimata)")
  meths <- meths[meths %in% ag$method]
  cols  <- c("BGCLM (A1)" = gray(0), "BGCLM" = gray(0.6),
             "gclm (C = I)" = "#C0504D", "gclm (C diag. stimata)" = "#4F81BD")[meths]
  labs  <- c("BGCLM (A1)" = "BGCLM (A1)", "BGCLM" = "BGCLM (unrestricted)",
             "gclm (C = I)" = "gclm (C = I)", "gclm (C diag. stimata)" = "gclm (C diag.)")[meths]

  png(file.path(FIG, "auc_coverage.png"), width = 2100, height = 1000, res = 200)
  par(mfrow = c(2, length(RUN_SCEN)), mar = c(3.4, 3.6, 2.2, 0.6), mgp = c(2.2, 0.7, 0))

  for (s in RUN_SCEN) {
    d <- ag[ag$scen == s, ]
    plot(NA, xlim = range(N_VALS), ylim = c(0.4, 1), log = "x",
         xlab = "n", ylab = "AUC", main = s)
    abline(h = 0.5, col = gray(0.75), lty = 3)
    for (m in meths) {
      dm <- d[d$method == m, ]
      dm <- dm[order(dm$n), ]
      lines(dm$n, dm$auc, col = cols[m], lwd = 2, type = "b", pch = 16, cex = 0.8)
    }
    if (s == RUN_SCEN[1])
      legend("bottomright", col = cols[meths], lwd = 2, legend = labs,
             pch = 16, bty = "n", cex = 0.62, seg.len = 1.6)
  }
  for (s in RUN_SCEN) {
    plot(NA, log = "x", xlim = range(N_VALS), ylim = c(0, 1), xlab = "n",
         ylab = "95% coverage (true edges)", main = s)
    abline(h = 0.95, col = "#C0504D", lty = 2)
    for (m in intersect(c("BGCLM (A1)", "BGCLM"), meths)) {
      d <- ag[ag$scen == s & ag$method == m, ]; d <- d[order(d$n), ]
      lines(d$n, d$cov_e, type = "b", pch = 16, lwd = 2, col = cols[m])
    }
  }
  dev.off()
} else {
  cat("salto le curve: manca", RDS, "- esegui prima compare_gclm.R\n")
}

cat("figure in", FIG, "\n")
