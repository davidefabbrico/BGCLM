# ============================================================
# Riassunto del confronto sul grafo di classe (class_comparison.R).
# Uso:  Rscript inst/simulations/class_collect.R
# Scrive inst/simulations/results/data/class_summary.csv e stampa le tabelle.
# ============================================================
IN <- "inst/simulations/results/data"
fl <- Sys.glob(file.path(IN, "class_*_*.rds")); fl <- fl[!grepl("summary", fl)]
tab <- do.call(rbind, lapply(fl, function(f) do.call(rbind, lapply(readRDS(f), `[[`, "row"))))
num <- c("auc_pair", "err_Sig", "cov_e", "cl_tpr_sk", "cl_fpr_sk", "cl_orient", "cl_overorient", "cl_shd", "cl_exact",
         "p_true_graph", "p_true_class", "time")
tab$cl_exact <- as.numeric(as.logical(as.character(tab$cl_exact))); tab$nontrivial <- as.logical(as.character(tab$nontrivial))
tab[num] <- lapply(tab[num], function(x) suppressWarnings(as.numeric(as.character(x))))
mn <- aggregate(tab[num], tab[c("scen", "n", "method")], function(x) mean(x, na.rm = TRUE))
sdv <- aggregate(tab["auc_pair"], tab[c("scen", "n", "method")], sd)
mn$auc_pair_sd <- sdv$auc_pair[match(paste(mn$scen, mn$n, mn$method), paste(sdv$scen, sdv$n, sdv$method))]
ord <- c("gershgorin", "weak", "random", "modules")
mn <- mn[order(match(mn$scen, ord), mn$n, mn$method), ]
write.csv(mn, file.path(IN, "class_summary.csv"), row.names = FALSE)
options(width = 200)
cat("== AUC sulle coppie (media, sd) ==\n")
print(reshape(mn[c("scen", "n", "method", "auc_pair", "auc_pair_sd")], idvar = c("scen", "n"), timevar = "method", direction = "wide"), digits = 2, row.names = FALSE)
cat("\n== Grafo di classe, n = 2000 ==\n")
print(mn[mn$n == 2000, c("scen", "method", "cl_tpr_sk", "cl_fpr_sk", "cl_orient", "cl_overorient", "cl_shd", "cl_exact", "err_Sig", "p_true_class")], digits = 2, row.names = FALSE)
cat("\n== Tutte le n: classe esatta, distanza, probabilita' della classe vera ==\n")
print(mn[c("scen", "n", "method", "cl_exact", "cl_shd", "cl_orient", "cl_overorient", "p_true_graph", "p_true_class")], digits = 2, row.names = FALSE)
cat("\n== modules, per tipo di grafo vero (n = 2000) ==\n")
m <- tab[tab$scen == "modules" & tab$n == 2000, ]
print(aggregate(m[c("cl_exact", "cl_shd", "cl_overorient", "p_true_class")], m[c("template", "method")], mean), digits = 2)
cat("\n== random: repliche con classe vera non banale ==\n")
r <- tab[tab$scen == "random" & tab$nontrivial == TRUE, ]
print(r[c("n", "rep", "method", "cl_exact", "cl_overorient", "p_true_graph", "p_true_class")], row.names = FALSE, digits = 2)
