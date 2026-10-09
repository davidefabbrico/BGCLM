# Fase 2.2 - A2: censimento dei grafi equivalenti SEMPLICI e realizzabili nella
# classe Gershgorin (a_i > 0), per i 10 grafi veri di gershgorin e di weak.
# Si esaminano tutti i supporti semplici con al piu' KMAX archi; per verificare
# A2 bastano quelli con al piu' 5 archi (il vero ne ha 5), gli altri sono
# riportati per completezza.
# Uso: Rscript inst/simulations/identifiability/02_census_simple_gershgorin.R
source("inst/simulations/identifiability/common_ident.R")
KMAX <- 7
sup <- simple_supports(KMAX)
cat("supporti semplici esaminati:", length(sup), "\n\n")
all <- list(); summ <- list()
for (scen in names(SEEDS)) for (r in 1:10) {
  tr <- true_model(scen, r); kv <- sum(tr$delta)
  cs <- census_simple(tr$B, tr$Sigma, tr$delta, sup)
  eq <- cs[!cs$vero, , drop = FALSE]
  g  <- eq[eq$gershgorin, , drop = FALSE]
  all[[length(all) + 1]] <- if (nrow(eq)) cbind(scenario = scen, rep = r, eq) else NULL
  summ[[length(summ) + 1]] <- data.frame(scenario = scen, rep = r, grafo_vero = edges_lab(tr$delta),
    archi_vero = kv, vero_trovato = sum(cs$vero),
    equivalenti_stabili = sum(eq$stabile), di_cui_gershgorin = nrow(g),
    gersh_meno_archi = sum(g$archi < kv), gersh_stessi_archi = sum(g$archi == kv),
    gersh_piu_archi = sum(g$archi > kv),
    stabili_meno_o_stessi = sum(eq$stabile & eq$archi <= kv))
}
summ <- do.call(rbind, summ); all <- do.call(rbind, all)
write.csv(summ, "inst/simulations/identifiability/output/02_census_summary.csv", row.names = FALSE)
if (!is.null(all)) write.csv(all, "inst/simulations/identifiability/output/02_census_equivalents.csv", row.names = FALSE)
print(summ, row.names = FALSE)
cat(sprintf("\nviolazioni di A2 (equivalenti Gershgorin con archi <= del vero): %d\n",
            sum(summ$gersh_meno_archi + summ$gersh_stessi_archi)))
cat(sprintf("equivalenti stabili (anche fuori dalla classe) con archi <= del vero: %d\n",
            sum(summ$stabili_meno_o_stessi)))
if (!is.null(all)) cat(sprintf("residuo massimo su Sigma degli equivalenti stabili: %.1e\n", max(all$residuo_Sigma, na.rm = TRUE)))
