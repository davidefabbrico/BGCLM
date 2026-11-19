# Fase 2.3 - Genericita': per ogni grafo vero si ridisegnano NDRAW volte i
# valori degli archi (ampiezze uniformi nell'intervallo dello scenario, segni
# casuali, diagonale di Gershgorin con a_i = 0.5) e si ripete il censimento dei
# grafi semplici. Se gli stessi equivalenti compaiono per tutti i valori, si
# tratta di una proprieta' del grafo (inclusione di modelli), non di valori
# particolari.
# Uso: Rscript inst/simulations/identifiability/03_genericity.R
source("inst/simulations/identifiability/common_ident.R")
NDRAW <- 20
RANGE <- list(gershgorin = c(0.4, 0.8), weak = c(0.10, 0.25), random = c(0.4, 0.8))
sup <- simple_supports(5 + 2)
res <- list(); sets_ok <- list()
for (scen in names(SEEDS)) for (r in 1:10) {
  tr <- true_model(scen, r); S <- tr$delta == 1; kv <- sum(S)
  sets <- list()
  for (dd in seq_len(NDRAW)) {
    set.seed(5e4 + 1000 * match(scen, names(SEEDS)) + 50 * r + dd)
    Bo <- matrix(0, p, p)
    Bo[S] <- sample(c(-1, 1), kv, TRUE) * runif(kv, RANGE[[scen]][1], RANGE[[scen]][2])
    B <- build_B(Bo, rep(0.5, p)); Sg <- lyapunov(B, diag(p))
    cs <- census_simple(B, Sg, tr$delta, sup)
    eq <- cs[!cs$vero, , drop = FALSE]
    sets[[dd]] <- sort(eq$grafo[eq$stabile])
    res[[length(res) + 1]] <- data.frame(scenario = scen, rep = r, estrazione = dd,
      equivalenti_stabili = sum(eq$stabile), gershgorin = sum(eq$gershgorin),
      viol_A2 = sum(eq$gershgorin & eq$archi <= kv), viol_stabili = sum(eq$stabile & eq$archi <= kv))
  }
  sets_ok[[length(sets_ok) + 1]] <- data.frame(scenario = scen, rep = r,
    stesso_insieme_in_tutte = all(sapply(sets, identical, sets[[1]])))
}
res <- do.call(rbind, res); sets_ok <- do.call(rbind, sets_ok)
write.csv(res, "inst/simulations/identifiability/output/03_genericity.csv", row.names = FALSE)
ag <- aggregate(cbind(equivalenti_stabili, gershgorin, viol_A2, viol_stabili) ~ scenario + rep, res,
                function(z) sprintf("%g-%g", min(z), max(z)))
ag <- merge(ag, sets_ok)
names(ag)[3:6] <- paste0(names(ag)[3:6], " (min-max)")
print(ag[order(ag$scenario, ag$rep), ], row.names = FALSE)
cat(sprintf("\nestrazioni totali: %d;  violazioni di A2: %d;  equivalenti stabili con archi <= del vero: %d\n",
            nrow(res), sum(res$viol_A2), sum(res$viol_stabili)))
cat(sprintf("grafi per cui l'insieme degli equivalenti e' lo stesso in tutte le %d estrazioni: %d su %d\n",
            NDRAW, sum(sets_ok$stesso_insieme_in_tutte), nrow(sets_ok)))
