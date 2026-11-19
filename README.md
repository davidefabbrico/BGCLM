# BGCLM — pacchetto R

Bayesian Graphical Continuous Lyapunov Models. Il paper è in `Desktop/Davide/Ricerca/BGCLM/paper`.

## Codice
| file | contenuto |
|---|---|
| `R/gclm_mcmc.R` | `gclm_mcmc()` (unica interfaccia ai campionatori), `restrict_simple()`, `select_graph()` (regola per coppia) |
| `R/class_graph.R` | `clique_modules()`, `model_class()`, `class_graph()`, `class_metrics()`: classe di equivalenza a posteriori |
| `R/plot.R`, `R/utils.R` | grafici e funzioni di servizio |
| `src/gclm_mcmc_pt.cpp` | campionatore del paper: modello A1 (grafi semplici), augmentazione con eps, parallel tempering |
| `src/gclm_mcmc_simple.cpp` | modello A1, catena singola (`simple = TRUE` senza `pt_betas`) |
| `src/gclm_mcmc.cpp` | modello libero (`simple = FALSE`), non usato nel paper, serve ai controlli |
| `inst/tinytest/` | test |

## Simulazioni (`inst/simulations/`)
| script | cosa produce nel paper |
|---|---|
| `compare_gclm.R`, `gclm_vh_sim.R`, `plot_results.R` | Sez. 4: AUC, TPR/FPR, errore su Sigma, copertura; figure in `figures/` |
| `simple_vs_free.R` | recupero del grafo esatto e orientamento (Tabella "graph") |
| `class_comparison.R`, `class_collect.R` | confronto sul grafo di classe; output in `results/data/class_*` |
| `equal_diag.R` | variante con diagonale uguale (non nel paper) |
| `sampler_checks/` | Sez. 4.4 (mixing): taratura di eps, scala delle temperature, catene multiple, controlli sul prior |
| `identifiability/01`–`06` | Sez. 4.3: rango dato il supporto, censimento degli equivalenti, genericità, dimensione, quasi-equivalenza |
| `identifiability/07`–`12` | Sez. 2.4–2.6 e appendice: classi di equivalenza esaustive (p <= 5), informazione sull'orientamento, validazione delle regole, controlli delle dimostrazioni, censimento esatto (Python), riduzione della congettura |

Output: `inst/simulations/*.rds`, `inst/simulations/results/data/`, `inst/simulations/identifiability/output/`.

## `_archivio/`
Copia del 9 ottobre degli script e dei risultati (`results_snapshot_2026-10-09`, superata) e la vecchia `.Rhistory`. Esclusa da git e dalla build.
