# Risultati di simulazione — BGCLM vs gclm

Confronto fra il campionatore bayesiano **BGCLM** (pacchetto in sviluppo) e la
stima LASSO di Varando et al. (**gclm**, CRAN 0.0.1) per i Graphical Continuous
Lyapunov Models. Tutti i numeri sono riprodotti dai file in `data/`.

> **Versione dell'8 ottobre 2026 (A1 + parallel tempering).** Sostituisce del
> tutto la versione precedente (BGCLM-s, 5 scenari, senza tempering), salvata in
> `backup/2026-10-08_pre_A1A2/results_old/`. Non mescolare i numeri.
> Aggiornamento: per BGCLM (A1) il grafo è selezionato con la regola per coppia
> (§2); cambiano solo TPR/FPR, grafo esatto e orientamento di A1, il resto è
> identico (stesse catene, semi fissi).

---

## 1. Modello e metodi

```
B Σ + Σ Bᵀ + C = 0          y_l | B ~ iid N_p(0, Σ(B, C)),   C = I
```

Convenzione: `B[j,i] ≠ 0 ⟺ arco i → j`. Diagonale di Gershgorin
`B_ii = −(Σ_{j≠i} |B_ij| + a_i)`, `a_i ~ Gamma(2, 2)`.

| metodo | descrizione |
|---|---|
| **BGCLM (A1)** | modello sotto l'Assunzione 1 (grafo semplice, niente coppie i ⇄ j). Prior sulle coppie: nessun arco con prob. 1−ρ, ciascuna direzione ρ/2, ρ ~ Beta(1,1). Augmentazione computazionale ε = 0.1 (la catena può passare per coppie reciproche; si tengono solo i draw semplici, la posteriori ristretta è esatta). `gclm_mcmc(..., simple = TRUE)` |
| **BGCLM** | modello libero, spike-and-slab Bernoulli(π) su ogni arco, π ~ Beta(1,1). `gclm_mcmc(..., simple = FALSE)` |
| **gclm (C = I)** | path LASSO con C fissata, selezione BIC |
| **gclm (C diag.)** | path LASSO con C diagonale stimata (λ_C = 0.01, non tarato) |

Entrambe le versioni di BGCLM usano **parallel tempering** con K = 10 livelli,
β geometriche da 1 a 0.01 (`pt_betas = 0.01^(seq(0, 1, length.out = 10))`),
10.000 iterazioni, 2.000 burn-in, thin 2 → 4.000 draw a β = 1. BGCLM (A1) ne
tiene il 39–53%.

## 2. Disegno

| | |
|---|---|
| nodi | p = 5 (20 archi diretti possibili) |
| archi veri | 5 (ciclo 1→2→3→1 con segni +,−,+ più 2 archi casuali) |
| n | {100, 500, 2000} |
| repliche | 10 per cella, semi fissi |
| path LASSO | 30 λ log-spaziati da 0.6 a 0.004, warm start; BIC con df = archi + p (archi + 2p − 1 con C stimata) |

| scenario | costruzione del vero |
|---|---|
| `gershgorin` | \|B_ji\| ∈ [0.4, 0.8], a_i = 0.5 |
| `weak` | come sopra con \|B_ji\| ∈ [0.10, 0.25] |
| `two_cycle` | come `gershgorin` più l'inverso di uno dei due archi casuali: il vero ha **una coppia reciproca** (6 archi) e viola A1 |

Righe: 360 (`data/raw.csv` = 3 scenari × 3 n × 10 repliche × 4 metodi).

**Regola di selezione** (`select_graph()` nel pacchetto). BGCLM libero: freccia
i → j se MPPI_ji > 0.5. BGCLM (A1): gli stati di una coppia si escludono a
vicenda, quindi MPPI_ij + MPPI_ji è la probabilità che la coppia abbia un arco; la
coppia è selezionata se questa supera 0.5 e l'arco va nel verso con MPPI maggiore.
La regola per freccia sottoselezionerebbe quando il verso è incerto (due MPPI di
0.4 → nessun arco). gclm: archi non nulli del modello scelto dal BIC.

**Metriche.** `auc` (punteggio = MPPI per BGCLM, λ di ingresso nel path per
gclm; non dipende da una soglia); `tpr`/`fpr` del grafo selezionato sulle 20
frecce; `tpr_sk`/`fpr_sk` gli stessi tassi sulle 10 coppie dello scheletro (arco
non orientato); `err_Sig` = ‖Σ̂ − Σ‖_F / ‖Σ‖_F; `cov_e`/`cov_0` = copertura degli
intervalli al 95% su archi veri / non archi (solo BGCLM); `tpr_2ciclo` = frazione
dei due archi della coppia reciproca selezionati (solo `two_cycle`).

---

## 3. AUC (media, sd fra parentesi)

| scenario | n | BGCLM (A1) | BGCLM | gclm (C=I) | gclm (C diag.) |
|---|---|---|---|---|---|
| gershgorin | 100 | **0.90** (0.08) | 0.89 (0.07) | 0.66 (0.10) | 0.72 (0.10) |
| gershgorin | 500 | **0.98** (0.04) | **0.98** (0.04) | 0.76 (0.09) | 0.83 (0.08) |
| gershgorin | 2000 | **0.98** (0.04) | **0.98** (0.04) | 0.77 (0.12) | 0.76 (0.12) |
| weak | 100 | 0.69 (0.13) | **0.70** (0.12) | 0.62 (0.15) | 0.64 (0.12) |
| weak | 500 | **0.81** (0.06) | **0.81** (0.05) | 0.72 (0.10) | 0.72 (0.13) |
| weak | 2000 | **0.91** (0.09) | 0.89 (0.08) | 0.75 (0.12) | 0.78 (0.10) |

Nei due scenari corretti le due versioni di BGCLM sono quasi uguali e sopra gclm
in tutte le celle; in `weak` a n = 100 il margine è entro una sd.

## 4. Grafo selezionato ed errore su Σ, n = 2000

| scenario | | BGCLM (A1) | BGCLM | gclm (C=I) | gclm (C diag.) |
|---|---|---|---|---|---|
| gershgorin | TPR / FPR frecce | 0.88 / **0.04** | **0.94** / 0.06 | 0.70 / 0.19 | 0.68 / 0.25 |
| gershgorin | TPR / FPR scheletro | 1.00 / 0.00 | 1.00 / 0.00 | 1.00 / 0.22 | 1.00 / 0.32 |
| gershgorin | err Σ | 0.046 | 0.046 | 0.053 | 0.056 |
| weak | TPR / FPR frecce | 0.64 / 0.12 | **0.74** / 0.17 | 0.66 / **0.11** | 0.72 / 0.13 |
| weak | TPR / FPR scheletro | **0.98** / 0.02 | 0.88 / **0.00** | 0.96 / 0.02 | 0.98 / 0.04 |
| weak | err Σ | 0.049 | 0.048 | 0.050 | 0.049 |

La TPR per freccia favorisce il modello che può selezionare entrambe le direzioni:
con il verso incerto il modello libero tiene le due frecce e guadagna un vero
positivo (0.20 di TPR su 5 archi veri) al prezzo di un falso positivo (0.07 di FPR
su 15 non archi); A1 deve sceglierne una, e se sbaglia paga un falso negativo e un
falso positivo. Sullo scheletro A1 trova tutte le coppie vere in `gershgorin` e
0.98 in `weak` (libero 0.88). In `weak` la distanza fra 0.98 sullo scheletro e
0.64 sulle frecce è l'orientamento. Nel paper lo scheletro è riportato solo per
le due versioni di BGCLM; quello di gclm è in `data/summary_mean.csv`.

Effetto della regola per coppia su A1 rispetto alla regola per freccia (stesse
catene): a n = 2000 identico in `gershgorin`, in `weak` TPR/FPR da 0.58/0.10 a
0.64/0.12. A n piccolo seleziona di più: `gershgorin` n = 100 da 0.42/0.03 a
0.64/0.17, n = 500 da 0.92/0.03 a 0.94/0.07.

## 5. Copertura al 95%

| versione | scenario | archi veri n=100 | 500 | 2000 | non archi |
|---|---|---|---|---|---|
| BGCLM (A1) | gershgorin | 0.74 | 1.00 | 0.94 | 1.00 |
| BGCLM (A1) | weak | 0.60 | 0.74 | 0.94 | 1.00 |
| BGCLM | gershgorin | 0.74 | 1.00 | 0.94 | 1.00 |
| BGCLM | weak | 0.60 | 0.80 | 0.98 | 1.00 |

50 intervalli per cella (5 archi × 10 repliche): errore Monte Carlo ≈ 0.03
vicino a 0.95. Sui non archi l'intervallo collassa spesso in [0, 0].

## 6. Recupero del grafo esatto — `data/simple_vs_free.csv`

Esperimento separato (`scripts/simple_vs_free.R`), 10 repliche per cella.
Valori libero (regola per freccia) / A1 (regola per coppia).

| scenario | n | P(grafo vero \| y) | grafo esatto selezionato | orientamento | FPR |
|---|---|---|---|---|---|
| gershgorin | 500 | 0.04 / **0.09** | 0.1 / **0.3** | 0.62 / **0.94** | 0.14 / **0.07** |
| gershgorin | 2000 | 0.15 / **0.30** | 0.6 / **0.7** | 0.82 / **0.88** | 0.06 / **0.04** |
| gershgorin | 10000 | 0.32 / **0.49** | 0.6 / 0.6 | 0.88 / 0.88 | **0.05** / 0.06 |
| weak | 500 | 0.00 / 0.00 | 0.0 / 0.0 | 0.14 / **0.28** | **0.09** / 0.10 |
| weak | 2000 | 0.01 / **0.02** | 0.0 / **0.1** | 0.38 / **0.64** | 0.17 / **0.12** |
| weak | 10000 | 0.09 / **0.14** | 0.2 / **0.5** | 0.70 / **0.88** | 0.08 / **0.05** |

"Grafo esatto" = frazione di repliche in cui il grafo con MPPI > 0.5 coincide
col vero; "orientamento" = frazione di archi veri selezionati senza il loro
inverso. La probabilità del grafo vero cresce con n in entrambe le versioni.

Nella replica delle figure (`gershgorin`, rep 1) a n = 500 e 2000 A1 seleziona lo
scheletro vero con il 3-ciclo invertito (1→3→2→1), MPPI 0.55–0.76 contro
0.24–0.45. Le catene con tempering concordano: è una proprietà della posteriori
per quei dati, non un problema di mixing.

## 7. Vero con una coppia reciproca — `two_cycle`, n = 2000

| | BGCLM (A1) | BGCLM | gclm (C=I) | gclm (C diag.) |
|---|---|---|---|---|
| AUC | 0.83 | **0.94** | 0.68 | 0.68 |
| TPR / FPR | 0.72 / 0.16 | **0.92** / **0.15** | 0.58 / 0.29 | 0.60 / 0.30 |
| TPR / FPR scheletro | 0.98 / 0.34 | 0.98 / **0.10** | 0.98 / 0.50 | 0.96 / 0.40 |
| TPR sulla coppia reciproca | 0.45 | **0.80** | 0.45 | 0.50 |
| err Σ | 0.055 | 0.054 | 0.061 | 0.063 |
| copertura archi veri | 0.72 | 0.88 | – | – |

Quando A1 è falsa A1 può tenere al più un arco della coppia: perde AUC e
copertura, e sullo scheletro aggiunge coppie false (FPR 0.34 contro 0.10), a
compensare l'arco che non può selezionare. Σ resta ben ricostruita.

## 8. Verifica numerica delle assunzioni — `data/ident_*.csv`

Sui 10 grafi veri di `gershgorin` e i 10 di `weak` (`scripts/identifiability/`):

- **A1 / rango:** dato il supporto vero, il sistema lineare in K (B' = B + KΣ⁻¹,
  K antisimmetrica) ha rango pieno p(p−1)/2 = 10 in 20 casi su 20.
- **A2 / censimento:** tutti i 41.384 supporti semplici con ≤ 7 archi. Nessuna
  altra deriva semplice con ≤ 5 archi riproduce Σ, dentro o fuori la classe di
  Gershgorin: **A2 vale in 20 casi su 20**.

| scenario | grafi con equivalenti | equivalenti per grafo | di cui Gershgorin | con ≤ 5 archi | dimensione del modello |
|---|---|---|---|---|---|
| gershgorin | 7 su 10 | 7–62 | 1–17 | 0 | 11–12 |
| weak | 4 su 10 | 62 | 60–62 | 0 | 11 |

Gli equivalenti hanno 6–7 archi e dimensione 11–12 contro 10 del vero: modelli
più grandi che contengono il vero. Ridisegnando 20 volte i valori di B per ogni
grafo si ottiene sempre lo stesso insieme di equivalenti (0 violazioni su 400).

## 9. Mixing — `data/pt_multichain.csv`

Per ogni replica di `gershgorin` a n = 2000 e 10000, 4 catene con tempering (2 dal
grafo vuoto, 2 dal vero):

| modello | n | max diff. MPPI vuoto vs vero | R̂ (E, loglik) | tasso di scambio minimo |
|---|---|---|---|---|
| A1 | 2000 | 0.08 | ≤ 1.01 | 0.25 |
| A1 | 10000 | 0.07 | ≤ 1.00 | 0.25 |
| libero | 2000 | 0.05 | ≤ 1.00 | 0.24 |
| libero | 10000 | 0.16 | ≤ 1.00 | 0.24 |

Le probabilità a posteriori del grafo vero dalle due partenze coincidono.

## 10. Costo

| metodo | secondi per fit (min – mediana – max) |
|---|---|
| BGCLM (A1), PT K = 10 | 7.2 – 14.6 – 23.2 |
| BGCLM, PT K = 10 | 6.2 – 14.4 – 24.3 |
| gclm (C = I) | 0.001 – 0.004 – 0.014 |
| gclm (C diag.) | 0.004 – 0.039 – 0.986 |

Il tempering moltiplica per K il costo per iterazione (una soluzione di
Lyapunov O(p³) per proposta).

---

## 11. Avvertenze

1. L'identificabilità vale **sotto A1 + A2**. A2 non è verificabile dai dati; i
   numeri del §8 dicono solo che vale nei grafi simulati (p = 5, C nota).
2. Identificabile non vuol dire concentrazione rapida: gli equivalenti sono
   modelli con 1–2 parametri in più, penalizzati lentamente. In `weak` il grafo
   esatto si recupera solo a n = 10000, e non sempre.
3. Se il vero ha coppie reciproche A1 è sbagliata e peggiora (§7).
4. Errore su Σ piccolo ovunque: non confonderlo con il recupero di B.
5. 10 repliche, sd grandi (`data/summary_sd.csv`); TPR/FPR di gclm dipendono
   dalla regola BIC; "gclm C diag." non è tarato.
6. Tutto a p = 5.

## 12. Contenuto della cartella

```
README.md
data/
  raw.csv                      360 righe, una per (scenario, n, replica, metodo)
  summary_mean.csv / summary_sd.csv   medie e sd per cella
  compare_gclm_results.rds     come raw.csv, formato R
  compare_gclm_results_pip.rds MPPI di BGCLM (A1) e libero per ogni replica, con il grafo vero
  simple_vs_free.csv           §6, 120 righe
  simple_vs_free_pip.rds       MPPI corrispondenti
  pt_multichain.csv            §9, 40 righe
  ident_01_rank_given_support.csv    §8, rango
  ident_02_census_summary.csv        §8, censimento per grafo
  ident_02_census_equivalents.csv    §8, elenco degli equivalenti
  ident_03_genericity.csv            §8, ridisegno dei valori
  ident_04_model_dimension.csv       §8, dimensione dei modelli
figures/
  graph_correct_spec.png       grafo a posteriori BGCLM (A1), gershgorin e weak, n = 500
  graphs_by_n.png              BGCLM (A1), gershgorin, n = 100 / 500 / 2000
  auc_coverage.png             AUC e copertura vs n, 3 scenari
scripts/
  compare_gclm.R               simulazione principale (§3–5, §7, §10)
  simple_vs_free.R             §6
  pt_multichain.R              §9
  plot_results.R               figure
  identifiability/             §8
```

**Figure del grafo.** Grigio e spessore = MPPI (nero = 1, non disegnato sotto
0.10). Alone arancione = arco vero: arco scuro senza alone = falso positivo,
alone senza arco = falso negativo. Tratteggio = coefficiente negativo.

**Riproduzione** (dalla radice del pacchetto; richiede `BGCLM` installato e
`gclm` 0.0.1):

```sh
Rscript inst/simulations/compare_gclm.R                 # 3 scenari; oppure uno per processo: ... compare_gclm.R weak
Rscript inst/simulations/simple_vs_free.R                # oppure: ... simple_vs_free.R gershgorin 2000
Rscript inst/simulations/sampler_checks/pt_multichain.R A1 2000      # A1|libero, 2000|10000
Rscript inst/simulations/identifiability/01_rank_given_support.R  # poi 02, 03, 04
Rscript inst/simulations/plot_results.R
```

Con un argomento `compare_gclm.R` scrive `compare_gclm_results_<scenario>.rds`;
`data/compare_gclm_results.rds` è l'unione dei tre. I semi sono fissi.
