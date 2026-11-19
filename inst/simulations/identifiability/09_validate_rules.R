# ============================================================
# Validazione numerica delle due regole su grafi diversi da quelli usati per
# formularle.
#
# Regola 1 (equivalenza esatta): due grafi semplici hanno lo stesso modello se e
#   solo se hanno lo stesso scheletro e uno si ottiene dall'altro riorientando
#   cricche-modulo (vedi 07_equivalence_classes.R). Su p = 6, 7, grafi casuali e
#   grafi con una cricca-modulo inserita. Per ogni grafo si confronta la classe
#   numerica (sistema lineare della fibra su tutti gli orientamenti dello
#   scheletro) con la chiusura della regola, e si cercano equivalenti con
#   scheletro diverso fra 300 grafi casuali con lo stesso numero di archi e fra
#   quelli che spostano un solo arco (in entrambi i versi).
# Regola 2 (archi deboli): con B = -(margini) + w U, un altro orientamento dello
#   stesso scheletro ha KL* = O(w^6) se ha lo stesso tipo (collider, fork,
#   catena in un verso o nell'altro) in ogni tripla non schermata a - k - c
#   (a, c non adiacenti), altrimenti O(w^4); KL* = 0 se e' equivalente (regola 1).
#   Si stima la pendenza di log KL* fra w = 0.05 e w = 0.1.
#
# Uso:  Rscript inst/simulations/identifiability/09_validate_rules.R regola1 <p>
#       Rscript inst/simulations/identifiability/09_validate_rules.R regola2
# Scrive inst/simulations/identifiability/output/09_regola1_p<p>.csv, 09_regola2.csv
# ============================================================
suppressMessages(library(BGCLM))
a <- commandArgs(TRUE); MODE <- a[1]
p <- if (MODE == "regola1") as.integer(a[2]) else 5
set.seed(900 + p + (MODE == "regola2"))
OFF <- which(diag(p) == 0); PAIRS <- which(upper.tri(diag(p)), arr.ind = TRUE); NP <- nrow(PAIRS)

lab <- function(d) { e <- which(d == 1, arr.ind = TRUE); e <- e[order(e[, 2], e[, 1]), , drop = FALSE]
  if (!nrow(e)) "vuoto" else paste(sprintf("%d->%d", e[, 2], e[, 1]), collapse = " ") }
skel <- function(d) (d | t(d)) * 1L
reach <- function(d) { R <- t(d) > 0; for (k in seq_len(p)) R <- R | (R[, k] %o% R[k, ]); R }
from_state <- function(s) { d <- matrix(0L, p, p)                     # stato per coppia: 0, 1 = r -> c, 2 = c -> r
  d[cbind(PAIRS[s == 1, 2], PAIRS[s == 1, 1])] <- 1L; d[cbind(PAIRS[s == 2, 1], PAIRS[s == 2, 2])] <- 1L; d }
orientations <- function(d) {                                          # tutti gli orientamenti dello scheletro
  q <- which(skel(d)[PAIRS] == 1)
  lapply(0:(2^length(q) - 1), function(o) { s <- integer(NP); s[q] <- ifelse(bitwAnd(o, 2^(seq_along(q) - 1)) > 0, 2L, 1L); from_state(s) })
}
cliques <- function(s) { out <- list()
  ext <- function(Q, cand) for (v in cand) { Q2 <- c(Q, v); if (length(Q2) >= 2) out[[length(out) + 1]] <<- Q2
    ext(Q2, cand[cand > v & s[v, cand] == 1]) }
  ext(integer(0), seq_len(p)); out }
module_moves <- function(x) {
  out <- list(); s <- skel(x); R <- reach(x) | diag(p) > 0
  for (T_ in cliques(s)) {
    pa_o <- lapply(T_, function(v) setdiff(which(x[v, ] == 1), T_)); ch_o <- lapply(T_, function(v) setdiff(which(x[, v] == 1), T_))
    if (!all(vapply(pa_o, setequal, TRUE, pa_o[[1]])) || !all(vapply(ch_o, setequal, TRUE, ch_o[[1]]))) next
    if (length(pa_o[[1]]) && length(ch_o[[1]]) && any(x[ch_o[[1]], pa_o[[1]]] == 0)) next
    oth <- setdiff(seq_len(p), c(T_, pa_o[[1]], ch_o[[1]]))
    trek <- function(k) any(vapply(T_, function(v) any(R[, k] & R[, v]), TRUE))
    if (!all(vapply(oth, function(k) all(s[k, T_] == 1) || !trek(k), TRUE))) next
    pr <- which(upper.tri(diag(length(T_))), arr.ind = TRUE)
    for (o in 0:(2^nrow(pr) - 1)) { x2 <- x; x2[T_, T_] <- 0L
      for (t in seq_len(nrow(pr))) { a_ <- T_[pr[t, 1]]; b_ <- T_[pr[t, 2]]
        if (bitwAnd(o, 2^(t - 1)) > 0) x2[b_, a_] <- 1L else x2[a_, b_] <- 1L }
      out[[length(out) + 1]] <- x2 }
  }
  out
}
closure <- function(x) {                                               # visita in ampiezza, grafi visti in una tabella hash
  seen <- new.env(hash = TRUE); assign(lab(x), TRUE, seen); queue <- list(x); i <- 1
  while (i <= length(queue)) { y <- queue[[i]]; i <- i + 1
    for (z in module_moves(y)) { l <- lab(z); if (!exists(l, seen, inherits = FALSE)) { assign(l, TRUE, seen); queue[[length(queue) + 1]] <- z } } }
  ls(seen) }
# firma della regola 2: tipo di ogni tripla non schermata a - k - c (C collider, F fork, K catena)
signature <- function(d) { s <- skel(d); out <- character(0)
  for (k in seq_len(p)) { nb <- which(s[k, ] == 1); if (length(nb) < 2) next
    for (u in combn(nb, 2, simplify = FALSE)) if (s[u[1], u[2]] == 0) {
      i1 <- d[k, u[1]] == 1; i2 <- d[k, u[2]] == 1
      out <- c(out, sprintf("%d-%d-%d:%s", u[1], k, u[2], if (i1 && i2) "C" else if (!i1 && !i2) "F" else "K")) } }
  sort(out) }
fiber <- function(S) { Si <- solve(S)
  sapply(seq_len(NP), function(q) { K <- matrix(0, p, p); K[PAIRS[q, 1], PAIRS[q, 2]] <- 1; K[PAIRS[q, 2], PAIRS[q, 1]] <- -1
    as.vector(K %*% Si)[OFF] }) }
consistent <- function(B, M, d2) { Z <- which(d2[OFF] == 0); f <- .lm.fit(M[Z, , drop = FALSE], -B[OFF][Z])
  sqrt(sum(f$residuals^2)) < 1e-8 * max(1, sqrt(sum(B^2))) }
# si escludono i grafi con cricche di 5 o piu' nodi: la loro classe ha 2^10 o piu' membri
# e la chiusura e' troppo lenta; il caso e' coperto dal censimento esaustivo con p = 5
max_clique <- function(d) max(c(1, lengths(cliques(skel(d)))))
rand_graph <- function(e) repeat { s <- integer(NP); s[sample(NP, e)] <- sample(1:2, e, TRUE); d <- from_state(s)
  if (max_clique(d) <= 4) return(d) }
plant_module <- function() {                                           # cricca-modulo T con genitori e figli comuni
  repeat {
    v <- sample(p); t <- sample(2:min(4, p - 2), 1); T_ <- v[1:t]; rest <- v[-(1:t)]
    np <- sample(0:min(2, length(rest)), 1); pa <- rest[seq_len(np)]; rest <- setdiff(rest, pa)
    nc <- sample(0:min(2, length(rest)), 1); ch <- rest[seq_len(nc)]; oth <- setdiff(rest, ch)
    d <- matrix(0L, p, p)
    for (u in combn(T_, 2, simplify = FALSE)) if (runif(1) < 0.5) d[u[2], u[1]] <- 1L else d[u[1], u[2]] <- 1L
    d[T_, pa] <- 1L; d[ch, T_] <- 1L; d[ch, pa] <- 1L
    for (u in combn(c(pa, ch, oth), 2, simplify = FALSE))               # archi casuali fuori da T
      if (skel(d)[u[1], u[2]] == 0 && runif(1) < 0.25) if (runif(1) < 0.5) d[u[2], u[1]] <- 1L else d[u[1], u[2]] <- 1L
    if (sum(d) <= 11 && max_clique(d) <= 4 && length(closure(d)) > 1) return(d)
  }
}
draw_B <- function(d) build_B(d * sample(c(-1, 1), p^2, TRUE) * runif(p^2, 0.3, 1.5), runif(p, 0.2, 1))

if (MODE == "regola1") {
  N_G <- 100; res <- list(); t0 <- proc.time()[3]
  for (g in seq_len(2 * N_G)) {
    tipo <- if (g <= N_G) "casuale" else "modulo inserito"
    message(sprintf("grafo %d (%s), %.0f s", g, tipo, proc.time()[3] - t0))
    d <- if (tipo == "casuale") rand_graph(sample(5:min(10, NP - 1), 1)) else plant_module()
    message("  ", lab(d))
    B <- draw_B(d); M <- fiber(lyapunov(B, diag(p))); e <- sum(d)
    ors <- orientations(d); ok <- vapply(ors, function(x) consistent(B, M, x), TRUE); num <- vapply(ors[ok], lab, "")
    rule <- closure(d)
    # scheletro diverso: 300 grafi casuali con e archi, e tutti gli spostamenti di un arco su una coppia vuota
    oth <- replicate(300, { repeat { x <- rand_graph(e); if (any(skel(x) != skel(d))) break }; x }, simplify = FALSE)
    st <- integer(NP); st[d[cbind(PAIRS[, 2], PAIRS[, 1])] == 1] <- 1L; st[d[cbind(PAIRS[, 1], PAIRS[, 2])] == 1] <- 2L
    for (i in which(st > 0)) for (j in which(st == 0)) for (o in 1:2) { s2 <- st; s2[i] <- 0L; s2[j] <- o; oth[[length(oth) + 1]] <- from_state(s2) }
    nsk <- sum(vapply(oth, function(x) consistent(B, M, x), TRUE))
    res[[g]] <- data.frame(p = p, tipo = tipo, archi = e, grafo = lab(d), orientamenti = length(ors),
                           classe_numerica = length(num), classe_regola = length(rule), coincidono = setequal(num, rule),
                           provati_scheletro_diverso = length(oth), equivalenti_scheletro_diverso = nsk,
                           stessa_firma_regola2 = all(vapply(ors[ok], function(x) identical(signature(x), signature(d)), TRUE)))
  }
  res <- do.call(rbind, res)
  write.csv(res, sprintf("inst/simulations/identifiability/output/09_regola1_p%d.csv", p), row.names = FALSE)
  cat(sprintf("Regola 1, p = %d: %d grafi, %.0f s\n", p, nrow(res), proc.time()[3] - t0))
  print(aggregate(cbind(grafi = 1, classe_non_banale = classe_numerica > 1, regola_giusta = coincidono,
                        provati_scheletro_diverso, equivalenti_scheletro_diverso, stessa_firma_regola2) ~ tipo, res, sum))
  cat("dimensione della classe numerica (grafi con modulo inserito):\n"); print(table(res$classe_numerica[res$tipo != "casuale"]))
  if (any(!res$coincidono)) { cat("\nDISCORDANZE:\n"); print(res[!res$coincidono, c("grafo", "classe_numerica", "classe_regola")], row.names = FALSE) }
}

if (MODE == "regola2") {
  kl <- function(S0, S) 0.5 * (sum(diag(solve(S, S0))) - nrow(S) + determinant(S)$modulus - determinant(S0)$modulus)
  # parametri fuori diagonale in unita' di w (il gradiente numerico di optim usa
  # passi assoluti); partenza che riproduce il primo ordine: l'arco i -> j del
  # vero, girato in j -> i, prende il valore B_ji d_j / d_i
  min_kl <- function(S0, Bt, d2, w) {
    idx <- which(d2 == 1); dd <- -diag(Bt)
    f <- function(th) { B <- matrix(0, p, p); B[idx] <- w * th[seq_along(idx)]; diag(B) <- -exp(th[-seq_along(idx)])
      if (!all(is.finite(B)) || max(Re(eigen(B, only.values = TRUE)$values)) >= -1e-10) return(1e6)
      v <- tryCatch(as.numeric(kl(S0, lyapunov(B, diag(p)))), error = function(e) 1e6); if (is.finite(v)) v else 1e6 }
    rc <- arrayInd(idx, c(p, p))
    first <- ifelse(Bt[idx] != 0, Bt[idx], Bt[cbind(rc[, 2], rc[, 1])] * dd[rc[, 2]] / dd[rc[, 1]])
    warm <- c(first / w, log(dd))
    min(sapply(1:3, function(s) { th0 <- if (s == 1) warm else warm + rnorm(length(warm), 0, 0.2)
      optim(th0, f, method = "BFGS", control = list(maxit = 20000, reltol = 1e-16, ndeps = rep(1e-6, length(th0))))$value }))
  }
  W <- c(0.05, 0.1); res <- list(); t0 <- proc.time()[3]; ng <- 0
  while (ng < 30) {
    d <- rand_graph(sample(4:6, 1)); ors <- orientations(d); sg <- signature(d); cl <- closure(d)
    cat_ <- vapply(ors, function(x) if (lab(x) %in% cl) "esatto (regola 1)" else if (identical(signature(x), sg))
      "stesso tipo nelle triple (w^6)" else "tipo diverso (w^4)", "")
    if (!any(cat_ == "stesso tipo nelle triple (w^6)") && runif(1) < 0.7) next   # privilegia i grafi con casi w^6
    ng <- ng + 1
    pick <- unlist(lapply(split(seq_along(ors), cat_), function(ix) ix[ix != match(lab(d), vapply(ors, lab, ""))][seq_len(min(3, length(ix)))]))
    pick <- pick[!is.na(pick)]
    U <- d * sample(c(-1, 1), p^2, TRUE) * runif(p^2, 0.5, 1.5); mg <- runif(p, 0.3, 1)
    for (h in pick) for (mt in c("diversi", "uguali")) {               # margini casuali o tutti 0.5
      m_ <- if (mt == "uguali") rep(0.5, p) else mg
      k <- vapply(W, function(w) { Bt <- build_B(w * U, m_); max(min_kl(lyapunov(Bt, diag(p)), Bt, ors[[h]], w), 0) }, 1)
      res[[length(res) + 1]] <- data.frame(grafo = lab(d), alternativa = lab(ors[[h]]), previsione = cat_[h], margini = mt,
                                           kl_005 = k[1], kl_01 = k[2], pendenza = log(k[2] / k[1]) / log(2))
    }
  }
  res <- do.call(rbind, res)
  write.csv(res, "inst/simulations/identifiability/output/09_regola2.csv", row.names = FALSE)
  cat(sprintf("Regola 2: %d grafi, %d alternative, %.0f s\n", ng, nrow(res), proc.time()[3] - t0))
  res$osservato <- ifelse(res$kl_01 < 1e-11, "KL = 0", ifelse(res$pendenza > 5, "pendenza > 5 (w^6)", "pendenza <= 5 (w^4)"))
  for (mt in c("uguali", "diversi")) { z <- res[res$margini == mt, ]
    cat(sprintf("\nmargini %s:\n", mt)); print(table(previsione = z$previsione, osservato = z$osservato))
    print(tapply(z$pendenza, z$previsione, function(v) round(quantile(v, c(0, 0.1, 0.5, 0.9, 1), na.rm = TRUE), 2))) }
}
