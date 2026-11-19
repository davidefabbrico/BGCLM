# ============================================================
# Classi di equivalenza di modello fra grafi semplici, esaustive su p nodi.
#
# Per ogni grafo semplice G (un rappresentante per classe di isomorfismo) si
# estrae un drift generico B su G, Sigma = Sigma(B, I), e per ogni grafo
# semplice G' con lo stesso numero di archi si verifica se esiste K antisimmetrica
# con B + K Sigma^-1 supportata su G' (sistema lineare in K, Prop. "Fibre").
# Due grafi semplici hanno dimensione e + p, e le chiusure sono irriducibili,
# quindi con probabilita' uno questo accade se e solo se le chiusure di Zariski
# coincidono. Archi diversi -> dimensione diversa -> mai equivalenti.
#
# Ogni classe viene confrontata con la chiusura di G sotto due mosse note:
#   (a) riorientare a piacere gli archi di una componente connessa con
#       scheletro completo (arco isolato, triangolo, K4: modello saturo);
#   (b) invertire un arco super-covered di una componente aciclica, con la
#       definizione 3.6 di Amendola et al. (2025) come parafrasata nel paper.
# "spiegata" = la classe numerica coincide con la chiusura sotto (a) e (b).
# Si verifica anche la mossa (c), che le unifica (vedi module_moves).
# Per le classi non banali si stima anche la frazione di drift in B (Gershgorin)
# per cui un equivalente cade in B, cioe' A2 fallisce.
#
# Uso:  Rscript inst/simulations/identifiability/07_equivalence_classes.R <p>
# Scrive inst/simulations/identifiability/output/07_classes_p<p>.csv
# ============================================================
suppressMessages(library(BGCLM))
p <- as.integer(commandArgs(TRUE)[1]); set.seed(700 + p)
OFF <- which(diag(p) == 0); PAIRS <- which(upper.tri(diag(p)), arr.ind = TRUE); NP <- nrow(PAIRS)
N_GERSH <- if (p <= 4) 30 else 10

# tutti i grafi semplici: stato per coppia 0 = nessuno, 1 = r -> c, 2 = c -> r
st <- as.matrix(expand.grid(rep(list(0:2), NP))); dimnames(st) <- NULL
to_mat <- function(s) { d <- matrix(0L, p, p)
  d[cbind(PAIRS[s == 1, 2], PAIRS[s == 1, 1])] <- 1L; d[cbind(PAIRS[s == 2, 1], PAIRS[s == 2, 2])] <- 1L; d }
SUP <- t(apply(st, 1, function(s) to_mat(s)[OFF] == 1))       # supporti sulle posizioni OFF
NE <- rowSums(st > 0)

# forma canonica per isomorfismo: minimo del codice binario su tutte le permutazioni
perms <- function(v) if (length(v) == 1) matrix(v) else do.call(rbind, lapply(seq_along(v), function(i) cbind(v[i], perms(v[-i]))))
PM <- perms(seq_len(p)); w <- 2^(seq_along(OFF) - 1)
code <- sapply(seq_len(nrow(PM)), function(k) {
  idx <- matrix(seq_len(p^2), p)[PM[k, ], PM[k, ]][OFF]                # posizione OFF dopo la permutazione
  as.vector((SUP[, match(idx, OFF)] * 1) %*% w) })
canon <- apply(code, 1, min); reps <- which(!duplicated(canon))

lab <- function(d) { e <- which(d == 1, arr.ind = TRUE); e <- e[order(e[, 2], e[, 1]), , drop = FALSE]
  if (!nrow(e)) "vuoto" else paste(sprintf("%d->%d", e[, 2], e[, 1]), collapse = " ") }
reach <- function(d) { A <- t(d) > 0; R <- A; for (k in seq_len(p)) R <- R | (R[, k] %o% R[k, ]); R }   # R[u,v]: u ~> v
acyclic <- function(d) !any(diag(reach(d)))
skel <- function(d) (d | t(d)) * 1L
comps <- function(d) { s <- skel(d) | diag(p) > 0; R <- s; for (k in seq_len(p)) R <- R | (R[, k] %o% R[k, ]); unique(lapply(seq_len(p), function(i) which(R[i, ]))) }

fiber <- function(Sigma) { Si <- solve(Sigma)
  sapply(seq_len(NP), function(q) { K <- matrix(0, p, p); K[PAIRS[q, 1], PAIRS[q, 2]] <- 1; K[PAIRS[q, 2], PAIRS[q, 1]] <- -1
    as.vector(K %*% Si)[OFF] }) }
kfull <- function(k) { K <- matrix(0, p, p); K[PAIRS] <- k; K - t(K) }
draw_B <- function(d) { Bo <- d * sample(c(-1, 1), p^2, TRUE) * runif(p^2, 0.3, 1.5); build_B(Bo, runif(p, 0.2, 1)) }
dominant <- function(B) all(-diag(B) - (rowSums(abs(B)) - abs(diag(B))) > 1e-10)
# B' = B + K Sigma^-1 supportato su sup, se esiste; NULL altrimenti
competitor <- function(B, M, Si, sup) {
  Z <- which(!sup); b <- B[OFF]
  f <- .lm.fit(M[Z, , drop = FALSE], -b[Z])
  if (sqrt(sum(f$residuals^2)) > 1e-8 * max(1, sqrt(sum(b^2)))) return(NULL)
  B + kfull(f$coefficients) %*% Si
}

# arco i -> j super-covered (parafrasi della Def. 3.6 nel paper)
supercovered <- function(d, i, j) {
  pa <- function(v) which(d[v, ] == 1); ch <- function(v) which(d[, v] == 1)
  if (!setequal(pa(j), c(pa(i), i)) || !setequal(ch(i), c(ch(j), j))) return(FALSE)
  if (!all(d[ch(i), pa(j), drop = FALSE] == 1 | outer(ch(i), pa(j), `==`))) return(FALSE)
  R <- reach(d) | diag(p) > 0; anc <- function(v) R[, v]
  trek <- function(a, b) any(anc(a) & anc(b)); s <- skel(d)
  all(sapply(setdiff(seq_len(p), c(i, j)), function(k) (s[k, i] && s[k, j]) || (!trek(k, i) && !trek(k, j))))
}
moves <- function(d) {
  out <- list()
  for (cc in comps(d)) {
    if (length(cc) < 2) next
    s <- skel(d)[cc, cc]
    if (all(s[upper.tri(s)] == 1)) {                                   # (a) componente con scheletro completo
      pr <- which(upper.tri(s), arr.ind = TRUE)
      for (o in 0:(2^nrow(pr) - 1)) { d2 <- d; d2[cc, cc] <- 0L
        for (t in seq_len(nrow(pr))) { a <- cc[pr[t, 1]]; b <- cc[pr[t, 2]]
          if (bitwAnd(o, 2^(t - 1)) > 0) d2[b, a] <- 1L else d2[a, b] <- 1L }
        out[[length(out) + 1]] <- d2 }
    } else {
      sub <- matrix(0L, p, p); sub[cc, cc] <- d[cc, cc]
      if (!acyclic(sub)) next
      for (e in which(d[cc, cc] == 1)) {                               # (b) inversione super-covered
        j <- cc[(e - 1) %% length(cc) + 1]; i <- cc[(e - 1) %/% length(cc) + 1]
        if (supercovered(d, i, j)) { d2 <- d; d2[j, i] <- 0L; d2[i, j] <- 1L; out[[length(out) + 1]] <- d2 }
      }
    }
  }
  out
}
closure_moves <- function(d, mv = moves) {
  seen <- lab(d); queue <- list(d)
  while (length(queue)) { x <- queue[[1]]; queue <- queue[-1]
    for (y in mv(x)) { l <- lab(y); if (!l %in% seen) { seen <- c(seen, l); queue[[length(queue) + 1]] <- y } } }
  seen
}
# (c) congettura: una cricca T i cui nodi hanno gli stessi genitori e gli stessi
# figli fuori da T ("cricca-modulo"), con ogni genitore esterno che punta a ogni
# figlio esterno e ogni altro nodo adiacente a tutta T o senza trek con T, si
# riorienta liberamente (anche in cicli). |T| = 2 e' l'arco super-covered,
# T componente completa e' il caso (a).
module_moves <- function(d) {
  out <- list(); s <- skel(d); R <- reach(d) | diag(p) > 0
  for (m in 3:(2^p - 1)) {
    T_ <- which(bitwAnd(m, 2^(seq_len(p) - 1)) > 0); if (length(T_) < 2) next
    if (any(s[T_, T_][upper.tri(diag(length(T_)))] == 0)) next
    pa_o <- lapply(T_, function(v) setdiff(which(d[v, ] == 1), T_)); ch_o <- lapply(T_, function(v) setdiff(which(d[, v] == 1), T_))
    if (!all(vapply(pa_o, setequal, TRUE, pa_o[[1]])) || !all(vapply(ch_o, setequal, TRUE, ch_o[[1]]))) next
    if (length(pa_o[[1]]) && length(ch_o[[1]]) && any(d[ch_o[[1]], pa_o[[1]]] == 0)) next
    oth <- setdiff(seq_len(p), c(T_, pa_o[[1]], ch_o[[1]]))
    trek <- function(k) any(vapply(T_, function(v) any(R[, k] & R[, v]), TRUE))
    if (!all(vapply(oth, function(k) all(s[k, T_] == 1) || !trek(k), TRUE))) next
    pr <- which(upper.tri(diag(length(T_))), arr.ind = TRUE)
    for (o in 0:(2^nrow(pr) - 1)) { d2 <- d; d2[T_, T_] <- 0L
      for (t in seq_len(nrow(pr))) { a <- T_[pr[t, 1]]; b <- T_[pr[t, 2]]
        if (bitwAnd(o, 2^(t - 1)) > 0) d2[b, a] <- 1L else d2[a, b] <- 1L }
      out[[length(out) + 1]] <- d2 }
  }
  out
}

t0 <- proc.time()[3]; res <- vector("list", length(reps))
for (r in seq_along(reps)) {
  g <- reps[r]; d <- to_mat(st[g, ]); B <- draw_B(d); S <- lyapunov(B, diag(p)); Si <- solve(S); M <- fiber(S)
  cand <- which(NE == NE[g]); eq <- cand[vapply(cand, function(h) !is.null(competitor(B, M, Si, SUP[h, ])), TRUE)]
  eqlab <- vapply(eq, function(h) lab(to_mat(st[h, ])), "")
  cl <- if (length(eq) > 1) closure_moves(d) else lab(d)
  clm <- closure_moves(d, module_moves)
  sameskel <- all(vapply(eq, function(h) identical(skel(to_mat(st[h, ])), skel(d)), TRUE))
  ndag <- sum(vapply(eq, function(h) acyclic(to_mat(st[h, ])), TRUE))
  fail <- NA
  if (length(eq) > 1) fail <- mean(replicate(N_GERSH, {                # A2 fallisce: un equivalente in B
    B2 <- draw_B(d); S2 <- lyapunov(B2, diag(p)); Si2 <- solve(S2); M2 <- fiber(S2)
    any(vapply(setdiff(eq, g), function(h) { Bp <- competitor(B2, M2, Si2, SUP[h, ]); !is.null(Bp) && dominant(Bp) }, TRUE)) }))
  cat_ <- if (length(eq) == 1) "identificabile" else if (setequal(eqlab, cl)) "spiegata" else
          if (all(cl %in% eqlab)) "altro (classe piu' grande)" else "chiusura piu' grande della classe"
  res[[r]] <- data.frame(p = p, grafo = lab(d), archi = NE[g], aciclico = acyclic(d), classe = length(eq),
                         membri_aciclici = ndag, chiusura_mosse = length(cl), stesso_scheletro = sameskel,
                         esito = cat_, chiusura_moduli = length(clm), moduli_spiegano = setequal(eqlab, clm),
                         frazione_A2_fallisce_in_B = fail,
                         membri = if (length(eq) > 1) paste(eqlab, collapse = " | ") else "")
}
res <- do.call(rbind, res)
dir.create("inst/simulations/identifiability/output", showWarnings = FALSE)
write.csv(res, sprintf("inst/simulations/identifiability/output/07_classes_p%d.csv", p), row.names = FALSE)
cat(sprintf("p = %d: %d grafi semplici, %d classi di isomorfismo, %.0f s\n", p, nrow(st), length(reps), proc.time()[3] - t0))
print(table(res$esito))
cat("classi non banali con scheletri diversi:", sum(res$classe > 1 & !res$stesso_scheletro), "\n")
cat(sprintf("cricche-modulo (c): classe = chiusura in %d classi su %d (non banali: %d su %d)\n",
            sum(res$moduli_spiegano), nrow(res), sum(res$moduli_spiegano & res$classe > 1), sum(res$classe > 1)))
nt <- res[res$classe > 1, ]
print(nt[order(nt$archi), c("grafo", "archi", "aciclico", "classe", "membri_aciclici", "chiusura_mosse", "esito", "frazione_A2_fallisce_in_B")], row.names = FALSE)
