# ============================================================
# Classi di equivalenza di modello e grafo di classe.
#
# Due grafi semplici hanno lo stesso modello (stesso insieme di Sigma) se si
# passa dall'uno all'altro riorientando cricche-modulo (Teorema "Clique-modules"
# del paper; caratterizzazione completa per p <= 5 e per grafi senza triangoli).
# Il grafo di classe tiene lo scheletro, orienta gli archi il cui verso e' lo
# stesso in tutta la classe e lascia non orientati gli altri.
#
# Convenzione: d[j, i] = 1 se i -> j. Nel grafo di classe un arco orientato
# i -> j ha cg[j, i] = 1, cg[i, j] = 0; un arco non orientato ha entrambe a 1.
# ============================================================

.graph_key <- function(d) paste0("g", paste(which(d == 1), collapse = ","))

.reach <- function(d) {                    # R[u, v]: u raggiunge v (u = v incluso)
  p <- nrow(d); R <- t(d) > 0 | diag(p) > 0
  for (k in seq_len(p)) R <- R | (R[, k] %o% R[k, ] > 0)
  R
}

.cliques <- function(s) {                  # tutte le cricche con almeno 2 nodi
  out <- list()
  ext <- function(Q, cand) for (v in cand) {
    Q2 <- c(Q, v); if (length(Q2) >= 2) out[[length(out) + 1]] <<- Q2
    ext(Q2, cand[cand > v & s[v, cand] == 1])
  }
  ext(integer(0), seq_len(nrow(s)))
  out
}

#' Cricche-modulo di un grafo semplice
#'
#' T e' una cricca-modulo se (M1) i suoi nodi sono tutti adiacenti, (M2) hanno
#' gli stessi genitori e gli stessi figli fuori da T, (M3) ogni genitore punta a
#' ogni figlio, (M4) nessun altro nodo ha un trek con T.
#' @param d matrice p x p 0/1, d[j, i] = 1 se i -> j
#' @return lista di vettori di nodi
clique_modules <- function(d) {
  p <- nrow(d); s <- (d | t(d)) * 1L; R <- .reach(d); out <- list()
  for (T_ in .cliques(s)) {
    pa <- lapply(T_, function(v) setdiff(which(d[v, ] == 1), T_))
    ch <- lapply(T_, function(v) setdiff(which(d[, v] == 1), T_))
    if (!all(vapply(pa, setequal, TRUE, pa[[1]])) || !all(vapply(ch, setequal, TRUE, ch[[1]]))) next
    P_ <- pa[[1]]; C_ <- ch[[1]]
    if (length(P_) && length(C_) && any(d[C_, P_] == 0)) next
    oth <- setdiff(seq_len(p), c(T_, P_, C_))
    trek <- function(k) any(vapply(T_, function(v) any(R[, k] & R[, v]), TRUE))
    if (any(vapply(oth, trek, TRUE))) next
    out[[length(out) + 1]] <- T_
  }
  out
}

.orient_all <- function(d, T_) {           # tutti gli orientamenti semplici dentro T
  pr <- which(upper.tri(diag(length(T_))), arr.ind = TRUE)
  lapply(0:(2^nrow(pr) - 1), function(o) {
    x <- d; x[T_, T_] <- 0L
    for (t in seq_len(nrow(pr))) { a <- T_[pr[t, 1]]; b <- T_[pr[t, 2]]
      if (bitwAnd(o, 2^(t - 1)) > 0) x[b, a] <- 1L else x[a, b] <- 1L }
    x
  })
}

#' Classe di equivalenza di modello di un grafo semplice
#'
#' Chiusura di d sotto il riorientamento delle cricche-modulo.
#' @param d matrice p x p 0/1 di un grafo semplice
#' @param max_size numero massimo di grafi (protezione per classi enormi)
#' @return lista di matrici
model_class <- function(d, max_size = 5000) {
  d <- (d != 0) * 1L; diag(d) <- 0L
  stopifnot(!any(d == 1 & t(d) == 1))
  seen <- new.env(hash = TRUE); assign(.graph_key(d), TRUE, seen)
  queue <- list(d); i <- 1
  while (i <= length(queue)) {
    y <- queue[[i]]; i <- i + 1
    for (T_ in clique_modules(y)) for (z in .orient_all(y, T_)) {
      k <- .graph_key(z)
      if (!exists(k, seen, inherits = FALSE)) {
        assign(k, TRUE, seen); queue[[length(queue) + 1]] <- z
        if (length(queue) > max_size) stop("classe con piu' di max_size grafi")
      }
    }
  }
  queue
}

#' Grafo di classe
#'
#' Per un grafo semplice: scheletro, archi orientati quando il verso e' lo
#' stesso in tutta la classe di modello, non orientati altrimenti. Per un grafo
#' con coppie reciproche (per esempio la stima di gclm): le coppie reciproche
#' diventano archi non orientati e la classe si calcola sul resto del grafo.
#' @param d matrice p x p 0/1, d[j, i] = 1 se i -> j
#' @return matrice p x p 0/1 del grafo di classe
class_graph <- function(d) {
  d <- (d != 0) * 1L; diag(d) <- 0L
  rec <- d == 1 & t(d) == 1
  ds <- d; ds[rec] <- 0L
  cl <- model_class(ds)
  cg <- ds
  for (x in cl) cg <- cg | x            # unione: un arco che cambia verso compare in entrambi i sensi
  cg <- cg * 1L; cg[rec] <- 1L
  cg
}

#' Metriche sul grafo di classe
#'
#' Ogni coppia ha uno stato: nessun arco, i -> j, j -> i, non orientato.
#' @param est grafo di classe stimato (da class_graph)
#' @param truth grafo di classe vero (da class_graph)
#' @return data.frame con tpr e fpr sullo scheletro; orient: frazione degli archi
#'   veri orientati nella classe che la stima orienta nello stesso verso;
#'   overorient: frazione degli archi veri non orientati nella classe che la stima
#'   orienta; shd: coppie con stato diverso; exact: shd == 0
class_metrics <- function(est, truth) {
  up <- upper.tri(est)
  state <- function(g) ifelse(!(g | t(g))[up], 0L, ifelse(g[up] & t(g)[up], 3L, ifelse(t(g)[up] == 1, 1L, 2L)))
  se <- state(est); st <- state(truth)
  data.frame(tpr_sk = mean(se[st > 0] > 0), fpr_sk = mean(se[st == 0] > 0),
             orient = if (any(st %in% 1:2)) mean(se[st %in% 1:2] == st[st %in% 1:2]) else NA_real_,
             overorient = if (any(st == 3)) mean(se[st == 3] %in% 1:2) else NA_real_,
             shd = sum(se != st), exact = all(se == st))
}
