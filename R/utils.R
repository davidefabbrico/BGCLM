# ============================================================
# Utilita' per il modello GCLM: costruzione di B, equazione di Lyapunov,
# generazione dei dati, diagnostica della catena e valutazione del fit.
#
# Convenzione sugli indici: B[j, i] != 0  <=>  arco i -> j
# ============================================================

#' Risolve l'equazione di Lyapunov B Sigma + Sigma B^T + C = 0
#'
#' @param B matrice p x p, deve essere stabile (autovalori con parte reale < 0)
#' @param C matrice p x p simmetrica definita positiva (default: identita')
#' @return la matrice Sigma p x p
lyapunov <- function(B, C = diag(nrow(B))) {
  B <- as.matrix(B); C <- as.matrix(C)
  stopifnot(nrow(B) == ncol(B), dim(B) == dim(C))
  lyapunov_cpp(B, C)
}

#' Costruisce B dagli off-diagonali imponendo la diagonale di Gershgorin
#'
#' B_ii = -( sum_{j != i} |B_ij| + a_i ), che garantisce autovalori con parte
#' reale <= -a_i < 0, quindi stabilita' per qualunque valore degli off-diagonali.
#'
#' @param Boff matrice p x p degli off-diagonali (la diagonale viene ignorata)
#' @param a vettore di p margini di stabilita', tutti > 0
build_B <- function(Boff, a) {
  B <- as.matrix(Boff)
  diag(B) <- 0
  stopifnot(length(a) == nrow(B), all(a > 0))
  diag(B) <- -(rowSums(abs(B)) + a)
  B
}

#' Genera un modello GCLM vero, sparso e con un ciclo
#'
#' Il grafo contiene sempre il ciclo 1 -> 2 -> 3 -> 1 (il modello ammette cicli,
#' non e' un DAG). Gli archi extra escludono le coppie reciproche i <-> j, che
#' sono poco identificabili: scambiare le due direzioni cambia poco Sigma.
#'
#' @param p numero di nodi (>= 3)
#' @param n_edge numero di archi veri
#' @param a_true margine di stabilita' comune a tutti i nodi
#' @param seed seme per la riproducibilita'
make_true_model <- function(p = 5, n_edge = 5, a_true = 0.5, seed = 1) {
  stopifnot(p >= 3, n_edge >= 3)
  set.seed(seed)
  Boff <- matrix(0, p, p)
  Boff[rbind(c(2, 1), c(3, 2), c(1, 3))] <- c(0.6, -0.5, 0.7)
  for (k in seq_len(n_edge - 3)) {
    free <- which(diag(p) == 0 & Boff == 0 & t(Boff) == 0, arr.ind = TRUE)
    if (!nrow(free)) break
    Boff[free[sample(nrow(free), 1), , drop = FALSE]] <- sample(c(-1, 1), 1) * runif(1, 0.4, 0.8)
  }
  a <- rep(a_true, p)
  B <- build_B(Boff, a)
  list(Boff = Boff, a = a, B = B, delta = (Boff != 0) * 1L,
       Sigma = lyapunov(B, diag(p)))
}

#' Estrae n osservazioni iid dall'equilibrio N_p(0, Sigma)
rmvn0 <- function(n, Sigma) matrix(rnorm(n * ncol(Sigma)), n) %*% chol(Sigma)

#' Dimensione campionaria efficace (somma delle autocorrelazioni fino al primo lag negativo)
ess <- function(x) {
  x <- as.numeric(x); m <- length(x)
  if (stats::sd(x) < 1e-12) return(0)
  a <- stats::acf(x, lag.max = min(m - 1, 500), plot = FALSE)$acf[-1]
  k <- which(a < 0)[1]
  k <- if (is.na(k)) length(a) else k - 1
  if (k < 1) return(m)
  m / (1 + 2 * sum(a[1:k]))
}

#' Area sotto la curva ROC (statistica di Mann-Whitney)
auc <- function(score, label) {
  n1 <- sum(label == 1); n0 <- sum(label == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  (sum(rank(score)[label == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

#' Confronta un fit con il modello vero sugli archi (solo off-diagonali)
#'
#' @param fit output di gclm_mcmc
#' @param truth output di make_true_model
#' @param thr soglia della regola di selezione, vedi [select_graph()]
eval_fit <- function(fit, truth, thr = 0.5) {
  p   <- ncol(truth$B)
  off <- which(diag(p) == 0)
  lab <- truth$delta[off]
  pip <- fit$PIP[off]
  sel <- select_graph(fit, thr = thr)[off]
  qs  <- apply(fit$B, c(1, 2), stats::quantile, probs = c(0.025, 0.975))
  data.frame(
    auc     = auc(pip, lab),
    tpr     = if (any(lab == 1)) mean(sel[lab == 1] == 1) else NA_real_,
    fpr     = if (any(lab == 0)) mean(sel[lab == 0] == 1) else NA_real_,
    rmse_B  = sqrt(mean((fit$B_mean[off] - truth$B[off])^2)),
    err_Sig = norm(fit$Sigma_mean - truth$Sigma, "F") / norm(truth$Sigma, "F"),
    cov95   = mean(truth$B[off] >= qs[1, , ][off] & truth$B[off] <= qs[2, , ][off]),
    n_edge  = mean(apply(fit$delta, 3, sum))
  )
}

#' Tabella degli archi ordinata per probabilita' a posteriori di inclusione
#'
#' post_m e' la media marginale (schiacciata verso 0 dal peso dello spike),
#' post_in la media condizionata all'inclusione: e' quest'ultima a stimare
#' l'ampiezza dell'arco.
edge_table <- function(fit, truth = NULL) {
  p   <- nrow(fit$PIP)
  off <- which(diag(p) == 0, arr.ind = TRUE)
  qs  <- apply(fit$B, c(1, 2), stats::quantile, probs = c(0.025, 0.975))
  out <- data.frame(
    arco    = sprintf("%d -> %d", off[, 2], off[, 1]),
    PIP     = round(fit$PIP[off], 3),
    post_m  = round(fit$B_mean[off], 3),
    post_in = round(fit$B_incl[off], 3),
    lo      = round(qs[1, , ][off], 3),
    hi      = round(qs[2, , ][off], 3)
  )
  if (!is.null(truth)) out <- cbind(out[1], vero = round(truth$B[off], 3), out[-1])
  out[order(-out$PIP), ]
}
