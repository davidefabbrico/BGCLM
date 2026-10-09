# ============================================================
# Funzioni comuni per la verifica numerica delle assunzioni A1 e A2.
#
# Convenzione: l'arco i -> j e' B_ji, cioe' B[j, i].
# Ogni deriva con la stessa Sigma (C = I) e' B' = B + K Sigma^-1, K antisimmetrica
# (Prop. "Fibre"): imporre zeri fuori diagonale e' un sistema lineare nei
# p(p-1)/2 parametri di K.
# ============================================================
suppressMessages(library(BGCLM))
src <- readLines("inst/simulations/compare_gclm.R")
eval(parse(text = paste(src[1:(grep("^GEN <- list", src) - 1)], collapse = "\n")))

p <- P
OFF <- which(diag(p) == 0)
PAIRS <- which(upper.tri(diag(p)), arr.ind = TRUE)          # coppie non ordinate (r < c)
SEEDS <- c(gershgorin = 1000, weak = 3000)                   # seme = SEEDS[scen] + 10 * rep

true_model <- function(scen, rep) {
  set.seed(SEEDS[[scen]] + 10 * rep)
  get(paste0("gen_", scen))(p, N_EDGE)
}

edges_lab <- function(d) {
  w <- which(d == 1, arr.ind = TRUE); w <- w[order(w[, 2], w[, 1]), , drop = FALSE]
  paste(sprintf("%d->%d", w[, 2], w[, 1]), collapse = " ")
}

# colonna q = parte fuori diagonale di K_q Sigma^-1, con K_q = E_rc - E_cr
fiber_matrix <- function(Sigma) {
  Si <- solve(Sigma)
  sapply(seq_len(nrow(PAIRS)), function(q) {
    K <- matrix(0, p, p); K[PAIRS[q, 1], PAIRS[q, 2]] <- 1; K[PAIRS[q, 2], PAIRS[q, 1]] <- -1
    as.vector(K %*% Si)[OFF] })
}
fiber_full <- function(Sigma, k) {        # matrice K Sigma^-1 completa per un vettore k
  K <- matrix(0, p, p)
  for (q in seq_along(k)) { K[PAIRS[q, 1], PAIRS[q, 2]] <- k[q]; K[PAIRS[q, 2], PAIRS[q, 1]] <- -k[q] }
  K %*% solve(Sigma)
}

is_dominant <- function(B) all(-diag(B) - (rowSums(abs(B)) - abs(diag(B))) > 1e-10)
is_stable   <- function(B) max(Re(eigen(B, only.values = TRUE)$values)) < -1e-10

# tutti i supporti SEMPLICI con al piu' kmax archi: per ogni insieme di coppie,
# ogni orientamento. Restituisce una lista di vettori logici sulle posizioni OFF.
simple_supports <- function(kmax) {
  out <- list()
  for (k in 1:kmax) {
    cmb <- combn(nrow(PAIRS), k)
    for (s in seq_len(ncol(cmb))) for (o in 0:(2^k - 1)) {
      d <- matrix(FALSE, p, p)
      for (t in seq_len(k)) {
        r <- PAIRS[cmb[t, s], 1]; c <- PAIRS[cmb[t, s], 2]
        if (bitwAnd(o, 2^(t - 1)) > 0) d[c, r] <- TRUE else d[r, c] <- TRUE
      }
      out[[length(out) + 1]] <- d[OFF]
    }
  }
  out
}

# censimento: supporti semplici per cui esiste B' = B + K Sigma^-1 con supporto
# esattamente S'. Per ognuno: stabile? nella classe Gershgorin (a_i > 0)?
census_simple <- function(B, Sigma, delta, supports) {
  M <- fiber_matrix(Sigma); b <- B[OFF]; out <- list()
  for (s in supports) {
    Z <- which(!s)
    fit <- lm.fit(M[Z, , drop = FALSE], -b[Z])
    if (sum(fit$residuals^2) > 1e-20) next
    k <- fit$coefficients; k[is.na(k)] <- 0
    Bp <- B + fiber_full(Sigma, k)
    d <- (abs(Bp) > 1e-9) * 1L; diag(d) <- 0L
    if (!identical(as.vector(d[OFF] == 1), s)) next          # supporto esatto
    out[[length(out) + 1]] <- data.frame(
      grafo = edges_lab(d), archi = sum(d), vero = all(d == delta),
      rango_Z = qr(M[Z, , drop = FALSE])$rank,
      stabile = is_stable(Bp), gershgorin = is_dominant(Bp),
      residuo_Sigma = if (is_stable(Bp)) max(abs(lyapunov(Bp, diag(p)) - Sigma)) else NA)
  }
  if (length(out)) do.call(rbind, out) else NULL
}
