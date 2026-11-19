# ============================================================
# Controlli numerici delle formule usate nelle dimostrazioni
# (Ricerca/BGCLM/paper/dimostrazioni.tex).
# 1. Termine di secondo ordine di Sigma sulle coppie non adiacenti (Lemma 2):
#    l'errore deve scalare come w.
# 2. Costruzione esplicita del Teorema 1 (cricca-modulo): B' ha supporto in G',
#    la stessa Sigma ed e' stabile.
# 3. Corollario 2: KL* / w^4 per la catena 1 -> 2 -> 3 invertita, contro la formula.
# Uso:  Rscript inst/simulations/identifiability/10_check_proofs.R
# ============================================================
suppressMessages(library(BGCLM)); set.seed(11)
cat("== 1. Termine di secondo ordine sulle coppie non adiacenti ==\n")
p <- 6; err <- c()
for (r in 1:20) {
  d <- matrix(0L, p, p); pr <- which(upper.tri(diag(p)), arr.ind = TRUE); s <- sample(nrow(pr), 7)
  for (q in s) if (runif(1) < .5) d[pr[q, 2], pr[q, 1]] <- 1L else d[pr[q, 1], pr[q, 2]] <- 1L
  U <- d * rnorm(p^2); D <- runif(p, 0.3, 1.5); Dl <- rnorm(p)            # anche una parte diagonale O(w)
  S0 <- diag(1 / (2 * D)); S1 <- (U %*% S0 + S0 %*% t(U)) / outer(D, D, "+")
  S2f <- matrix(0, p, p)                                                   # formula: tipi della tripla
  sk <- (d | t(d))
  for (i in 1:p) for (j in 1:p) if (i != j && !sk[i, j]) for (k in which(sk[i, ] & sk[j, ])) {
    m_ik <- U[i, k] / D[k] + U[k, i] / D[i]; m_kj <- U[k, j] / D[j] + U[j, k] / D[k]
    c_ <- (if (U[i, k] != 0 && U[k, j] != 0) D[k] / (2 * (D[k] + D[j])) else 0) +      # k -> i, j -> k
          (if (U[k, i] != 0 && U[j, k] != 0) D[k] / (2 * (D[k] + D[i])) else 0) +      # i -> k, k -> j
          (if (U[i, k] != 0 && U[j, k] != 0) D[k] / (2 * (D[k] + D[j])) + D[k] / (2 * (D[k] + D[i])) else 0)  # fork
    S2f[i, j] <- S2f[i, j] + m_ik * m_kj * c_ / (D[i] + D[j]) }
  for (w in c(1e-3, 1e-4)) { B <- -diag(D) + w * (U + diag(Dl)); S <- lyapunov(B, diag(p))
    nad <- !sk & diag(p) == 0; err <- rbind(err, c(w, max(abs((S[nad] / w^2) - S2f[nad])))) }
}
for (w in c(1e-3, 1e-4)) cat(sprintf("w = %.0e: errore massimo |Sigma_ij / w^2 - formula| sulle coppie non adiacenti, 20 grafi: %.1e\n", w, max(err[err[, 1] == w, 2])))

cat("\n== 2. Costruzione esplicita per una cricca-modulo ==\n")
build <- function(B, T_, P_, Ch, newT) {          # newT: matrice 0/1 dell'orientamento nuovo dentro T
  p <- nrow(B); S <- lyapunov(B, diag(p)); Sset <- c(P_, T_); Sss <- solve(S[Sset, Sset])
  t <- length(T_); pr <- which(upper.tri(diag(t)), arr.ind = TRUE)
  # incognite K_T (t(t-1)/2); vincolo: zero nella posizione non-arco del nuovo orientamento
  f <- function(k) { K <- matrix(0, t, t); K[pr] <- k; K <- K - t(K)
    Bn <- B; Bn[T_, Sset] <- B[T_, Sset] + cbind(matrix(0, t, length(P_)), K) %*% Sss; Bn }
  A <- sapply(seq_len(nrow(pr)), function(q) { e <- numeric(nrow(pr)); e[q] <- 1
    Bn <- f(e) - f(0 * e); sapply(seq_len(nrow(pr)), function(r) { a <- pr[r, 1]; b <- pr[r, 2]
      if (newT[b, a] == 1) Bn[T_[a], T_[b]] else Bn[T_[b], T_[a]] }) })      # posizione da annullare
  B0 <- f(numeric(nrow(pr))); b0 <- sapply(seq_len(nrow(pr)), function(r) { a <- pr[r, 1]; b <- pr[r, 2]
      if (newT[b, a] == 1) B0[T_[a], T_[b]] else B0[T_[b], T_[a]] })
  k <- solve(matrix(A, nrow(pr)), -b0); K <- matrix(0, t, t); K[pr] <- k; K <- K - t(K)
  Bn <- f(k); beta <- Sss %*% S[Sset, Ch, drop = FALSE]; bT <- beta[(length(P_) + 1):length(Sset), , drop = FALSE]
  for (ci in seq_along(Ch))                                                # righe dei figli comuni
    Bn[Ch[ci], Sset] <- B[Ch[ci], Sset] - c(rep(0, length(P_)), K %*% bT[, ci]) %*% Sss
  list(B = Bn, S = S)
}
mk <- function(p, e) { d <- matrix(0L, p, p); d[e[, 2:1, drop = FALSE]] <- 1L; d }
E <- function(...) matrix(c(...), ncol = 2, byrow = TRUE)
casi <- list(
  list(nome = "triangolo {2,3,4} -> figlio 1, 5 -> 1 (5 senza trek)", p = 5, G = E(2,3, 3,4, 2,4, 2,1, 3,1, 4,1, 5,1), T = c(2,3,4), P = integer(0), Ch = 1, cic = TRUE),
  list(nome = "genitori 4,5 -> triangolo {1,2,3} (in 3-ciclo)", p = 5, G = E(1,2, 2,3, 1,3, 4,1,4,2,4,3, 5,1,5,2,5,3), T = 1:3, P = 4:5, Ch = integer(0), cic = TRUE),
  list(nome = "genitore 6, figli 4,5, T = {1,2,3}, 7 -> 4 e 6 -> 4,5", p = 7, G = E(1,2, 2,3, 1,3, 6,1,6,2,6,3, 1,4,2,4,3,4, 1,5,2,5,3,5, 6,4,6,5, 4,5, 7,4), T = 1:3, P = 6, Ch = 4:5, cic = TRUE),
  list(nome = "super-covered 1 -> 4 (random rep 6)", p = 5, G = E(1,4, 3,1, 3,4, 5,1, 5,4), T = c(1,4), P = c(3,5), Ch = integer(0), cic = FALSE))
for (cs in casi) { res <- c()
  for (r in 1:10) {
    d <- mk(cs$p, cs$G); B <- build_B(d * sample(c(-1, 1), cs$p^2, TRUE) * runif(cs$p^2, .3, 1.5), runif(cs$p, .2, 1))
    t <- length(cs$T)
    newT <- if (cs$cic && t == 3) { m <- matrix(0, 3, 3); m[2, 1] <- m[3, 2] <- m[1, 3] <- 1; m } else { m <- 1 - (d[cs$T, cs$T]); diag(m) <- 0; m }  # 3-ciclo, o tutto invertito
    out <- build(B, cs$T, cs$P, cs$Ch, newT)
    Gp <- d; Gp[cs$T, cs$T] <- newT; Bn <- out$B
    fuori <- max(abs(Bn[Gp == 0 & diag(cs$p) == 0]))
    res <- rbind(res, c(fuori = fuori, errore_Sigma = max(abs(lyapunov(Bn, diag(cs$p)) - out$S)),
                        stabile = max(Re(eigen(Bn, only.values = TRUE)$values)) < 0)) }
  cat(sprintf("%s: max |B'| fuori da G' %.1e, max errore su Sigma %.1e, stabile in %d/10\n", cs$nome, max(res[, 1]), max(res[, 2]), sum(res[, 3])))
}

cat("\n== 3. Catena invertita con diagonale diversa: KL / w^4 contro (d1 - d3)^2 ==\n")
kl <- function(S0, S) 0.5 * (sum(diag(solve(S, S0))) - nrow(S) + determinant(S)$modulus - determinant(S0)$modulus)
min_kl <- function(S0, sup, warm) { idx <- which(sup == 1)
  f <- function(th) { B <- matrix(0, 3, 3); B[idx] <- th[seq_along(idx)]; diag(B) <- -exp(th[-seq_along(idx)])
    if (max(Re(eigen(B, only.values = TRUE)$values)) >= 0) return(1e6); as.numeric(kl(S0, lyapunov(B, diag(3)))) }
  optim(warm, f, method = "BFGS", control = list(maxit = 20000, reltol = 1e-16, ndeps = rep(1e-7, length(warm))))$value }
w <- 0.02; tab <- c()
for (a1 in c(0.5, 0.55, 0.6, 0.7, 0.9)) {
  a <- c(a1, 0.5, 0.5); Bo <- matrix(0, 3, 3); Bo[2, 1] <- w; Bo[3, 2] <- w; B <- build_B(Bo, a); S0 <- lyapunov(B, diag(3))
  sup <- matrix(0, 3, 3); sup[1, 2] <- 1; sup[2, 3] <- 1; dd <- -diag(B)
  warm <- c(w * dd[1] / dd[2], w * dd[2] / dd[3], log(dd))
  tab <- rbind(tab, c(d1_meno_d3 = dd[1] - dd[3], kl_su_w4 = min_kl(S0, sup, warm) / w^4, previsto = dd[3] * (dd[1] - dd[3])^2 / (2 * dd[1] * (dd[1] + dd[2])^2 * (dd[2] + dd[3])^2 * (dd[1] + dd[3])^2))) }
print(signif(tab, 3))
