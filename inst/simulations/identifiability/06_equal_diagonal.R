# Diagonale uguale (B_ii = -d per ogni i, d ignoto) al posto di A2.
# Ogni grafo semplice e' contenuto in un torneo T (un verso per ogni coppia),
# quindi basta controllare i 2^(p(p-1)/2) tornei: per ciascuno si cerca K
# antisimmetrica con B + K Sigma^-1 nullo sulle entrate escluse da T e
# diag(K Sigma^-1) costante (= 0, perche' tr(K Sigma^-1) = 0).
#  - T contiene il grafo vero: sistema omogeneo, serve rango pieno (solo K = 0);
#  - T non lo contiene: il sistema deve essere incompatibile.
# Lo stesso conto senza il vincolo sulla diagonale mostra che ogni torneo ammette
# una deriva equivalente. Tutti i grafi semplici su p = 5 nodi, a meno di
# isomorfismo, NDRAW estrazioni dei valori per grafo.
# Uso: Rscript inst/simulations/identifiability/06_equal_diagonal.R
source("inst/simulations/identifiability/common_ident.R")
NDRAW <- 3
m <- nrow(PAIRS)                                         # 10 coppie (r < c)

# ---- classi di isomorfismo dei grafi semplici ----------------------------
# stato della coppia q = (r, c): 0 nessun arco, 1 arco r -> c (B[c, r]), 2 arco c -> r
states <- as.matrix(expand.grid(rep(list(0:2), m)))
perms <- as.matrix(expand.grid(rep(list(1:p), p))); perms <- perms[apply(perms, 1, function(x) length(unique(x)) == p), ]
code <- function(S) as.vector(S %*% 3^(0:(m - 1)))
canon <- code(states)
for (k in seq_len(nrow(perms))) {
  s <- perms[k, ]; S2 <- matrix(0L, nrow(states), m)
  for (q in seq_len(m)) {
    a <- s[PAIRS[q, 1]]; b <- s[PAIRS[q, 2]]
    q2 <- which(PAIRS[, 1] == min(a, b) & PAIRS[, 2] == max(a, b))
    st <- states[, q]; if (a > b) st <- ifelse(st == 1, 2L, ifelse(st == 2, 1L, 0L))
    S2[, q2] <- st
  }
  canon <- pmin(canon, code(S2))
}
reps <- which(!duplicated(canon))
cat("grafi semplici etichettati:", nrow(states), "| classi di isomorfismo:", length(reps), "\n")

to_delta <- function(st) { d <- matrix(0L, p, p)
  for (q in seq_len(m)) { r <- PAIRS[q, 1]; c <- PAIRS[q, 2]
    if (st[q] == 1) d[c, r] <- 1L else if (st[q] == 2) d[r, c] <- 1L }
  d }

# ---- tornei: entrate (posizioni in OFF) che T obbliga a zero --------------
tours <- as.matrix(expand.grid(rep(list(1:2), m)))
zero_idx <- lapply(seq_len(nrow(tours)), function(t) {
  dz <- to_delta(ifelse(tours[t, ] == 1, 2L, 1L))       # il verso escluso
  which(dz[OFF] == 1) })
tour_delta <- lapply(seq_len(nrow(tours)), function(t) to_delta(tours[t, ]))

diag_rows <- function(Sigma) { Si <- solve(Sigma)
  sapply(seq_len(m), function(q) { K <- matrix(0, p, p); K[PAIRS[q, 1], PAIRS[q, 2]] <- 1
    K[PAIRS[q, 2], PAIRS[q, 1]] <- -1; diag(K %*% Si) }) }

res <- list()
for (g in reps) for (dd in seq_len(NDRAW)) {
  delta <- to_delta(states[g, ]); set.seed(1e5 + 7 * g + dd)
  Bo <- matrix(0, p, p); e <- which(delta == 1)
  Bo[e] <- sample(c(-1, 1), length(e), TRUE) * runif(length(e), 0.4, 0.8)
  d <- max(rowSums(abs(Bo))) + runif(1, 0.2, 1)           # stabile per Gershgorin
  B <- Bo - d * diag(p); Sigma <- lyapunov(B, diag(p))
  M <- fiber_matrix(Sigma); Dm <- diag_rows(Sigma); Dd <- Dm[-1, ] - matrix(Dm[1, ], p - 1, m, byrow = TRUE)
  b <- B[OFF]
  min_res <- Inf; min_sv_cont <- Inf; n_free_ok <- 0
  for (t in seq_along(zero_idx)) {
    Z <- zero_idx[[t]]; contains <- all(tour_delta[[t]][delta == 1] == 1)
    A <- rbind(M[Z, ], Dd); y <- c(-b[Z], rep(0, p - 1))
    if (contains) {                                        # omogeneo: serve rango pieno
      min_sv_cont <- min(min_sv_cont, min(svd(A)$d))
    } else {                                               # deve essere incompatibile
      qa <- qr(A); r <- sqrt(sum(qr.resid(qa, y)^2)) / sqrt(sum(y^2))
      min_res <- min(min_res, r)
    }
    qf <- qr(M[Z, ]); rf <- sqrt(sum(qr.resid(qf, -b[Z])^2))   # senza vincolo sulla diagonale
    n_free_ok <- n_free_ok + (rf < 1e-10)
  }
  res[[length(res) + 1]] <- data.frame(classe = match(g, reps), estrazione = dd, archi = sum(delta),
    d = d, d_da_Sigma = sum(diag(solve(Sigma))) / (2 * p),
    min_sv_tornei_contenenti = min_sv_cont, min_residuo_rel_altri = min_res,
    tornei_compatibili_diag_libera = n_free_ok)
}
res <- do.call(rbind, res)
write.csv(res, "inst/simulations/identifiability/output/06_equal_diagonal.csv", row.names = FALSE)
cat("casi (classe x estrazione):", nrow(res), "\n")
cat("max |d - tr(Sigma^-1)/(2p)|:", signif(max(abs(res$d - res$d_da_Sigma)), 3), "\n")
cat("diagonale uguale - min valore singolare (tornei che contengono il vero):", signif(min(res$min_sv_tornei_contenenti), 3), "\n")
cat("diagonale uguale - min residuo relativo (tornei che non lo contengono):", signif(min(res$min_residuo_rel_altri), 3), "\n")
cat("diagonale libera - tornei con una deriva equivalente: min", min(res$tornei_compatibili_diag_libera),
    "max", max(res$tornei_compatibili_diag_libera), "su", nrow(tours), "\n")
print(aggregate(cbind(min_residuo_rel_altri, min_sv_tornei_contenenti) ~ archi, res, min), digits = 3)
