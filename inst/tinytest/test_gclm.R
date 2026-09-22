## --- equazione di Lyapunov ------------------------------------------
set.seed(1)
p <- 6
B <- matrix(runif(p * p, -1, 1), p)
diag(B) <- 0
B <- build_B(B, rep(0.5, p))
S <- lyapunov(B, diag(p))

# risolve davvero B S + S B^T + C = 0
expect_true(max(abs(B %*% S + S %*% t(B) + diag(p))) < 1e-10)
expect_true(isSymmetric(S, tol = 1e-10))
expect_true(min(eigen(S, symmetric = TRUE, only.values = TRUE)$values) > 0)

# coincide con la soluzione via vettorizzazione di Kronecker
Sk <- matrix(solve(kronecker(diag(p), B) + kronecker(B, diag(p)), -as.vector(diag(p))), p, p)
expect_true(max(abs(S - Sk)) < 1e-8)

# tiene anche con B difettiva (blocco di Jordan), dove la via spettrale fallisce
J <- matrix(c(-1, 1, 0, 0, -1, 1, 0, 0, -1), 3, 3, byrow = TRUE)
SJ <- lyapunov(J, diag(3))
expect_true(max(abs(J %*% SJ + SJ %*% t(J) + diag(3))) < 1e-10)

## --- diagonale di Gershgorin ----------------------------------------
# B costruita cosi' e' sempre stabile, qualunque siano gli off-diagonali
set.seed(2)
for (k in 1:20) {
  M <- matrix(runif(25, -3, 3), 5)
  Bk <- build_B(M, runif(5, 0.1, 1))
  expect_true(max(Re(eigen(Bk)$values)) < 0)
}
expect_error(build_B(matrix(0, 3, 3), c(1, 1, -1)))   # a_i deve essere > 0

## --- il campionatore gira e restituisce oggetti coerenti -------------
tr <- make_true_model(p = 4, n_edge = 3, seed = 5)
expect_equal(dim(tr$Sigma), c(4L, 4L))
expect_true(all(diag(tr$B) < 0))

set.seed(3)
Y   <- rmvn0(200, tr$Sigma)
fit <- gclm_mcmc(Y, n_iter = 600, burnin = 300, thin = 3, verbose = FALSE)

expect_equal(dim(fit$B), c(4L, 4L, 100L))
expect_equal(dim(fit$PIP), c(4L, 4L))
expect_true(all(fit$PIP >= 0 & fit$PIP <= 1))
expect_true(all(diag(fit$PIP) == 0))             # nessun cappio
expect_true(all(is.finite(fit$loglik)))
expect_true(all(fit$tau > 0), all(fit$a > 0))

# la diagonale campionata rispetta il vincolo di Gershgorin
kk <- dim(fit$B)[3]
expect_true(all(sapply(seq_len(kk), function(s) {
  Bs <- fit$B[, , s]
  all(-diag(Bs) - (rowSums(abs(Bs)) - abs(diag(Bs))) > 0)
})))

# stesso seme, stessa catena
set.seed(99); f1 <- gclm_mcmc(Y, n_iter = 400, burnin = 200, thin = 2, verbose = FALSE)
set.seed(99); f2 <- gclm_mcmc(Y, n_iter = 400, burnin = 200, thin = 2, verbose = FALSE)
expect_identical(f1$PIP, f2$PIP)
expect_identical(f1$loglik, f2$loglik)

## --- input non validi ------------------------------------------------
expect_error(gclm_mcmc(Y, n_iter = 100, burnin = 200))          # burnin >= n_iter
expect_error(gclm_mcmc(Y, C = matrix(1, 4, 4)))                 # C non definita positiva
