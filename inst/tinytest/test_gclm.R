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

# ogni draw e' stabile (Gershgorin per costruzione)
expect_true(all(apply(fit$B, 3, function(Bs) max(Re(eigen(Bs, only.values = TRUE)$values)) < 0)))

## --- restrizione ai grafi senza 2-cicli -------------------------------
set.seed(4)
fl <- gclm_mcmc(Y, n_iter = 3000, burnin = 1000, thin = 2, verbose = FALSE)
fs <- restrict_simple(fl, min_draws = 1)
expect_true(fs$simple && fs$n_kept >= 1 && fs$frac_kept <= 1)
expect_true(all(apply(fs$delta, 3, function(d) !any(d == 1 & t(d) == 1))))
expect_equal(fs$PIP, apply(fs$delta, c(1, 2), mean))
expect_equal(dim(fs$B)[3], fs$n_kept)
expect_error(restrict_simple(fs))                                   # gia' ristretto

## --- campionatore per grafi semplici (prior sulle coppie) -------------
set.seed(6)
fsi <- gclm_mcmc(Y, n_iter = 3000, burnin = 1000, thin = 2, simple = TRUE, verbose = FALSE)
expect_true(fsi$simple && !is.null(fsi$rho))
expect_true(all(apply(fsi$delta, 3, function(d) !any(d == 1 & t(d) == 1))))   # nessun 2-ciclo
expect_true(all(apply(fsi$B, 3, function(Bs) max(Re(eigen(Bs, only.values = TRUE)$values)) < 0)))
expect_error(restrict_simple(fsi))
# prior sulle coppie con verosimiglianza spenta: E ~ beta-binomiale(P, a, b),
# stato della coppia: nessuno b/(a+b), ciascuna direzione a/(2(a+b))
q <- 4; Pq <- q * (q - 1) / 2
set.seed(7)
fsp <- gclm_mcmc(matrix(rnorm(40 * q), 40), n_iter = 60000, burnin = 5000, thin = 5,
                 a_rho = 2, b_rho = 3, simple = TRUE, prior_only = TRUE, verbose = FALSE)
Eq <- apply(fsp$delta, 3, sum)
expect_true(abs(mean(Eq) - Pq * 2 / 5) < 0.12)
expect_true(abs(mean(fsp$rho) - 2 / 5) < 0.02)
expect_true(abs(mean(fsp$a) - 1) < 0.05)
expect_true(all(apply(fsp$delta, 3, function(d) !any(d == 1 & t(d) == 1))))

## --- recupero della prior ---------------------------------------------
# Con verosimiglianza piatta la catena deve restituire la prior: e' il test dei
# rapporti di accettazione e dei jacobiani di tutte le mosse. Margini a_i ~
# Gamma(2, 2) indipendenti dagli archi; con pi ~ Beta(1, 1) il numero di archi
# e' uniforme su 0..K, quindi ha media K / 2; ristretta ai grafi senza 2-cicli,
# P(k) e' proporzionale a choose(m, k) 2^k B(k + 1, K - k + 1), m = K / 2 coppie.
q <- 4; K <- q * (q - 1); m <- K / 2
set.seed(5)
fp <- gclm_mcmc(matrix(rnorm(40 * q), 40), n_iter = 40000, burnin = 5000, thin = 5,
                prior_only = TRUE, verbose = FALSE)
a  <- as.vector(fp$a)
expect_true(abs(mean(a) - 1) < 0.05)
expect_true(abs(var(a) - 0.5) < 0.08)
expect_true(abs(mean(apply(fp$delta, 3, sum)) - K / 2) < 0.4)
rs <- t(apply(fp$B, 3, function(b) rowSums(abs(b)) - abs(diag(b))))
expect_true(max(abs(diag(cor(fp$a, rs)))) < 0.05)                  # a_i indipendente dagli archi
w  <- sapply(0:m, function(k) choose(m, k) * 2^k * beta(k + 1, K - k + 1))
ps <- restrict_simple(fp)
expect_true(abs(mean(apply(ps$delta, 3, sum)) - sum((0:m) * w) / sum(w)) < 0.15)
expect_true(abs(mean(ps$a) - 1) < 0.05)

## --- parallel tempering -----------------------------------------------
# con verosimiglianza spenta tutte le temperature hanno lo stesso bersaglio:
# gli scambi sono sempre accettati e la catena fredda restituisce la prior
bt <- c(1, 0.3, 0.1)
set.seed(8)
pa <- gclm_mcmc(matrix(rnorm(40 * q), 40), n_iter = 30000, burnin = 3000, thin = 3,
                a_rho = 2, b_rho = 3, simple = TRUE, eps_bridge = 0.1, pt_betas = bt,
                prior_only = TRUE, verbose = FALSE)
expect_true(all(pa$swap == 1))
expect_true(abs(mean(apply(pa$delta, 3, sum)) - Pq * 2 / 5) < 0.15)
expect_true(all(apply(pa$delta, 3, function(d) !any(d == 1 & t(d) == 1))))
set.seed(9)
pf <- gclm_mcmc(matrix(rnorm(40 * q), 40), n_iter = 30000, burnin = 3000, thin = 3,
                pt_betas = bt, prior_only = TRUE, verbose = FALSE)
expect_true(abs(mean(apply(pf$delta, 3, sum)) - q * (q - 1) / 2) < 0.3)   # pi ~ Beta(1, 1)
expect_error(gclm_mcmc(Y, pt_betas = c(0.5, 0.2)))                          # deve iniziare da 1
# con i dati: stabilita' e nessun 2-ciclo nella catena fredda di A1
set.seed(10)
pd <- gclm_mcmc(Y, n_iter = 1500, burnin = 500, thin = 2, simple = TRUE, pt_betas = bt, verbose = FALSE)
expect_true(all(apply(pd$B, 3, function(Bs) max(Re(eigen(Bs, only.values = TRUE)$values)) < 0)))
expect_true(all(apply(pd$delta, 3, function(d) !any(d == 1 & t(d) == 1))))

## --- input non validi ------------------------------------------------
expect_error(gclm_mcmc(Y, n_iter = 100, burnin = 200))          # burnin >= n_iter
expect_error(gclm_mcmc(Y, C = matrix(1, 4, 4)))                 # C non definita positiva
expect_error(gclm_mcmc(Y, simple = "si"))

# ---- regola di selezione per coppia ------------------------------------
fk <- list(PIP = matrix(0, 3, 3), simple = TRUE)
fk$PIP[2, 1] <- 0.4; fk$PIP[1, 2] <- 0.3    # coppia {1,2}: arco con prob. 0.7, verso 1 -> 2
fk$PIP[3, 1] <- 0.2; fk$PIP[1, 3] <- 0.2    # coppia {1,3}: prob. 0.4, non selezionata
sg <- select_graph(fk)
expect_equal(sum(sg), 1L)
expect_equal(sg[2, 1], 1L)
expect_equal(sum(select_graph(fk, rule = "arrow")), 0L)
fk$simple <- FALSE
expect_equal(select_graph(fk), select_graph(fk, rule = "arrow"))
