# ============================================================
# Quanta informazione c'e' sull'orientamento.
# (1) Catena vera 1 -> 2 -> 3 con archi di intensita' w. Per la catena invertita,
#     il fork e il collider: KL* = min KL(N(0, Sigma) || N(0, Sigma(B'))) sui drift
#     stabili con quel supporto (diagonale libera, C = I). Pendenza locale di
#     log KL* in log w e n per 10 nat (= 10 / KL*).
# (2) Arco isolato 1 -> 2 (classe di equivalenza esatta): MPPI del verso vero al
#     crescere di n, quando il drift invertito e' in B e quando non lo e'.
# Uso:  Rscript inst/simulations/identifiability/08_orientation_information.R
# Scrive inst/simulations/identifiability/output/08_*.csv
# ============================================================
suppressMessages(library(BGCLM))
kl <- function(S0, S) 0.5 * (sum(diag(solve(S, S0))) - nrow(S) + determinant(S)$modulus - determinant(S0)$modulus)
min_kl <- function(S0, sup, starts = 15) {
  p <- nrow(S0); idx <- which(sup == 1)
  f <- function(th) { B <- matrix(0, p, p); B[idx] <- th[seq_along(idx)]; diag(B) <- -exp(th[-seq_along(idx)])
    if (max(Re(eigen(B, only.values = TRUE)$values)) >= -1e-10) return(1e6)
    v <- tryCatch(as.numeric(kl(S0, lyapunov(B, diag(p)))), error = function(e) 1e6); if (is.finite(v)) v else 1e6 }
  min(sapply(seq_len(starts), function(s) { set.seed(s)
    optim(c(rnorm(length(idx), 0, 0.7), log(runif(p, 0.3, 2))), f, method = "BFGS",
          control = list(maxit = 5000, reltol = 1e-15))$value }))
}

# ---------------- (1) ----------------
# archi come righe (da, a); al secondo ordine conta, per ogni coppia con un vicino
# comune k, il tipo di k: collider 0, catena 1 (in entrambi i versi), fork 2
E <- function(...) matrix(c(...), ncol = 2, byrow = TRUE)
cases1 <- list(
  "catena 1>2>3" = list(true = E(1, 2, 2, 3), alt = list(catena_invertita = E(3, 2, 2, 1), fork = E(2, 1, 2, 3), collider = E(1, 2, 3, 2))),
  "cammino 1>2>3>4" = list(true = E(1, 2, 2, 3, 3, 4), alt = list(cammino_invertito = E(4, 3, 3, 2, 2, 1), inverti_solo_3_4 = E(1, 2, 2, 3, 4, 3))),
  "hub 1>2, 2>3, 2>4" = list(true = E(1, 2, 2, 3, 2, 4), alt = list(inverti_1_2 = E(2, 1, 2, 3, 2, 4), tutti_invertiti = E(2, 1, 3, 2, 4, 2))))
W <- c(0.1, 0.2, 0.4, 0.8, 1.6); r1 <- list()
for (cn in names(cases1)) for (w in W) {
  tr <- cases1[[cn]]$true; p <- max(tr); Bo <- matrix(0, p, p); Bo[tr[, 2:1, drop = FALSE]] <- w
  S0 <- lyapunov(build_B(Bo, rep(0.5, p)), diag(p))
  for (a in names(cases1[[cn]]$alt)) { sup <- matrix(0, p, p); sup[cases1[[cn]]$alt[[a]][, 2:1, drop = FALSE]] <- 1
    r1[[length(r1) + 1]] <- data.frame(vero = cn, w = w, alternativa = a, corr_12 = cov2cor(S0)[1, 2], kl = max(min_kl(S0, sup), 0)) }
}
r1 <- do.call(rbind, r1); r1$n_per_10_nat <- signif(10 / r1$kl, 3)
r1$pendenza <- ave(log(r1$kl), r1$vero, r1$alternativa, FUN = function(z) c(NA, diff(z) / diff(log(W))))
write.csv(r1, "inst/simulations/identifiability/output/08_kl_orientation.csv", row.names = FALSE)
cat("(1) KL* delle alternative sullo stesso scheletro, archi di intensita' w, margini 0.5\n")
print(transform(r1, kl = signif(kl, 3), corr_12 = round(corr_12, 2), pendenza = round(pendenza, 1)), row.names = FALSE)

# ---------------- (2) ----------------
cases <- list(rovescio_in_B = list(b = 0.6, a = c(1, 1)), rovescio_fuori_da_B = list(b = 2, a = c(0.3, 0.5)))
dominant <- function(B) all(-diag(B) - (rowSums(abs(B)) - abs(diag(B))) > 1e-10)
r2 <- list()
for (cs in names(cases)) {
  Bo <- matrix(0, 2, 2); Bo[2, 1] <- cases[[cs]]$b; B <- build_B(Bo, cases[[cs]]$a); S0 <- lyapunov(B, diag(2))
  Si <- solve(S0); k <- B[2, 1] / Si[1, 1]                            # B' = B + K Sigma^-1 con B'[2,1] = 0
  Bp <- B + matrix(c(0, -k, k, 0), 2) %*% Si
  cat(sprintf("\n(2) %s: B'[1,2] = %.2f, diagonale di B' = %.2f %.2f, B' in B: %s, residuo su Sigma %.1e\n", cs, Bp[1, 2],
              Bp[1, 1], Bp[2, 2], dominant(Bp), max(abs(lyapunov(Bp, diag(2)) - S0))))
  for (n in c(200, 2000, 20000)) for (s in 1:2) {
    set.seed(100 * s + n %% 97); Y <- matrix(rnorm(n * 2), n) %*% chol(S0)
    f <- gclm_mcmc(Y, n_iter = 6000, burnin = 2000, thin = 2, simple = TRUE, eps_bridge = 0.1,
                   pt_betas = exp(seq(0, log(0.01), length.out = 10)), verbose = FALSE)
    r2[[length(r2) + 1]] <- data.frame(caso = cs, n = n, seme = s, P_coppia = f$PIP[2, 1] + f$PIP[1, 2],
                                       P_verso_vero = f$PIP[2, 1] / (f$PIP[2, 1] + f$PIP[1, 2]))
  }
}
r2 <- do.call(rbind, r2)
write.csv(r2, "inst/simulations/identifiability/output/08_isolated_edge_mcmc.csv", row.names = FALSE)
print(transform(r2, P_coppia = round(P_coppia, 2), P_verso_vero = round(P_verso_vero, 2)), row.names = FALSE)
