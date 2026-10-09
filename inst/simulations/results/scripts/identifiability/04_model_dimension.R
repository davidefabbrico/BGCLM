# Fase 2 (complemento) - dimensione del modello di ogni equivalente trovato in
# 02_census: rango dello jacobiano di Sigma rispetto ai parametri liberi
# (archi del grafo + p margini), valutato nella deriva equivalente. Se e' piu'
# grande della dimensione del vero, il modello dell'equivalente contiene
# strettamente quello del vero: non e' equivalente, e' un modello piu' grande.
# Uso: Rscript inst/simulations/identifiability/04_model_dimension.R
source("inst/simulations/identifiability/common_ident.R")
eq <- read.csv("inst/simulations/identifiability/output/02_census_equivalents.csv")
parse_g <- function(g) { e <- do.call(rbind, lapply(strsplit(strsplit(g, " ")[[1]], "->"), as.integer))
  d <- matrix(0L, p, p); d[cbind(e[, 2], e[, 1])] <- 1L; d }
jdim <- function(B, d) { Sg <- lyapunov(B, diag(p)); A <- which(d == 1 | diag(p) == 1)
  J <- sapply(A, function(e) { E <- matrix(0, p, p); E[e] <- 1
    dS <- lyapunov(B, E %*% Sg + Sg %*% t(E)); dS[upper.tri(dS, diag = TRUE)] })
  qr(J, tol = 1e-8)$rank }
out <- list()
for (k in seq_len(nrow(eq))) {
  if (!eq$stabile[k]) next
  tr <- true_model(eq$scenario[k], eq$rep[k]); d <- parse_g(eq$grafo[k])
  M <- fiber_matrix(tr$Sigma); Z <- which(d[OFF] == 0)
  kk <- lm.fit(M[Z, , drop = FALSE], -tr$B[OFF][Z])$coefficients; kk[is.na(kk)] <- 0
  Bp <- tr$B + fiber_full(tr$Sigma, kk)
  out[[length(out) + 1]] <- data.frame(scenario = eq$scenario[k], rep = eq$rep[k],
    archi = eq$archi[k], gershgorin = eq$gershgorin[k], dim_equiv = jdim(Bp, d),
    dim_vero = jdim(tr$B, tr$delta))
}
out <- do.call(rbind, out)
write.csv(out, "inst/simulations/identifiability/output/04_model_dimension.csv", row.names = FALSE)
cat("equivalenti stabili esaminati:", nrow(out), "\n")
print(table(scenario = out$scenario, `dim equivalente - dim vero` = out$dim_equiv - out$dim_vero))
cat("dimensione del vero:", paste(unique(out$dim_vero), collapse = ", "), "\n")
