# Test dell'augmentazione con eps: con verosimiglianza spenta,
# - catena completa: stato della coppia nessuno/un arco/entrambi con
#   probabilita' (b/(a+b))(1-eps), (a/(a+b))(1-eps), eps (in media su rho);
# - draw tenuti (senza 2-cicli): E e stati delle coppie come nella prior del
#   modello A1, cioe' beta-binomiale(P, a, b) e b/(a+b), a/(2(a+b)) ciascuna.
suppressMessages(library(BGCLM))
p <- 5; P <- p * (p - 1) / 2; pr <- which(upper.tri(diag(p)), arr.ind = TRUE)
exactE <- function(a, b) sapply(0:P, function(k) choose(P, k) * beta(a + k, b + P - k) / beta(a, b))
pair_state <- function(D) t(apply(D, 3, function(d) {
  s <- apply(pr, 1, function(ij) { f <- d[ij[2], ij[1]]; g <- d[ij[1], ij[2]]; if (f && g) 3 else if (f) 1 else if (g) 2 else 0 })
  tabulate(s + 1, 4) / nrow(pr) }))
for (cfg in list(c(1, 1, 0.05), c(2, 5, 0.1))) {
  set.seed(21)
  raw <- gclm_simple_cpp(matrix(rnorm(250), 50), 220000L, 20000L, 10L, cfg[1], cfg[2], 2, 0.5, 2, 2,
                         diag(p), FALSE, TRUE, cfg[3], matrix(0, p, p), rep(1, p))
  st <- colMeans(pair_state(raw$delta))
  keep <- apply(raw$delta, 3, function(d) !any(d == 1 & t(d) == 1))
  E <- apply(raw$delta[, , keep], 3, sum); sk <- colMeans(pair_state(raw$delta[, , keep, drop = FALSE]))
  a <- cfg[1]; b <- cfg[2]; e <- cfg[3]
  cat(sprintf("\n=== a_rho=%g b_rho=%g eps=%g: draw tenuti %.1f%% (ESS di E sui tenuti %.0f) ===\n", a, b, e, 100 * mean(keep), ess(E)))
  cat(sprintf("catena completa  nessuno/i->j/j->i/entrambi: %.3f %.3f %.3f %.3f   attesi %.3f %.3f %.3f %.3f\n",
      st[1], st[2], st[3], st[4], b / (a + b) * (1 - e), a / (2 * (a + b)) * (1 - e), a / (2 * (a + b)) * (1 - e), e))
  cat(sprintf("draw tenuti      nessuno/i->j/j->i:          %.3f %.3f %.3f        attesi %.3f %.3f %.3f\n",
      sk[1], sk[2], sk[3], b / (a + b), a / (2 * (a + b)), a / (2 * (a + b))))
  print(round(rbind(E = 0:P, tenuti = tabulate(E + 1, P + 1) / length(E), esatta = exactE(a, b)), 3))
  cat(sprintf("media di E: %.3f (esatta %.3f);  margini: media %.3f var %.3f\n", mean(E), P * a / (a + b), mean(raw$a), var(as.vector(raw$a))))
}
