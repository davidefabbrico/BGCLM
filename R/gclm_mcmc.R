# ============================================================
# Interfaccia R al campionatore scritto in C++ (src/gclm_mcmc.cpp)
# ============================================================

#' MCMC per il modello GCLM bayesiano
#'
#' Campiona dalla posteriori di (B, delta, a, tau, pi_B) sotto
#'   y_l | B ~ iid N_p(0, Sigma(B, C)),  B Sigma + Sigma B^T + C = 0,
#' con spike-and-slab sugli off-diagonali di B e diagonale determinata da
#' B_ii = -( sum_{j != i} |B_ij| + a_i ) (stabilita' per Gershgorin).
#'
#' C e' tenuta fissa: Sigma(aB, aC) = Sigma(B, C) per ogni a > 0, quindi
#' (B, C) e' identificata solo a meno di scala e fissare C = I la determina.
#'
#' @param Y matrice n x p delle osservazioni all'equilibrio
#' @param n_iter numero totale di iterazioni
#' @param burnin iterazioni scartate (l'adattamento delle scale avviene qui)
#' @param thin passo di thinning sulle iterazioni conservate
#' @param a_pi,b_pi parametri di pi_B ~ Beta, che corregge per molteplicita'
#' @param a_tau,b_tau parametri di tau ~ InvGamma, varianza dello slab
#' @param a_shape,a_rate parametri di a_i ~ Gamma (forma e rate)
#' @param C matrice di covarianza del rumore, default l'identita'
#' @param verbose stampa l'avanzamento
#' @return lista con le catene (B, delta, a, tau, pi, loglik), le sintesi a
#'   posteriori (PIP, B_mean, B_incl, Sigma_mean) e i tassi di accettazione.
#'   B e delta sono array p x p x n_campioni.
gclm_mcmc <- function(Y,
                      n_iter = 20000, burnin = 10000, thin = 5,
                      a_pi = 1, b_pi = 1,
                      a_tau = 2, b_tau = 0.5,
                      a_shape = 2, a_rate = 2,
                      C = NULL, verbose = TRUE) {

  Y <- as.matrix(Y)
  storage.mode(Y) <- "double"
  p <- ncol(Y)
  if (is.null(C)) C <- diag(p)
  C <- as.matrix(C)

  stopifnot(nrow(Y) > 1, p > 1,
            all(is.finite(Y)),
            dim(C) == c(p, p),
            isSymmetric(unname(C), tol = 1e-10),
            n_iter > burnin, burnin >= 0, thin >= 1,
            (n_iter - burnin) %/% thin >= 1,
            a_pi > 0, b_pi > 0, a_tau > 0, b_tau > 0, a_shape > 0, a_rate > 0)
  if (min(eigen(C, symmetric = TRUE, only.values = TRUE)$values) <= 0)
    stop("C deve essere definita positiva")

  fit <- gclm_mcmc_cpp(Y, as.integer(n_iter), as.integer(burnin), as.integer(thin),
                       a_pi, b_pi, a_tau, b_tau, a_shape, a_rate, C, verbose)
  fit$n <- nrow(Y)
  fit$p <- p
  fit$call <- match.call()
  class(fit) <- "gclm"
  fit
}

#' @export
print.gclm <- function(x, ...) {
  ns <- length(x$loglik)
  cat("MCMC per GCLM bayesiano\n")
  cat(sprintf("  n = %d, p = %d, %d campioni conservati\n", x$n, x$p, ns))
  cat(sprintf("  archi attesi a posteriori: %.2f su %d\n",
              mean(apply(x$delta, 3, sum)), x$p * (x$p - 1)))
  cat(sprintf("  archi con PIP > 0.5: %d\n",
              sum(x$PIP[diag(x$p) == 0] > 0.5)))
  cat(sprintf("  accettazione: flip %.3f, rw %.3f, a %.3f\n",
              x$acc$flip, stats::median(x$acc$rw[diag(x$p) == 0]), mean(x$acc$a)))
  cat(sprintf("  log-verosimiglianza: media %.1f, ESS %.0f\n",
              mean(x$loglik), ess(x$loglik)))
  invisible(x)
}
