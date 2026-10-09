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
#' Due versioni del modello:
#' * `simple = FALSE`: archi indipendenti, gamma_ji ~ Bern(pi_B), pi_B ~
#'   Beta(a_pi, b_pi), 2-cicli ammessi (campionatore in src/gclm_mcmc.cpp);
#' * `simple = TRUE`: solo grafi semplici (nessuna coppia i -> j, j -> i), con
#'   prior sulle coppie non ordinate: P(nessun arco) = 1 - rho, P(i -> j) =
#'   P(j -> i) = rho / 2, rho ~ Beta(a_rho, b_rho) (src/gclm_mcmc_simple.cpp).
#'   Per il mixing si campiona un modello ausiliario in cui una coppia ha
#'   entrambi gli archi con probabilita' `eps_bridge` e si tengono i draw senza
#'   2-cicli: sui grafi semplici la prior ausiliaria differisce da quella del
#'   modello per il fattore costante (1 - eps)^P, quindi la posteriori ristretta
#'   e' esattamente quella del modello. Una mossa di inversione cambia
#'   l'orientamento di un arco in un passo.
#'
#' Convenzione: l'arco i -> j e' B_ji, cioe' `B[j, i]` (riga j).
#'
#' [restrict_simple()] da' invece la posteriori della prima versione condizionata
#' ai grafi semplici, tenendo i draw senza 2-cicli: e' un modello diverso, perche'
#' la prior implicita sulle coppie non e' quella della seconda versione.
#'
#' @param Y matrice n x p delle osservazioni all'equilibrio
#' @param n_iter numero totale di iterazioni
#' @param burnin iterazioni scartate (l'adattamento delle scale avviene qui)
#' @param thin passo di thinning sulle iterazioni conservate
#' @param a_pi,b_pi parametri di pi_B ~ Beta, che corregge per molteplicita'
#'   (solo `simple = FALSE`)
#' @param a_rho,b_rho parametri di rho ~ Beta, probabilita' che una coppia abbia
#'   un arco (solo `simple = TRUE`)
#' @param a_tau,b_tau parametri di tau ~ InvGamma, varianza dello slab
#' @param a_shape,a_rate parametri di a_i ~ Gamma (forma e rate)
#' @param C matrice di covarianza del rumore, default l'identita'
#' @param simple se TRUE usa il modello per grafi semplici (vedi sopra)
#' @param eps_bridge probabilita' a priori dello stato ausiliario "entrambi gli
#'   archi" (solo `simple = TRUE`); 0 vieta i 2-cicli anche durante il
#'   campionamento, che a n grande blocca la catena
#' @param pt_betas NULL (catena singola) oppure il vettore decrescente delle
#'   temperature inverse del parallel tempering, con primo elemento 1; si
#'   salvano i draw della catena con beta = 1
#' @param B_init deriva iniziale p x p nella classe Gershgorin (solo `simple =
#'   TRUE`); NULL parte dal grafo vuoto con margini 1
#' @param prior_only se TRUE la verosimiglianza e' piatta e la catena campiona
#'   la prior: serve solo a verificare la correttezza del campionatore
#' @param verbose stampa l'avanzamento
#' @return lista con le catene (B, delta, a, tau, pi, loglik), le sintesi a
#'   posteriori (PIP, B_mean, B_incl, Sigma_mean) e i tassi di accettazione.
#'   B e delta sono array p x p x n_campioni.
gclm_mcmc <- function(Y,
                      n_iter = 20000, burnin = 10000, thin = 5,
                      a_pi = 1, b_pi = 1, a_rho = 1, b_rho = 1,
                      a_tau = 2, b_tau = 0.5,
                      a_shape = 2, a_rate = 2,
                      C = NULL, simple = FALSE, eps_bridge = 0.05, B_init = NULL, pt_betas = NULL,
                      prior_only = FALSE, verbose = TRUE) {

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
            a_pi > 0, b_pi > 0, a_rho > 0, b_rho > 0,
            a_tau > 0, b_tau > 0, a_shape > 0, a_rate > 0,
            is.logical(simple), length(simple) == 1, is.logical(prior_only))
  if (min(eigen(C, symmetric = TRUE, only.values = TRUE)$values) <= 0)
    stop("C deve essere definita positiva")
  stopifnot(eps_bridge >= 0, eps_bridge < 1)
  Boff0 <- matrix(0, p, p); a0 <- rep(1, p)
  if (!is.null(B_init)) {
    stopifnot(dim(B_init) == c(p, p))
    Boff0 <- B_init; diag(Boff0) <- 0
    a0 <- -diag(B_init) - rowSums(abs(Boff0))
    if (any(a0 <= 0)) stop("B_init deve avere diagonale di Gershgorin (a_i > 0)")
  }

  if (!is.null(pt_betas))
    stopifnot(abs(pt_betas[1] - 1) < 1e-12, all(diff(pt_betas) < 0), all(pt_betas > 0))
  fit <- if (!is.null(pt_betas)) {
    out <- gclm_simple_pt_cpp(Y, as.integer(n_iter), as.integer(burnin), as.integer(thin),
             if (simple) a_rho else a_pi, if (simple) b_rho else b_pi,
             a_tau, b_tau, a_shape, a_rate, C, verbose, prior_only,
             eps_bridge, Boff0, a0, as.numeric(pt_betas), !simple)
    if (!simple) { out$pi <- out$rho; out$rho <- NULL }
    out$C <- C
    if (!simple) out <- .summarise_draws(out, rep(TRUE, length(out$loglik)))
    out
  } else if (simple)
    gclm_simple_cpp(Y, as.integer(n_iter), as.integer(burnin), as.integer(thin),
                    a_rho, b_rho, a_tau, b_tau, a_shape, a_rate, C, verbose, prior_only,
                    eps_bridge, Boff0, a0)
  else
    gclm_mcmc_cpp(Y, as.integer(n_iter), as.integer(burnin), as.integer(thin),
                  a_pi, b_pi, a_tau, b_tau, a_shape, a_rate, C, verbose, prior_only)
  fit$n <- nrow(Y)
  fit$p <- p
  fit$C <- C
  fit$simple <- FALSE
  fit$call <- match.call()
  class(fit) <- "gclm"
  if (simple) {
    fit$eps_bridge <- eps_bridge
    fit <- .subset_simple(fit, "prior sulle coppie", min_draws = 100)
  }
  fit
}

# tiene i draw senza coppie reciproche e ricalcola le sintesi su quei draw
.subset_simple <- function(fit, method, min_draws) {
  keep <- apply(fit$delta, 3, function(d) !any(d == 1 & t(d) == 1))
  if (!any(keep)) stop("nessun draw senza 2-cicli: aumentare n_iter")
  if (sum(keep) < min_draws)
    warning(sprintf("solo %d draw senza 2-cicli: le sintesi sono poco precise, aumentare n_iter",
                    sum(keep)))
  out <- .summarise_draws(fit, keep)
  out$simple        <- TRUE
  out$simple_method <- method
  out$keep          <- keep
  out$n_kept        <- sum(keep)
  out$frac_kept     <- mean(keep)
  out
}

# sintesi a posteriori (MPPI, medie, Sigma media) sui draw indicati da keep
.summarise_draws <- function(fit, keep) {
  out <- fit
  out$B      <- fit$B[, , keep, drop = FALSE]
  out$delta  <- fit$delta[, , keep, drop = FALSE]
  out$a      <- fit$a[keep, , drop = FALSE]
  out$tau    <- fit$tau[keep]
  if (!is.null(fit$pi))  out$pi  <- fit$pi[keep]
  if (!is.null(fit$rho)) out$rho <- fit$rho[keep]
  out$loglik <- fit$loglik[keep]
  out$PIP    <- apply(out$delta, c(1, 2), mean)
  out$B_mean <- apply(out$B, c(1, 2), mean)
  out$B_incl <- apply(out$B * out$delta, c(1, 2), sum) /
                pmax(apply(out$delta, c(1, 2), sum), 1)
  out$Sigma_mean <- Reduce(`+`, lapply(seq_len(sum(keep)), function(k)
    lyapunov(out$B[, , k], fit$C))) / sum(keep)
  out
}

#' Restringe un fit ai grafi senza 2-cicli
#'
#' Tiene i draw il cui grafo non contiene coppie reciproche i -> j, j -> i e
#' ricalcola le sintesi a posteriori su quei draw. Permette di avere dalla stessa
#' catena sia la posteriori libera sia quella ristretta.
#'
#' @param fit output di [gclm_mcmc()] con `simple = FALSE`
#' @param min_draws sotto questo numero di draw conservati viene dato un avviso
#' @return un oggetto "gclm" con le catene filtrate, `simple = TRUE`, il numero
#'   (`n_kept`) e la frazione (`frac_kept`) di draw conservati
restrict_simple <- function(fit, min_draws = 200) {
  stopifnot(inherits(fit, "gclm"), !isTRUE(fit$simple))
  .subset_simple(fit, "filtro a posteriori", min_draws)
}

#' Grafo selezionato da un fit
#'
#' Regola per freccia: arco i -> j presente se MPPI_ji > thr. Regola per coppia:
#' nel modello A1 gli stati di una coppia (nessun arco, i -> j, j -> i) si
#' escludono a vicenda, quindi MPPI_ij + MPPI_ji e' la probabilita' a posteriori
#' che la coppia abbia un arco; la coppia e' selezionata se questa supera thr e
#' l'arco e' orientato nel verso con MPPI maggiore. Con il verso incerto la regola
#' per freccia sottoseleziona (due MPPI di 0.4 non danno nessun arco).
#'
#' @param fit output di [gclm_mcmc()]
#' @param rule "auto" (coppia per i fit ristretti ai grafi semplici, freccia
#'   altrimenti), "pair" o "arrow"
#' @param thr soglia
#' @return matrice 0/1 p x p, con `[j, i] = 1` per l'arco i -> j
select_graph <- function(fit, rule = c("auto", "pair", "arrow"), thr = 0.5) {
  rule <- match.arg(rule)
  if (rule == "auto") rule <- if (isTRUE(fit$simple)) "pair" else "arrow"
  M <- fit$PIP; diag(M) <- 0
  if (rule == "arrow") return((M > thr) * 1L)
  ((M + t(M) > thr) & (M > t(M) | (M == t(M) & lower.tri(M)))) * 1L
}

#' @export
print.gclm <- function(x, ...) {
  ns <- length(x$loglik)
  cat("MCMC per GCLM bayesiano",
      if (isTRUE(x$simple)) sprintf("(grafi semplici, %s)", x$simple_method), "\n")
  cat(sprintf("  n = %d, p = %d, %d campioni conservati", x$n, x$p, ns))
  if (!is.null(x$frac_kept)) cat(sprintf(" (%.0f%% dei draw)", 100 * x$frac_kept))
  cat("\n")
  cat(sprintf("  archi attesi a posteriori: %.2f su %d\n",
              mean(apply(x$delta, 3, sum)), x$p * (x$p - 1)))
  cat(sprintf("  archi selezionati (regola per %s): %d\n",
              if (isTRUE(x$simple)) "coppia" else "freccia", sum(select_graph(x))))
  cat(sprintf("  accettazione: nascita/morte %.3f, inversione %.3f, rw %.3f, a %.3f\n",
              x$acc$flip, x$acc$rev, stats::median(x$acc$rw[diag(x$p) == 0]), mean(x$acc$a)))
  cat(sprintf("  log-verosimiglianza: media %.1f, ESS %.0f\n",
              mean(x$loglik), ess(x$loglik)))
  invisible(x)
}
