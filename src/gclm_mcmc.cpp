// ============================================================
// Bayesian Graphical Continuous Lyapunov Model (GCLM) - MCMC
//
//   B Sigma + Sigma B^T + C = 0        ->  Sigma(B, C)
//   y_l | B ~ iid N_p(0, Sigma(B, C))
//
//   B_ji | delta_ji ~ delta_ji N(0, tau) + (1 - delta_ji) delta_0
//   delta_ji ~ Bern(pi_B)
//   B_ii = -( sum_{j != i} |B_ij| + a_i ),  a_i > 0   (stabilita' per Gershgorin)
//
// Parametrizzazione: si campiona in (B_off, d) con d_i = -B_ii, tenendo
// a_i = d_i - sum_j |B_ij| > 0 come vincolo di supporto. E' lo stesso modello
// (mappa a -> d uno-a-uno a d fissato, jacobiano 1), ma in coordinate (B_off, a)
// aggiungere un arco sposterebbe anche B_ii: la mossa perturba tutta la riga di
// Sigma e viene quasi sempre rifiutata. In (B_off, d) le mosse sugli archi sono
// locali, B_ii resta fermo, e la catena mescola.
//
// Oltre a nascita/morte e random walk c'e' una mossa che inverte un arco
// (j -> i diventa i -> j) in un solo passo. Senza di essa l'unico modo di
// cambiare orientamento e' passare per il 2-ciclo i <-> j, che a n grande e'
// raro: la catena resterebbe sull'orientamento scelto all'inizio.
//
// prior_only = true rende piatta la verosimiglianza: la catena deve allora
// restituire la prior, ed e' il test di correttezza dei rapporti MH.
//
// L'RNG e' quello di R: set.seed() dal lato R rende la catena riproducibile.
// ============================================================

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace arma;

// Lyapunov via Bartels-Stewart (LAPACK dtrsyl): sylvester risolve A X + X B + C = 0,
// quindi sylvester(B, B^T, C) e' esattamente l'equazione che serve. O(p^3) e stabile
// anche quando B e' difettiva, caso in cui la via spettrale non ha base completa.
static inline bool lyap_solve(const mat& B, const mat& C, mat& S) {
  if (!sylvester(S, B, B.t(), C)) return false;
  if (!S.is_finite()) return false;
  S = 0.5 * (S + S.t());
  return true;
}

static inline double gauss_loglik(const mat& S, const mat& Shat,
                                  double n, double p) {
  mat R;
  if (!chol(R, S)) return -datum::inf;                 // S = R^T R, R triangolare sup.
  const double ldet = 2.0 * accu(log(R.diag()));
  mat Ri;
  if (!inv(Ri, trimatu(R))) return -datum::inf;
  const double tr = accu((Ri * Ri.t()) % Shat);        // tr(S^{-1} Shat)
  return -0.5 * n * (p * std::log(2.0 * datum::pi) + ldet + tr);
}

// Risolve B Sigma + Sigma B^T + C = 0. Esposta a R: serve per costruire Sigma
// da un B noto (generazione dei dati) senza rifare il solver in R.
// [[Rcpp::export]]
Rcpp::NumericMatrix lyapunov_cpp(const arma::mat& B, const arma::mat& C) {
  mat S;
  if (!lyap_solve(B, C, S))
    Rcpp::stop("equazione di Lyapunov non risolvibile: B e' stabile?");
  return Rcpp::wrap(S);
}

// [[Rcpp::export]]
Rcpp::List gclm_mcmc_cpp(const arma::mat& Y,
                         const int n_iter, const int burnin, const int thin,
                         const double a_pi, const double b_pi,
                         const double a_tau, const double b_tau,
                         const double a_shape, const double a_rate,
                         const arma::mat& C, const bool verbose,
                         const bool prior_only = false) {

  const uword n = Y.n_rows, p = Y.n_cols;
  const double nd = (double) n, pd = (double) p;
  const mat Shat = Y.t() * Y / nd;
  const double a_scale = 1.0 / a_rate;                 // R::dgamma/rgamma usano la scala
  auto loglik = [&](const mat& S) -> double {
    return prior_only ? 0.0 : gauss_loglik(S, Shat, nd, pd);
  };

  // archi (i, j) con i != j
  std::vector<uword> ei, ej;
  for (uword i = 0; i < p; ++i)
    for (uword j = 0; j < p; ++j)
      if (i != j) { ei.push_back(i); ej.push_back(j); }
  const uword K = ei.size();
  std::vector<uword> ord(K);
  for (uword k = 0; k < K; ++k) ord[k] = k;

  mat  Boff(p, p, fill::zeros);
  umat delta(p, p, fill::zeros);
  vec  rsum(p, fill::zeros);                           // sum_{j != i} |B_ij|
  vec  avec(p, fill::ones);
  vec  dvec = rsum + avec;
  double tau = 0.25, pi_B = 0.5;

  mat Bcur = Boff;
  Bcur.diag() = -dvec;
  mat Sigcur;
  if (!lyap_solve(Bcur, C, Sigcur))
    Rcpp::stop("Lyapunov non risolvibile nello stato iniziale");
  double llcur = loglik(Sigcur);

  mat s_rw(p, p, fill::value(0.2)), acc_rw(p, p, fill::zeros), n_rw(p, p, fill::zeros);
  vec s_a(p, fill::value(0.3)),     acc_a(p, fill::zeros),     n_a(p, fill::zeros);
  double acc_flip = 0.0, n_flip = 0.0, acc_rev = 0.0, n_rev = 0.0;

  const int n_save = (n_iter - burnin) / thin;
  cube Bsave(p, p, n_save), Dsave(p, p, n_save);
  mat  Asave(n_save, p);
  vec  tausave(n_save), pisave(n_save), llsave(n_save);
  mat  Sigsum(p, p, fill::zeros), Bsum(p, p, fill::zeros);
  mat  Dsum(p, p, fill::zeros), BDsum(p, p, fill::zeros);

  const int report = std::max(1, n_iter / 10);
  mat Stry;

  for (int it = 1; it <= n_iter; ++it) {

    if (it % 100 == 0) Rcpp::checkUserInterrupt();

    // ordine di scansione degli archi casuale (Fisher-Yates con l'RNG di R)
    for (uword k = K; k > 1; --k)
      std::swap(ord[k - 1], ord[(uword) std::floor(R::unif_rand() * (double) k)]);

    for (uword t = 0; t < K; ++t) {
      const uword i = ei[ord[t]], j = ej[ord[t]];

      // ---- (0) inversione di un arco della coppia {i, j} --------------
      // La stessa mossa propone entrambe le direzioni, cosi' il kernel e'
      // reversibile da solo: se fosse attaccata solo alla coppia ordinata (i, j)
      // la catena non conserverebbe la posteriori. Due varianti: scambio con la
      // stessa ampiezza (involuzione, jacobiano 1) e scambio con ampiezza nuova
      // dallo slab (slab e proposta si cancellano). Le diagonali restano ferme e
      // i margini delle due righe assorbono la variazione.
      for (int mv = 0; mv < 2; ++mv) {
        uword ro, co;
        if (delta(i, j) == 1 && delta(j, i) == 0) { ro = i; co = j; }
        else if (delta(j, i) == 1 && delta(i, j) == 0) { ro = j; co = i; }
        else break;
        const double b_old = Boff(ro, co);
        const double b_new = (mv == 0) ? b_old : R::rnorm(0.0, std::sqrt(tau));
        const double a_ro = avec(ro) + std::fabs(b_old);
        const double a_co = avec(co) - std::fabs(b_new);
        n_rev += 1.0;
        if (a_co <= 0.0) continue;
        Boff(ro, co) = 0.0; Bcur(ro, co) = 0.0; Boff(co, ro) = b_new; Bcur(co, ro) = b_new;
        const double ll2 = lyap_solve(Bcur, C, Stry) ? loglik(Stry) : -datum::inf;
        const double lr2 = R::dgamma(a_ro, a_shape, a_scale, 1) - R::dgamma(avec(ro), a_shape, a_scale, 1)
                         + R::dgamma(a_co, a_shape, a_scale, 1) - R::dgamma(avec(co), a_shape, a_scale, 1);
        if (std::log(R::unif_rand()) < ll2 - llcur + lr2) {
          delta(ro, co) = 0; delta(co, ro) = 1;
          rsum(ro) -= std::fabs(b_old); avec(ro) = a_ro;
          rsum(co) += std::fabs(b_new); avec(co) = a_co;
          Sigcur = Stry; llcur = ll2; acc_rev += 1.0;
        } else {
          Boff(ro, co) = b_old; Bcur(ro, co) = b_old; Boff(co, ro) = 0.0; Bcur(co, ro) = 0.0;
        }
      }

      // ---- (a) birth / death ----------------------------------------
      // Proponendo dallo slab N(0, tau) i termini di slab e di proposta si
      // cancellano: resta lik x pi/(1-pi) x rapporto del prior su a_i.
      double bnew, logr;
      if (delta(i, j) == 0) {
        bnew = R::rnorm(0.0, std::sqrt(tau));
        logr =  std::log(pi_B) - std::log1p(-pi_B);
      } else {
        bnew = 0.0;
        logr = -std::log(pi_B) + std::log1p(-pi_B);
      }
      double rnew = rsum(i) - std::fabs(Boff(i, j)) + std::fabs(bnew);
      double anew = dvec(i) - rnew;
      n_flip += 1.0;
      if (anew > 0.0) {
        const double bold = Boff(i, j);
        Boff(i, j) = bnew; Bcur(i, j) = bnew;
        double lltry = lyap_solve(Bcur, C, Stry)
                       ? loglik(Stry) : -datum::inf;
        logr += R::dgamma(anew, a_shape, a_scale, 1)
              - R::dgamma(avec(i), a_shape, a_scale, 1);
        if (std::log(R::unif_rand()) < lltry - llcur + logr) {
          delta(i, j) = 1 - delta(i, j);
          rsum(i) = rnew; avec(i) = anew;
          Sigcur = Stry; llcur = lltry; acc_flip += 1.0;
        } else {
          Boff(i, j) = bold; Bcur(i, j) = bold;
        }
      }

      // ---- (b) random walk sull'ampiezza (solo se l'arco e' attivo) ----
      if (delta(i, j) == 1) {
        const double bold = Boff(i, j);
        bnew = bold + R::rnorm(0.0, s_rw(i, j));
        rnew = rsum(i) - std::fabs(bold) + std::fabs(bnew);
        anew = dvec(i) - rnew;
        n_rw(i, j) += 1.0;
        if (anew > 0.0) {
          Boff(i, j) = bnew; Bcur(i, j) = bnew;
          double lltry = lyap_solve(Bcur, C, Stry)
                         ? loglik(Stry) : -datum::inf;
          logr = R::dnorm4(bnew, 0.0, std::sqrt(tau), 1)
               - R::dnorm4(bold, 0.0, std::sqrt(tau), 1)
               + R::dgamma(anew, a_shape, a_scale, 1)
               - R::dgamma(avec(i), a_shape, a_scale, 1);
          if (std::log(R::unif_rand()) < lltry - llcur + logr) {
            rsum(i) = rnew; avec(i) = anew;
            Sigcur = Stry; llcur = lltry; acc_rw(i, j) += 1.0;
          } else {
            Boff(i, j) = bold; Bcur(i, j) = bold;
          }
        }
      }
    }

    // ---- margini di stabilita' a_i (random walk su log, con jacobiano) ----
    for (uword i = 0; i < p; ++i) {
      const double aold = avec(i), dold = dvec(i);
      const double anew = aold * std::exp(R::rnorm(0.0, s_a(i)));
      Bcur(i, i) = -(rsum(i) + anew);
      double lltry = lyap_solve(Bcur, C, Stry)
                     ? loglik(Stry) : -datum::inf;
      const double logr = R::dgamma(anew, a_shape, a_scale, 1)
                        - R::dgamma(aold, a_shape, a_scale, 1)
                        + std::log(anew) - std::log(aold);
      n_a(i) += 1.0;
      if (std::log(R::unif_rand()) < lltry - llcur + logr) {
        avec(i) = anew; dvec(i) = rsum(i) + anew;
        Sigcur = Stry; llcur = lltry; acc_a(i) += 1.0;
      } else {
        Bcur(i, i) = -dold;
      }
    }

    // ---- Gibbs coniugati: non richiedono Lyapunov -------------------
    const double K1 = (double) accu(delta);
    tau  = 1.0 / R::rgamma(a_tau + 0.5 * K1,
                           1.0 / (b_tau + 0.5 * accu(square(Boff))));
    pi_B = R::rbeta(a_pi + K1, b_pi + (double) K - K1);

    // ---- adattamento delle scale (solo in burn-in) ------------------
    if (it <= burnin && it % 50 == 0) {
      const double d = std::min(0.1, 1.0 / std::sqrt((double) it / 50.0));
      for (uword t = 0; t < K; ++t) {
        const uword i = ei[t], j = ej[t];
        const double r = n_rw(i, j) > 0.0 ? acc_rw(i, j) / n_rw(i, j) : 0.44;
        s_rw(i, j) = std::min(10.0, std::max(1e-4, s_rw(i, j) * std::exp(d * (r - 0.44))));
      }
      for (uword i = 0; i < p; ++i) {
        const double r = n_a(i) > 0.0 ? acc_a(i) / n_a(i) : 0.44;
        s_a(i) = std::min(10.0, std::max(1e-4, s_a(i) * std::exp(d * (r - 0.44))));
      }
      acc_rw.zeros(); n_rw.zeros(); acc_a.zeros(); n_a.zeros();
    }
    if (it == burnin) {
      acc_rw.zeros(); n_rw.zeros(); acc_a.zeros(); n_a.zeros();
      acc_flip = 0.0; n_flip = 0.0; acc_rev = 0.0; n_rev = 0.0;
    }

    if (it > burnin && (it - burnin) % thin == 0) {
      const int k = (it - burnin) / thin - 1;
      mat Dd = conv_to<mat>::from(delta);
      Bsave.slice(k) = Bcur; Dsave.slice(k) = Dd;
      Asave.row(k)   = avec.t();
      tausave(k) = tau; pisave(k) = pi_B; llsave(k) = llcur;
      Sigsum += Sigcur; Bsum += Bcur; Dsum += Dd; BDsum += Bcur % Dd;
    }

    if (verbose && it % report == 0)
      Rcpp::Rcout << "    iter " << it << " / " << n_iter
                  << "   loglik " << llcur
                  << "   archi " << accu(delta) << std::endl;
  }

  const double ns = (double) n_save;
  return Rcpp::List::create(
    Rcpp::_["B"]          = Bsave,
    Rcpp::_["delta"]      = Dsave,
    Rcpp::_["a"]          = Asave,
    Rcpp::_["tau"]        = tausave,
    Rcpp::_["pi"]         = pisave,
    Rcpp::_["loglik"]     = llsave,
    Rcpp::_["Sigma_mean"] = Sigsum / ns,
    Rcpp::_["PIP"]        = Dsum / ns,
    Rcpp::_["B_mean"]     = Bsum / ns,
    // media condizionata all'inclusione: la media marginale e' schiacciata verso
    // 0 dal peso dello spike, quindi non va letta come stima dell'ampiezza
    Rcpp::_["B_incl"]     = BDsum / clamp(Dsum, 1.0, datum::inf),
    Rcpp::_["acc"]        = Rcpp::List::create(
        Rcpp::_["flip"] = n_flip > 0 ? acc_flip / n_flip : NA_REAL,
        Rcpp::_["rev"]  = n_rev > 0 ? acc_rev / n_rev : NA_REAL,
        Rcpp::_["rw"]   = acc_rw / clamp(n_rw, 1.0, datum::inf),
        Rcpp::_["a"]    = acc_a / clamp(n_a, 1.0, datum::inf)));
}
