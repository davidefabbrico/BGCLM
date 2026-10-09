// ============================================================
// BGCLM - campionatore per grafi SEMPLICI (nessun 2-ciclo)
//
// Stesso modello di gclm_mcmc.cpp (Lyapunov, diagonale di Gershgorin, slab
// N(0, tau), margini a_i ~ Gamma) con una prior diversa sugli archi, definita
// sulle coppie non ordinate {i, j}:
//
//   stato della coppia in { nessun arco, i -> j, j -> i }
//   P(nessun arco) = 1 - rho,  P(i -> j) = P(j -> i) = rho / 2,  rho ~ Beta(a_rho, b_rho)
//
// Convenzione degli indici: l'arco i -> j e' B_ji (riga j), coerente con
// dX_j = sum_i B_ji X_i; una variazione di B_ji e' assorbita dal margine a_j.
//
// Perche' non un Bernoulli(pi) troncato ai grafi semplici. Con archi
// indipendenti Bern(pi) sulle p(p-1) coppie ordinate e il vincolo "nessun
// 2-ciclo", la costante di normalizzazione dipende da pi:
//   p(gamma | pi) = pi^E (1 - pi)^(2P - E) / (1 - pi^2)^P,   P = p(p-1)/2,
// quindi la full conditional di pi e' proporzionale a
//   Beta(pi) pi^E (1 - pi)^(2P - E) (1 - pi^2)^(-P),
// che non e' una Beta. Con la prior sulle coppie invece
//   p(gamma | rho) = (rho / 2)^E (1 - rho)^(P - E),
// e (senza augmentazione) la full conditional e' coniugata: rho | . ~ Beta(a_rho + E, b_rho + P - E).
//
// Augmentazione computazionale con eps. Con le nascite di 2-cicli vietate la
// catena non riesce a cambiare l'orientamento di piu' archi insieme e a n
// grande resta bloccata. Si campiona allora un modello ausiliario in cui la
// coppia puo' avere anche entrambi gli archi:
//   P(nessuno) = (1-rho)(1-eps),  P(i->j) = P(j->i) = rho (1-eps)/2,  P(entrambi) = eps,
// e si tengono solo i draw senza coppie reciproche (lo fa il wrapper R). Sui
// grafi semplici la prior ausiliaria e' (1-eps)^P (rho/2)^E (1-rho)^(P-E):
// il fattore (1-eps)^P non dipende da nulla, quindi la posteriori ausiliaria
// ristretta ai grafi semplici coincide con quella del modello senza 2-cicli.
// Le coppie con due archi hanno prior eps, indipendente da rho, per cui
//   rho | . ~ Beta(a_rho + n_singoli, b_rho + n_vuoti).
// Con eps = 0 si ha il campionatore senza augmentazione.
//
// Mosse per ogni coppia ordinata visitata (i, j), in ordine casuale:
//  (0) inversione: se la coppia {i, j} ha esattamente un arco, lo si sposta
//      nella direzione opposta con lo stesso valore. La mossa e' attaccata
//      alla coppia NON ordinata: dallo stato con l'arco in una direzione
//      propone l'altra e viceversa, cosi' il kernel e' reversibile da solo.
//      E' un'involuzione deterministica (jacobiano 1); numero di archi e
//      valore nello slab non cambiano, quindi prior su rho e slab si
//      semplificano. Restano verosimiglianza e margini delle due righe.
//  (a) nascita/morte dell'arco in (i, j), con i rapporti di prior sopra;
//      con eps = 0 la nascita di un 2-ciclo e' rifiutata senza risolvere Lyapunov.
//  (b) random walk sul valore, se l'arco e' attivo.
// Poi margini a_i (random walk su log) e Gibbs per tau e rho.
//
// prior_only = true rende piatta la verosimiglianza (test dei rapporti MH).
// ============================================================

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace arma;

static inline bool lyap_solve_s(const mat& B, const mat& C, mat& S) {
  if (!sylvester(S, B, B.t(), C)) return false;
  if (!S.is_finite()) return false;
  S = 0.5 * (S + S.t());
  return true;
}

static inline double gauss_loglik_s(const mat& S, const mat& Shat, double n, double p) {
  mat R;
  if (!chol(R, S)) return -datum::inf;
  const double ldet = 2.0 * accu(log(R.diag()));
  mat Ri;
  if (!inv(Ri, trimatu(R))) return -datum::inf;
  const double tr = accu((Ri * Ri.t()) % Shat);
  return -0.5 * n * (p * std::log(2.0 * datum::pi) + ldet + tr);
}

// [[Rcpp::export]]
Rcpp::List gclm_simple_cpp(const arma::mat& Y,
                           const int n_iter, const int burnin, const int thin,
                           const double a_rho, const double b_rho,
                           const double a_tau, const double b_tau,
                           const double a_shape, const double a_rate,
                           const arma::mat& C, const bool verbose,
                           const bool prior_only,
                           const double eps_bridge,
                           const arma::mat& Boff0, const arma::vec& a0) {

  const uword n = Y.n_rows, p = Y.n_cols;
  const double nd = (double) n, pd = (double) p;
  const mat Shat = Y.t() * Y / nd;
  const double a_scale = 1.0 / a_rate;
  const double Ppairs = 0.5 * pd * (pd - 1.0);
  auto loglik = [&](const mat& S) -> double {
    return prior_only ? 0.0 : gauss_loglik_s(S, Shat, nd, pd);
  };

  std::vector<uword> ei, ej;
  for (uword i = 0; i < p; ++i)
    for (uword j = 0; j < p; ++j)
      if (i != j) { ei.push_back(i); ej.push_back(j); }
  const uword K = ei.size();
  std::vector<uword> ord(K);
  for (uword k = 0; k < K; ++k) ord[k] = k;

  // stato iniziale (grafo vuoto con a = 1, oppure quello passato da R)
  mat  Boff = Boff0;
  Boff.diag().zeros();
  umat delta = conv_to<umat>::from(Boff != 0.0);
  vec  rsum = sum(abs(Boff), 1);
  vec  avec = a0;
  if (any(avec <= 0.0)) Rcpp::stop("i margini iniziali devono essere positivi");
  vec  dvec = rsum + avec;
  double tau = 0.25, rho = 0.5;

  mat Bcur = Boff;
  Bcur.diag() = -dvec;
  mat Sigcur;
  if (!lyap_solve_s(Bcur, C, Sigcur))
    Rcpp::stop("Lyapunov non risolvibile nello stato iniziale");
  double llcur = loglik(Sigcur);

  mat s_rw(p, p, fill::value(0.2)), acc_rw(p, p, fill::zeros), n_rw(p, p, fill::zeros);
  vec s_a(p, fill::value(0.3)),     acc_a(p, fill::zeros),     n_a(p, fill::zeros);
  double acc_flip = 0.0, n_flip = 0.0, acc_rev = 0.0, n_rev = 0.0, n_blocked = 0.0;

  const int n_save = (n_iter - burnin) / thin;
  cube Bsave(p, p, n_save), Dsave(p, p, n_save);
  mat  Asave(n_save, p);
  vec  tausave(n_save), rhosave(n_save), llsave(n_save);
  mat  Sigsum(p, p, fill::zeros), Bsum(p, p, fill::zeros);
  mat  Dsum(p, p, fill::zeros), BDsum(p, p, fill::zeros);

  const int report = std::max(1, n_iter / 10);
  mat Stry;

  for (int it = 1; it <= n_iter; ++it) {

    if (it % 100 == 0) Rcpp::checkUserInterrupt();

    for (uword k = K; k > 1; --k)
      std::swap(ord[k - 1], ord[(uword) std::floor(R::unif_rand() * (double) k)]);

    for (uword t = 0; t < K; ++t) {
      // coppia ordinata (i, j) del codice = voce B(i, j) = arco j -> i
      const uword i = ei[ord[t]], j = ej[ord[t]];

      // ---- (0) inversione dell'arco della coppia {i, j} ----------------
      {
        uword ro = 0, co = 0;
        bool has = true;
        if (delta(i, j) == 1 && delta(j, i) == 0) { ro = i; co = j; }
        else if (delta(j, i) == 1 && delta(i, j) == 0) { ro = j; co = i; }
        else has = false;
        if (has) {
          const double b = Boff(ro, co);
          const double a_ro = avec(ro) + std::fabs(b);   // la riga ro perde l'arco
          const double a_co = avec(co) - std::fabs(b);   // la riga co lo riceve
          n_rev += 1.0;
          if (a_co > 0.0) {
            Boff(ro, co) = 0.0; Bcur(ro, co) = 0.0; Boff(co, ro) = b; Bcur(co, ro) = b;
            const double ll2 = lyap_solve_s(Bcur, C, Stry) ? loglik(Stry) : -datum::inf;
            const double lr2 = R::dgamma(a_ro, a_shape, a_scale, 1) - R::dgamma(avec(ro), a_shape, a_scale, 1)
                             + R::dgamma(a_co, a_shape, a_scale, 1) - R::dgamma(avec(co), a_shape, a_scale, 1);
            if (std::log(R::unif_rand()) < ll2 - llcur + lr2) {
              delta(ro, co) = 0; delta(co, ro) = 1;
              rsum(ro) -= std::fabs(b); avec(ro) = a_ro;
              rsum(co) += std::fabs(b); avec(co) = a_co;
              Sigcur = Stry; llcur = ll2; acc_rev += 1.0;
            } else {
              Boff(ro, co) = b; Bcur(ro, co) = b; Boff(co, ro) = 0.0; Bcur(co, ro) = 0.0;
            }
          }
        }
      }

      // ---- (a) nascita / morte ----------------------------------------
      // Rapporti di prior fra gli stati della coppia (slab e proposta si
      // cancellano): nessuno -> un arco (rho/2)/(1-rho); un arco -> entrambi
      // eps / (rho (1-eps)/2); le morti hanno i reciproci. Con eps = 0 lo stato
      // "entrambi" ha massa nulla e la nascita e' rifiutata subito.
      double bnew, logr;
      bool blocked = false;
      if (delta(i, j) == 0) {
        bnew = R::rnorm(0.0, std::sqrt(tau));
        if (delta(j, i) == 1) {
          if (eps_bridge <= 0.0) { blocked = true; n_blocked += 1.0; logr = 0.0; }
          else logr = std::log(eps_bridge) - std::log(0.5 * rho * (1.0 - eps_bridge));
        } else {
          logr = std::log(0.5 * rho) - std::log1p(-rho);
        }
      } else {
        bnew = 0.0;
        if (delta(j, i) == 1)
          logr = std::log(0.5 * rho * (1.0 - eps_bridge)) - std::log(eps_bridge);
        else
          logr = std::log1p(-rho) - std::log(0.5 * rho);
      }
      if (!blocked) {
        double rnew = rsum(i) - std::fabs(Boff(i, j)) + std::fabs(bnew);
        double anew = dvec(i) - rnew;
        n_flip += 1.0;
        if (anew > 0.0) {
          const double bold = Boff(i, j);
          Boff(i, j) = bnew; Bcur(i, j) = bnew;
          double lltry = lyap_solve_s(Bcur, C, Stry) ? loglik(Stry) : -datum::inf;
          logr += R::dgamma(anew, a_shape, a_scale, 1) - R::dgamma(avec(i), a_shape, a_scale, 1);
          if (std::log(R::unif_rand()) < lltry - llcur + logr) {
            delta(i, j) = 1 - delta(i, j);
            rsum(i) = rnew; avec(i) = anew;
            Sigcur = Stry; llcur = lltry; acc_flip += 1.0;
          } else {
            Boff(i, j) = bold; Bcur(i, j) = bold;
          }
        }
      }

      // ---- (b) random walk sul valore (solo se l'arco e' attivo) ---------
      if (delta(i, j) == 1) {
        const double bold = Boff(i, j);
        bnew = bold + R::rnorm(0.0, s_rw(i, j));
        double rnew = rsum(i) - std::fabs(bold) + std::fabs(bnew);
        double anew = dvec(i) - rnew;
        n_rw(i, j) += 1.0;
        if (anew > 0.0) {
          Boff(i, j) = bnew; Bcur(i, j) = bnew;
          double lltry = lyap_solve_s(Bcur, C, Stry) ? loglik(Stry) : -datum::inf;
          logr = R::dnorm4(bnew, 0.0, std::sqrt(tau), 1) - R::dnorm4(bold, 0.0, std::sqrt(tau), 1)
               + R::dgamma(anew, a_shape, a_scale, 1) - R::dgamma(avec(i), a_shape, a_scale, 1);
          if (std::log(R::unif_rand()) < lltry - llcur + logr) {
            rsum(i) = rnew; avec(i) = anew;
            Sigcur = Stry; llcur = lltry; acc_rw(i, j) += 1.0;
          } else {
            Boff(i, j) = bold; Bcur(i, j) = bold;
          }
        }
      }
    }

    // ---- margini a_i -----------------------------------------------------
    for (uword i = 0; i < p; ++i) {
      const double aold = avec(i), dold = dvec(i);
      const double anew = aold * std::exp(R::rnorm(0.0, s_a(i)));
      Bcur(i, i) = -(rsum(i) + anew);
      double lltry = lyap_solve_s(Bcur, C, Stry) ? loglik(Stry) : -datum::inf;
      const double logr = R::dgamma(anew, a_shape, a_scale, 1) - R::dgamma(aold, a_shape, a_scale, 1)
                        + std::log(anew) - std::log(aold);
      n_a(i) += 1.0;
      if (std::log(R::unif_rand()) < lltry - llcur + logr) {
        avec(i) = anew; dvec(i) = rsum(i) + anew;
        Sigcur = Stry; llcur = lltry; acc_a(i) += 1.0;
      } else {
        Bcur(i, i) = -dold;
      }
    }

    // ---- Gibbs coniugati ---------------------------------------------------
    const double E = (double) accu(delta);
    tau = 1.0 / R::rgamma(a_tau + 0.5 * E, 1.0 / (b_tau + 0.5 * accu(square(Boff))));
    // le coppie con due archi hanno prior eps, che non dipende da rho
    const double n_both = 0.5 * (double) accu(delta % delta.t());
    const double n_single = E - 2.0 * n_both;
    rho = R::rbeta(a_rho + n_single, b_rho + Ppairs - n_single - n_both);

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
      acc_flip = 0.0; n_flip = 0.0; acc_rev = 0.0; n_rev = 0.0; n_blocked = 0.0;
    }

    if (it > burnin && (it - burnin) % thin == 0) {
      const int k = (it - burnin) / thin - 1;
      mat Dd = conv_to<mat>::from(delta);
      Bsave.slice(k) = Bcur; Dsave.slice(k) = Dd;
      Asave.row(k) = avec.t();
      tausave(k) = tau; rhosave(k) = rho; llsave(k) = llcur;
      Sigsum += Sigcur; Bsum += Bcur; Dsum += Dd; BDsum += Bcur % Dd;
    }

    if (verbose && it % report == 0)
      Rcpp::Rcout << "    iter " << it << " / " << n_iter
                  << "   loglik " << llcur << "   archi " << accu(delta) << std::endl;
  }

  const double ns = (double) n_save;
  return Rcpp::List::create(
    Rcpp::_["B"]          = Bsave,
    Rcpp::_["delta"]      = Dsave,
    Rcpp::_["a"]          = Asave,
    Rcpp::_["tau"]        = tausave,
    Rcpp::_["rho"]        = rhosave,
    Rcpp::_["loglik"]     = llsave,
    Rcpp::_["Sigma_mean"] = Sigsum / ns,
    Rcpp::_["PIP"]        = Dsum / ns,
    Rcpp::_["B_mean"]     = Bsum / ns,
    Rcpp::_["B_incl"]     = BDsum / clamp(Dsum, 1.0, datum::inf),
    Rcpp::_["acc"]        = Rcpp::List::create(
        Rcpp::_["flip"]    = n_flip > 0 ? acc_flip / n_flip : NA_REAL,
        Rcpp::_["rev"]     = n_rev > 0 ? acc_rev / n_rev : NA_REAL,
        Rcpp::_["blocked"] = n_blocked,
        Rcpp::_["rw"]      = acc_rw / clamp(n_rw, 1.0, datum::inf),
        Rcpp::_["a"]       = acc_a / clamp(n_a, 1.0, datum::inf)));
}
