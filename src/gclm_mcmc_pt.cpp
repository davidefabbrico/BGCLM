// ============================================================
// BGCLM - parallel tempering per il modello A1 (grafi semplici) e per il
// modello libero (free_model = true)
//
// Stesso modello e stesse mosse di gclm_mcmc_simple.cpp (prior sulle coppie,
// augmentazione eps, inversione di un arco, nascita/morte, random walk,
// margini, Gibbs per tau e rho). Si fanno girare K catene con bersaglio
//   pi_k(theta) proporzionale a L(B)^beta_k p(theta),   1 = beta_1 > ... > beta_K > 0,
// e dopo ogni sweep si propone lo scambio degli stati fra livelli adiacenti,
// accettato con probabilita' min{1, exp((beta_k - beta_{k+1}) (l_{k+1} - l_k))},
// dove l e' la log-verosimiglianza non temperata. Ogni kernel lascia invariata
// la propria pi_k e lo scambio lascia invariato il prodotto delle pi_k, quindi
// la catena fredda (beta = 1) campiona la posteriori. Si salvano solo i suoi
// draw; il filtro sui grafi semplici lo fa il wrapper R.
//
// Con free_model = true la prior sugli archi e' quella del modello libero:
// indicatori indipendenti Bern(pi) sulle p(p-1) coppie ordinate, pi ~
// Beta(a, b), 2-cicli ammessi (a, b ed eps sono passati in a_rho, b_rho;
// eps e' ignorato). Nascita pi/(1-pi), morte il reciproco, inversione con
// rapporto di prior 1, pi | . ~ Beta(a + E, b + p(p-1) - E).
//
// Gli stati (parametri del modello) passano da un livello all'altro con gli
// scambi; le scale adattive del random walk restano legate al livello, perche'
// la scala giusta dipende dalla temperatura.
// ============================================================

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace arma;

static inline bool lyap_solve_t(const mat& B, const mat& C, mat& S) {
  if (!sylvester(S, B, B.t(), C)) return false;
  if (!S.is_finite()) return false;
  S = 0.5 * (S + S.t());
  return true;
}

static inline double gauss_loglik_t(const mat& S, const mat& Shat, double n, double p) {
  mat R;
  if (!chol(R, S)) return -datum::inf;
  const double ldet = 2.0 * accu(log(R.diag()));
  mat Ri;
  if (!inv(Ri, trimatu(R))) return -datum::inf;
  const double tr = accu((Ri * Ri.t()) % Shat);
  return -0.5 * n * (p * std::log(2.0 * datum::pi) + ldet + tr);
}

struct PTState {
  mat Boff, Bcur, Sigcur;
  umat delta;
  vec rsum, avec, dvec;
  double llcur, tau, rho;
};

struct PTTune {
  mat s_rw, acc_rw, n_rw;
  vec s_a, acc_a, n_a;
  double acc_flip = 0.0, n_flip = 0.0, acc_rev = 0.0, n_rev = 0.0;
};

// [[Rcpp::export]]
Rcpp::List gclm_simple_pt_cpp(const arma::mat& Y,
                              const int n_iter, const int burnin, const int thin,
                              const double a_rho, const double b_rho,
                              const double a_tau, const double b_tau,
                              const double a_shape, const double a_rate,
                              const arma::mat& C, const bool verbose,
                              const bool prior_only, const double eps_bridge,
                              const arma::mat& Boff0, const arma::vec& a0,
                              const arma::vec& betas, const bool free_model) {

  const uword n = Y.n_rows, p = Y.n_cols, KT = betas.n_elem;
  const double nd = (double) n, pd = (double) p;
  const mat Shat = Y.t() * Y / nd;
  const double a_scale = 1.0 / a_rate;
  const double Ppairs = 0.5 * pd * (pd - 1.0);
  if (KT < 1 || std::fabs(betas(0) - 1.0) > 1e-12) Rcpp::stop("betas deve iniziare con 1");
  auto loglik = [&](const mat& S) -> double {
    return prior_only ? 0.0 : gauss_loglik_t(S, Shat, nd, pd);
  };

  std::vector<uword> ei, ej;
  for (uword i = 0; i < p; ++i)
    for (uword j = 0; j < p; ++j)
      if (i != j) { ei.push_back(i); ej.push_back(j); }
  const uword K = ei.size();
  std::vector<uword> ord(K);
  for (uword k = 0; k < K; ++k) ord[k] = k;

  // stati e scale per livello; perm[k] = indice dello stato al livello k
  std::vector<PTState> st(KT);
  std::vector<PTTune>  tu(KT);
  std::vector<uword>   perm(KT);
  for (uword k = 0; k < KT; ++k) {
    PTState& s = st[k];
    s.Boff = Boff0; s.Boff.diag().zeros();
    s.delta = conv_to<umat>::from(s.Boff != 0.0);
    s.rsum = sum(abs(s.Boff), 1); s.avec = a0; s.dvec = s.rsum + s.avec;
    if (any(s.avec <= 0.0)) Rcpp::stop("i margini iniziali devono essere positivi");
    s.Bcur = s.Boff; s.Bcur.diag() = -s.dvec;
    if (!lyap_solve_t(s.Bcur, C, s.Sigcur)) Rcpp::stop("Lyapunov non risolvibile nello stato iniziale");
    s.llcur = loglik(s.Sigcur); s.tau = 0.25; s.rho = 0.5;
    tu[k].s_rw = mat(p, p, fill::value(0.2)); tu[k].acc_rw = mat(p, p, fill::zeros); tu[k].n_rw = mat(p, p, fill::zeros);
    tu[k].s_a = vec(p, fill::value(0.3));     tu[k].acc_a = vec(p, fill::zeros);     tu[k].n_a = vec(p, fill::zeros);
    perm[k] = k;
  }
  vec swap_acc(KT > 1 ? KT - 1 : 1, fill::zeros), swap_n(KT > 1 ? KT - 1 : 1, fill::zeros);

  const int n_save = (n_iter - burnin) / thin;
  cube Bsave(p, p, n_save), Dsave(p, p, n_save);
  mat  Asave(n_save, p);
  vec  tausave(n_save), rhosave(n_save), llsave(n_save);
  const int report = std::max(1, n_iter / 10);
  mat Stry;

  for (int it = 1; it <= n_iter; ++it) {
    if (it % 100 == 0) Rcpp::checkUserInterrupt();

    for (uword lev = 0; lev < KT; ++lev) {
      PTState& s = st[perm[lev]];
      PTTune&  T = tu[lev];
      const double beta = betas(lev);

      for (uword k = K; k > 1; --k)
        std::swap(ord[k - 1], ord[(uword) std::floor(R::unif_rand() * (double) k)]);

      for (uword t = 0; t < K; ++t) {
        const uword i = ei[ord[t]], j = ej[ord[t]];

        // (0) inversione dell'arco della coppia {i, j}, stesso valore
        {
          uword ro = 0, co = 0; bool has = true;
          if (s.delta(i, j) == 1 && s.delta(j, i) == 0) { ro = i; co = j; }
          else if (s.delta(j, i) == 1 && s.delta(i, j) == 0) { ro = j; co = i; }
          else has = false;
          if (has) {
            const double b = s.Boff(ro, co);
            const double a_ro = s.avec(ro) + std::fabs(b), a_co = s.avec(co) - std::fabs(b);
            T.n_rev += 1.0;
            if (a_co > 0.0) {
              s.Boff(ro, co) = 0.0; s.Bcur(ro, co) = 0.0; s.Boff(co, ro) = b; s.Bcur(co, ro) = b;
              const double ll2 = lyap_solve_t(s.Bcur, C, Stry) ? loglik(Stry) : -datum::inf;
              const double lr2 = R::dgamma(a_ro, a_shape, a_scale, 1) - R::dgamma(s.avec(ro), a_shape, a_scale, 1)
                               + R::dgamma(a_co, a_shape, a_scale, 1) - R::dgamma(s.avec(co), a_shape, a_scale, 1);
              if (std::log(R::unif_rand()) < beta * (ll2 - s.llcur) + lr2) {
                s.delta(ro, co) = 0; s.delta(co, ro) = 1;
                s.rsum(ro) -= std::fabs(b); s.avec(ro) = a_ro;
                s.rsum(co) += std::fabs(b); s.avec(co) = a_co;
                s.Sigcur = Stry; s.llcur = ll2; T.acc_rev += 1.0;
              } else {
                s.Boff(ro, co) = b; s.Bcur(ro, co) = b; s.Boff(co, ro) = 0.0; s.Bcur(co, ro) = 0.0;
              }
            }
          }
        }

        // (a) nascita / morte con la prior aumentata
        double bnew, logr = 0.0; bool blocked = false;
        if (free_model) {
          bnew = (s.delta(i, j) == 0) ? R::rnorm(0.0, std::sqrt(s.tau)) : 0.0;
          logr = (s.delta(i, j) == 0) ? std::log(s.rho) - std::log1p(-s.rho)
                                      : std::log1p(-s.rho) - std::log(s.rho);
        } else if (s.delta(i, j) == 0) {
          bnew = R::rnorm(0.0, std::sqrt(s.tau));
          if (s.delta(j, i) == 1) {
            if (eps_bridge <= 0.0) blocked = true;
            else logr = std::log(eps_bridge) - std::log(0.5 * s.rho * (1.0 - eps_bridge));
          } else logr = std::log(0.5 * s.rho) - std::log1p(-s.rho);
        } else {
          bnew = 0.0;
          if (s.delta(j, i) == 1) logr = std::log(0.5 * s.rho * (1.0 - eps_bridge)) - std::log(eps_bridge);
          else                    logr = std::log1p(-s.rho) - std::log(0.5 * s.rho);
        }
        if (!blocked) {
          const double rnew = s.rsum(i) - std::fabs(s.Boff(i, j)) + std::fabs(bnew);
          const double anew = s.dvec(i) - rnew;
          T.n_flip += 1.0;
          if (anew > 0.0) {
            const double bold = s.Boff(i, j);
            s.Boff(i, j) = bnew; s.Bcur(i, j) = bnew;
            const double lltry = lyap_solve_t(s.Bcur, C, Stry) ? loglik(Stry) : -datum::inf;
            logr += R::dgamma(anew, a_shape, a_scale, 1) - R::dgamma(s.avec(i), a_shape, a_scale, 1);
            if (std::log(R::unif_rand()) < beta * (lltry - s.llcur) + logr) {
              s.delta(i, j) = 1 - s.delta(i, j);
              s.rsum(i) = rnew; s.avec(i) = anew;
              s.Sigcur = Stry; s.llcur = lltry; T.acc_flip += 1.0;
            } else { s.Boff(i, j) = bold; s.Bcur(i, j) = bold; }
          }
        }

        // (b) random walk sul valore
        if (s.delta(i, j) == 1) {
          const double bold = s.Boff(i, j);
          bnew = bold + R::rnorm(0.0, T.s_rw(i, j));
          const double rnew = s.rsum(i) - std::fabs(bold) + std::fabs(bnew);
          const double anew = s.dvec(i) - rnew;
          T.n_rw(i, j) += 1.0;
          if (anew > 0.0) {
            s.Boff(i, j) = bnew; s.Bcur(i, j) = bnew;
            const double lltry = lyap_solve_t(s.Bcur, C, Stry) ? loglik(Stry) : -datum::inf;
            const double lr = R::dnorm4(bnew, 0.0, std::sqrt(s.tau), 1) - R::dnorm4(bold, 0.0, std::sqrt(s.tau), 1)
                            + R::dgamma(anew, a_shape, a_scale, 1) - R::dgamma(s.avec(i), a_shape, a_scale, 1);
            if (std::log(R::unif_rand()) < beta * (lltry - s.llcur) + lr) {
              s.rsum(i) = rnew; s.avec(i) = anew;
              s.Sigcur = Stry; s.llcur = lltry; T.acc_rw(i, j) += 1.0;
            } else { s.Boff(i, j) = bold; s.Bcur(i, j) = bold; }
          }
        }
      }

      // margini a_i
      for (uword i = 0; i < p; ++i) {
        const double aold = s.avec(i), dold = s.dvec(i);
        const double anew = aold * std::exp(R::rnorm(0.0, T.s_a(i)));
        s.Bcur(i, i) = -(s.rsum(i) + anew);
        const double lltry = lyap_solve_t(s.Bcur, C, Stry) ? loglik(Stry) : -datum::inf;
        const double lr = R::dgamma(anew, a_shape, a_scale, 1) - R::dgamma(aold, a_shape, a_scale, 1)
                        + std::log(anew) - std::log(aold);
        T.n_a(i) += 1.0;
        if (std::log(R::unif_rand()) < beta * (lltry - s.llcur) + lr) {
          s.avec(i) = anew; s.dvec(i) = s.rsum(i) + anew;
          s.Sigcur = Stry; s.llcur = lltry; T.acc_a(i) += 1.0;
        } else s.Bcur(i, i) = -dold;
      }

      // Gibbs (non dipendono dalla verosimiglianza, quindi nemmeno da beta)
      const double E = (double) accu(s.delta);
      s.tau = 1.0 / R::rgamma(a_tau + 0.5 * E, 1.0 / (b_tau + 0.5 * accu(square(s.Boff))));
      if (free_model) {
        s.rho = R::rbeta(a_rho + E, b_rho + 2.0 * Ppairs - E);
      } else {
        const double n_both = 0.5 * (double) accu(s.delta % s.delta.t());
        const double n_single = E - 2.0 * n_both;
        s.rho = R::rbeta(a_rho + n_single, b_rho + Ppairs - n_single - n_both);
      }

      // adattamento delle scale del livello (solo in burn-in)
      if (it <= burnin && it % 50 == 0) {
        const double d = std::min(0.1, 1.0 / std::sqrt((double) it / 50.0));
        for (uword t = 0; t < K; ++t) {
          const uword i = ei[t], j = ej[t];
          const double r = T.n_rw(i, j) > 0.0 ? T.acc_rw(i, j) / T.n_rw(i, j) : 0.44;
          T.s_rw(i, j) = std::min(10.0, std::max(1e-4, T.s_rw(i, j) * std::exp(d * (r - 0.44))));
        }
        for (uword i = 0; i < p; ++i) {
          const double r = T.n_a(i) > 0.0 ? T.acc_a(i) / T.n_a(i) : 0.44;
          T.s_a(i) = std::min(10.0, std::max(1e-4, T.s_a(i) * std::exp(d * (r - 0.44))));
        }
        T.acc_rw.zeros(); T.n_rw.zeros(); T.acc_a.zeros(); T.n_a.zeros();
      }
      if (it == burnin) {
        T.acc_rw.zeros(); T.n_rw.zeros(); T.acc_a.zeros(); T.n_a.zeros();
        T.acc_flip = 0.0; T.n_flip = 0.0; T.acc_rev = 0.0; T.n_rev = 0.0;
      }
    }

    // scambi fra livelli adiacenti
    for (uword lev = 0; lev + 1 < KT; ++lev) {
      const double l0 = st[perm[lev]].llcur, l1 = st[perm[lev + 1]].llcur;
      const double la = (betas(lev) - betas(lev + 1)) * (l1 - l0);
      swap_n(lev) += 1.0;
      if (std::log(R::unif_rand()) < la) { std::swap(perm[lev], perm[lev + 1]); swap_acc(lev) += 1.0; }
    }
    if (it == burnin) { swap_acc.zeros(); swap_n.zeros(); }

    if (it > burnin && (it - burnin) % thin == 0) {
      const int k = (it - burnin) / thin - 1;
      const PTState& s = st[perm[0]];
      Bsave.slice(k) = s.Bcur; Dsave.slice(k) = conv_to<mat>::from(s.delta);
      Asave.row(k) = s.avec.t();
      tausave(k) = s.tau; rhosave(k) = s.rho; llsave(k) = s.llcur;
    }

    if (verbose && it % report == 0)
      Rcpp::Rcout << "    iter " << it << " / " << n_iter
                  << "   loglik (beta=1) " << st[perm[0]].llcur << std::endl;
  }

  const PTTune& T0 = tu[0];
  return Rcpp::List::create(
    Rcpp::_["B"]      = Bsave,
    Rcpp::_["delta"]  = Dsave,
    Rcpp::_["a"]      = Asave,
    Rcpp::_["tau"]    = tausave,
    Rcpp::_["rho"]    = rhosave,
    Rcpp::_["loglik"] = llsave,
    Rcpp::_["betas"]  = betas,
    Rcpp::_["swap"]   = swap_acc / clamp(swap_n, 1.0, datum::inf),
    Rcpp::_["acc"]    = Rcpp::List::create(
        Rcpp::_["flip"] = T0.n_flip > 0 ? T0.acc_flip / T0.n_flip : NA_REAL,
        Rcpp::_["rev"]  = T0.n_rev > 0 ? T0.acc_rev / T0.n_rev : NA_REAL,
        Rcpp::_["rw"]   = T0.acc_rw / clamp(T0.n_rw, 1.0, datum::inf),
        Rcpp::_["a"]    = T0.acc_a / clamp(T0.n_a, 1.0, datum::inf)));
}
