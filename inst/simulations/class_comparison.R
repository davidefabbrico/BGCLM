# ============================================================
# Confronto fra BGCLM (modello A1) e gclm (procedura di Varando & Hansen, 200
# divisioni) sul grafo di classe: i dati identificano la classe di modello,
# quindi entrambi i metodi vengono portati al grafo di classe (class_graph) e
# valutati su cio' che e' identificabile.
#
# Scenari: gershgorin, weak, random (stessi dati e stessi semi di compare_gclm.R
# e gclm_vh_sim.R, quindi le metriche gia' pubblicate si riproducono) e modules,
# nuovo, in cui il grafo vero contiene una cricca-modulo: un arco super-covered,
# un triangolo con un figlio comune, o un arco isolato con un cammino.
#
# Uso:  Rscript inst/simulations/class_comparison.R <bgclm|gclm> [scenari...]
# Scrive inst/simulations/results/data/class_<metodo>_<scenario>.rds
# ============================================================
Sys.setenv(OMP_NUM_THREADS = "1")
suppressMessages({ library(BGCLM); library(gclm) })
src <- "inst/simulations/compare_gclm.R"
lines <- readLines(src)
eval(parse(text = paste(lines[1:(grep("^GEN <- list", lines) - 1)], collapse = "\n")))
for (e in parse(src)) if (is.call(e) && identical(e[[1]], as.name("<-")) && is.call(e[[3]]) &&
                          identical(e[[3]][[1]], as.name("function"))) eval(e)   # metrics, shat, safe_lyap, ...
N_REP <- as.integer(Sys.getenv("CLASS_N_REP", N_REP))                      # per prove brevi
N_VALS <- as.integer(strsplit(Sys.getenv("CLASS_N_VALS", paste(N_VALS, collapse = ",")), ",")[[1]])
a <- commandArgs(TRUE); METHOD <- a[1]
RUN <- if (length(a) > 1) a[-1] else c("gershgorin", "weak", "random", "modules")
SEED_SCEN <- c(SCENARIOS, "modules")                 # modules in fondo: i semi degli altri non cambiano
OUT <- "inst/simulations/results/data"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# scenario modules: tre modelli di grafo con una cricca-modulo, a rotazione sulle repliche
gen_modules <- function(p, n_edge, lo = 0.4, hi = 0.8, a = 0.5) {
  E <- function(...) matrix(c(...), ncol = 2, byrow = TRUE)          # righe (da, a)
  tpl <- (CUR_REP - 1) %% 3 + 1
  if (tpl == 1) {                                     # arco super-covered 1 - 4
    e <- E(1, 4, 3, 1, 3, 4, 5, 1, 5, 4); if (runif(1) < 0.5) e[1, ] <- c(4, 1)
  } else if (tpl == 2) {                              # triangolo {1,2,3} orientato a caso, figlio comune 4
    tri <- E(1, 2, 2, 3, 1, 3); fl <- runif(3) < 0.5; tri[fl, ] <- tri[fl, 2:1, drop = FALSE]
    e <- rbind(tri, E(1, 4, 2, 4, 3, 4, 5, 4))
  } else {                                            # arco isolato 1 - 2 e cammino 3 - 4 - 5
    e <- E(1, 2, 3, 4, 4, 5); fl <- runif(3) < 0.5; e[fl, ] <- e[fl, 2:1, drop = FALSE]
  }
  perm <- sample(p); e[] <- perm[e]
  Boff <- matrix(0, p, p); Boff[e[, 2:1, drop = FALSE]] <- sample(c(-1, 1), nrow(e), TRUE) * runif(nrow(e), lo, hi)
  B <- build_B(Boff, rep(a, p))
  list(B = B, C = diag(p), delta = (Boff != 0) * 1L, Sigma = lyapunov(B, diag(p)), template = tpl)
}
GEN2 <- list(gershgorin = gen_gershgorin, weak = gen_weak, random = gen_random, modules = gen_modules)

VH_LAMBDAS <- 6 * 10^seq(-4, 0, length.out = 100); N_SPLIT <- 200; SOGLIA <- 0.85
vh_select <- function(Y) {                           # come gclm_vh_sim.R, con anche la frequenza per coppia
  p <- ncol(Y); freq <- pfreq <- matrix(0, p, p)
  for (r in seq_len(N_SPLIT)) {
    ix <- sample(nrow(Y), floor(nrow(Y) / 2))
    Str <- cov(Y[-ix, ]); dd <- diag(1 / sqrt(diag(Str))); Ste <- dd %*% cov(Y[ix, ]) %*% dd
    path <- tryCatch(gclm.path(cov2cor(Str), lambdas = VH_LAMBDAS, B = -0.5 * solve(cov2cor(Str)), C = rep(1, p),
                               lambdac = 0.01, eps = 1e-6, maxIter = 1000), error = function(e) NULL)
    if (is.null(path)) next
    path <- Filter(function(z) !is.null(z$B) && all(is.finite(z$B)), path)
    refit <- lapply(path, function(z) gclm(cov2cor(Str), B = z$B, C = z$C, lambda = 0, lambdac = -1, eps = 1e-10, job = 10))
    mll <- sapply(refit, function(z) { P <- tryCatch(solve(z$Sigma), error = function(e) NULL)
      if (is.null(P) || !all(is.finite(P))) Inf else -determinant(P, logarithm = TRUE)$modulus + sum(Ste * P) })
    s <- (abs(refit[[which.min(mll)]]$B) > 1e-8) * 1; diag(s) <- 0
    freq <- freq + s; pfreq <- pfreq + ((s + t(s)) > 0)
  }
  list(freq = freq / N_SPLIT, pfreq = pfreq / N_SPLIT)
}

up <- upper.tri(diag(P))
for (SC in RUN) {
  res <- list(); grid <- expand.grid(rep = seq_len(N_REP), n = N_VALS)
  for (g in seq_len(nrow(grid))) {
    n <- grid$n[g]; rep_i <- grid$rep[g]; CUR_REP <- rep_i
    set.seed(1000 * match(SC, SEED_SCEN) + 10 * rep_i)
    truth <- GEN2[[SC]](P, N_EDGE); Y <- rmvn0(n, truth$Sigma)
    cg_true <- class_graph(truth$delta); sk_true <- (truth$delta | t(truth$delta))[up]
    t0 <- proc.time()[3]
    if (METHOD == "bgclm") {
      set.seed(77 + rep_i)
      f <- gclm_mcmc(Y, n_iter = N_ITER_PT, burnin = BURNIN_PT, thin = THIN_PT, simple = TRUE,
                     eps_bridge = EPS_BRIDGE, pt_betas = PT_BETAS, verbose = FALSE)
      sel <- select_graph(f)
      keys <- apply(f$delta, 3, function(d) paste0("g", paste(which(d == 1), collapse = ",")))
      cls <- vapply(model_class(truth$delta), function(d) paste0("g", paste(which(d == 1), collapse = ",")), "")
      qs <- apply(f$B, c(1, 2), stats::quantile, probs = c(0.025, 0.975))
      m <- metrics(f$B_mean, f$Sigma_mean, f$PIP[which(diag(P) == 0)], truth, qs = qs, sel = sel)
      extra <- data.frame(auc_pair = auc((f$PIP + t(f$PIP))[up], sk_true * 1),
                          p_true_graph = mean(keys == paste0("g", paste(which(truth$delta == 1), collapse = ","))),
                          p_true_class = mean(keys %in% cls))
      store <- list(PIP = f$PIP, draws = table(keys))
    } else {
      set.seed(4000 + 1000 * match(SC, SEED_SCEN) + 10 * rep_i + n)       # come gclm_vh_sim.R
      v <- vh_select(Y); sel <- (v$freq >= SOGLIA) * 1L
      fit <- tryCatch(gclm(shat(Y), B = -0.5 * diag(P) + 0.1 * sel, C = rep(1, P), lambda = 0, lambdac = 0.01,
                           eps = 1e-10, maxIter = 5000, job = 10), error = function(e) NULL)
      B_hat <- if (is.null(fit)) matrix(NA, P, P) else fit$B
      S_hat <- if (is.null(fit)) NULL else safe_lyap(fit$B, diag(fit$C, nrow = P))
      m <- metrics(B_hat, S_hat, v$freq[which(diag(P) == 0)], truth, sel = sel)
      extra <- data.frame(auc_pair = auc(v$pfreq[up], sk_true * 1), p_true_graph = NA_real_, p_true_class = NA_real_)
      store <- v
    }
    cm <- class_metrics(class_graph(sel), cg_true)
    names(cm) <- paste0("cl_", names(cm))
    res[[g]] <- list(row = cbind(scen = SC, n = n, rep = rep_i, method = METHOD,
                                 template = if (is.null(truth$template)) NA else truth$template,
                                 nontrivial = any(cg_true & t(cg_true)), m, extra, cm, time = proc.time()[3] - t0),
                     sel = sel, truth = truth$delta, store = store)
    message(sprintf("%s %s n = %d rep %d: %.0f s", METHOD, SC, n, rep_i, proc.time()[3] - t0))
  }
  saveRDS(res, file.path(OUT, sprintf("class_%s_%s.rds", METHOD, SC)))
  tab <- do.call(rbind, lapply(res, `[[`, "row"))
  print(aggregate(tab[c("auc_pair", "cl_tpr_sk", "cl_fpr_sk", "cl_orient", "cl_overorient", "cl_shd", "cl_exact", "p_true_class")],
                  tab["n"], mean, na.rm = TRUE), digits = 3)
}
