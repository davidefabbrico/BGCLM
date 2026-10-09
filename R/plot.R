# ============================================================
# Disegno del grafo a posteriori.
#
# Convenzione (come in utils.R): PIP[j, i] e' la probabilita' dell'arco i -> j.
# ============================================================

# Bezier quadratica da P0 a P2, curvata verso la normale sinistra della corda.
# Il segno fisso della curvatura fa sì che i -> j e j -> i si separino invece di
# sovrapporsi.
bezier_arc <- function(p0, p2, curv, n = 120) {
  d <- p2 - p0
  len <- sqrt(sum(d^2))
  if (len < 1e-9) return(NULL)
  nrm <- c(-d[2], d[1]) / len
  p1 <- (p0 + p2) / 2 + curv * len * nrm
  t <- seq(0, 1, length.out = n)
  cbind((1 - t)^2 * p0[1] + 2 * (1 - t) * t * p1[1] + t^2 * p2[1],
        (1 - t)^2 * p0[2] + 2 * (1 - t) * t * p1[2] + t^2 * p2[2])
}

arrow_head <- function(a, b, size, col) {
  u <- b - a
  u <- u / sqrt(sum(u^2))
  v <- c(-u[2], u[1])
  polygon(c(b[1], b[1] - size * u[1] + 0.45 * size * v[1],
            b[1] - size * u[1] - 0.45 * size * v[1]),
          c(b[2], b[2] - size * u[2] + 0.45 * size * v[2],
            b[2] - size * u[2] - 0.45 * size * v[2]),
          col = col, border = col)
}

#' Disegna il grafo a posteriori con gli archi sfumati per PIP
#'
#' Ogni arco e' colorato su una scala di grigi: nero = PIP 1 (inclusione certa),
#' grigio chiaro = PIP vicino alla soglia. Lo spessore segue la stessa scala,
#' cosi' l'incertezza resta leggibile anche in stampa in bianco e nero.
#'
#' @param x oggetto di classe "gclm" (da gclm_mcmc) oppure una matrice di PIP
#' @param truth opzionale: output di un generatore o matrice delta 0/1 del vero.
#'   Se fornito, gli archi veri ricevono un alone colorato sotto l'arco stimato:
#'   alone senza arco scuro = falso negativo, arco scuro senza alone = falso
#'   positivo.
#' @param thr PIP sotto la quale l'arco non viene disegnato
#' @param sign se TRUE, tratteggia gli archi con B_incl negativo
#' @param labels etichette dei nodi (default 1..p)
#' @param main titolo
#' @param legend se TRUE aggiunge la barra di riferimento dei PIP
#' @param curv curvatura degli archi
#' @param accent colore dell'alone degli archi veri
plot_gclm_graph <- function(x, truth = NULL, thr = 0.10, sign = TRUE,
                            labels = NULL, main = NULL, legend = TRUE,
                            curv = 0.16, accent = "#E39B2E") {

  PIP <- if (inherits(x, "gclm")) x$PIP else as.matrix(x)
  Bin <- if (inherits(x, "gclm")) x$B_incl else NULL
  p   <- nrow(PIP)
  if (is.null(labels)) labels <- seq_len(p)

  delta <- NULL
  if (!is.null(truth)) delta <- if (is.list(truth)) truth$delta else as.matrix(truth)

  ang <- seq(pi / 2, by = -2 * pi / p, length.out = p)
  pos <- cbind(cos(ang), sin(ang))
  r   <- if (p <= 6) 0.145 else if (p <= 12) 0.115 else 0.085

  op <- par(mar = c(if (legend) 2.6 else 0.6, 0.6, if (is.null(main)) 0.6 else 2.2, 0.6))
  on.exit(par(op), add = TRUE)
  plot(NA, xlim = c(-1.3, 1.3), ylim = c(-1.28, 1.3), asp = 1,
       axes = FALSE, xlab = "", ylab = "", main = main)

  ed <- which(diag(p) == 0 & PIP > thr, arr.ind = TRUE)   # riga j, colonna i
  if (nrow(ed)) ed <- ed[order(PIP[ed]), , drop = FALSE]  # i PIP alti sopra

  for (k in seq_len(nrow(ed))) {
    j <- ed[k, 1]; i <- ed[k, 2]
    arc <- bezier_arc(pos[i, ], pos[j, ], curv)
    keep <- sqrt(rowSums((arc - rep(pos[i, ], each = nrow(arc)))^2)) > r &
            sqrt(rowSums((arc - rep(pos[j, ], each = nrow(arc)))^2)) > r
    arc <- arc[keep, , drop = FALSE]
    if (nrow(arc) < 3) next

    pip <- PIP[j, i]
    col <- gray(0.90 * (1 - pip))
    lwd <- 0.5 + 3.5 * pip
    lty <- if (sign && !is.null(Bin) && Bin[j, i] < 0) 2 else 1

    if (!is.null(delta) && delta[j, i] == 1)
      lines(arc, col = accent, lwd = lwd + 3.6, lty = 1)
    lines(arc, col = col, lwd = lwd, lty = lty)
    arrow_head(arc[nrow(arc) - 2, ], arc[nrow(arc), ], size = 0.075, col = col)
  }

  # archi veri mancati del tutto: restano solo come alone
  if (!is.null(delta)) {
    miss <- which(diag(p) == 0 & delta == 1 & PIP <= thr, arr.ind = TRUE)
    for (k in seq_len(nrow(miss))) {
      j <- miss[k, 1]; i <- miss[k, 2]
      arc <- bezier_arc(pos[i, ], pos[j, ], curv)
      keep <- sqrt(rowSums((arc - rep(pos[i, ], each = nrow(arc)))^2)) > r &
              sqrt(rowSums((arc - rep(pos[j, ], each = nrow(arc)))^2)) > r
      arc <- arc[keep, , drop = FALSE]
      if (nrow(arc) >= 3) lines(arc, col = accent, lwd = 4.3)
    }
  }

  symbols(pos[, 1], pos[, 2], circles = rep(r, p), inches = FALSE, add = TRUE,
          bg = "white", fg = gray(0.25))
  text(pos[, 1], pos[, 2], labels, cex = if (p <= 8) 1.05 else 0.8)

  if (legend) {
    xs <- seq(-0.55, 0.55, length.out = 100)
    yb <- -1.22; yt <- -1.15
    rect(xs[-100], yb, xs[-1], yt,
         col = gray(0.90 * (1 - seq(0, 1, length.out = 99))), border = NA)
    rect(-0.55, yb, 0.55, yt, border = gray(0.4))
    text(c(-0.55, 0.55), yb - 0.055, c("0", "1"), cex = 0.72)
    text(0, yt + 0.06, "MPPI", cex = 0.78)
  }
  invisible(pos)
}
