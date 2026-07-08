#' @export
print.ivc <- function(x, ...) {
  cat("IV-Compass (IVC) exogeneity diagnostic\n")
  cat(sprintf("  tau_IV   = %.4f   (SE = %s)\n", x$tau_IV,
              ifelse(is.finite(x$se_tau_IV), sprintf("%.4f", x$se_tau_IV), "NA")))
  cat(sprintf("  tau_comp = %.4f   (SE = %s)\n", x$tau_comp,
              ifelse(is.finite(x$se_tau_comp), sprintf("%.4f", x$se_tau_comp), "NA")))
  cat(sprintf("  delta_hat= %.4f   (SE = %s)\n", x$delta_hat,
              ifelse(is.finite(x$se_delta_hat), sprintf("%.4f", x$se_delta_hat), "NA")))
  cat(sprintf("  Delta    = %.4f   (SE = %s)\n", x$Delta,
              ifelse(is.finite(x$se), sprintf("%.4f", x$se), "NA")))
  cat(sprintf("  %.0f%% CI = [%.4f, %.4f]   p = %s\n",
              100 * x$level, x$ci[1], x$ci[2],
              ifelse(is.finite(x$p_value), sprintf("%.4f", x$p_value), "NA")))
  cat(sprintf("  method = %s%s, n = %d\n", x$method,
              if (!is.null(x$n_clusters) && is.finite(x$n_clusters))
                sprintf(" (cluster-robust, G = %d clusters)", x$n_clusters)
              else "",
              x$n))
  cat(sprintf("  Conclusion: %s\n",
              if (isTRUE(x$reject))
                "interval excludes 0 -> evidence against exogeneity"
              else
                "interval contains 0 -> no evidence against exogeneity"))
  invisible(x)
}

#' @export
summary.ivc <- function(object, ...) {
  print(object, ...)
  cat("\nNote: non-rejection does not prove instrument validity. The diagnostic\n")
  cat("is most informative when the instrument is relevant and the pretest-\n")
  cat("posttest structure yields a stable compass correction.\n")
  invisible(object)
}

#' @export
print.ivc_power <- function(x, ...) {
  lab <- if (x$kappa == 0) "Type I error rate" else "Power"
  cat(sprintf("IVC power analysis (%s)\n", x$se))
  cat(sprintf("  design: n = %d, gamma = %.3f, kappa = %.3f\n", x$n, x$gamma, x$kappa))
  cat(sprintf("  %s = %.3f  (valid reps: %d/%d)\n",
              lab, x$rate, x$n_valid, x$n_rep))
  invisible(x)
}

#' Plot the IVC diagnostic
#'
#' Draws a forest plot with all three quantities on a common scale, each with
#' its own confidence interval: the conventional IV estimate, the
#' compass-corrected estimate, and their contrast \eqn{\Delta}. When the
#' interval for \eqn{\Delta} excludes zero, that row is shaded and colored to
#' flag evidence against exogeneity. A base-graphics plot is used so the
#' package has no plotting dependencies.
#'
#' @param x An object of class `"ivc"` from [compare_iv_cv()].
#' @param ... Further arguments passed to the underlying plotting calls.
#' @return `x`, invisibly. Called for the side effect of drawing a plot.
#' @examples
#' d <- ivc_simulate(n = 1000, gamma = 1, kappa = 0.3, seed = 1)
#' fit <- compare_iv_cv(d, "Y", "P", "A", "Z", se = "mom")
#' plot(fit)
#' @export
plot.ivc <- function(x, ...) {
  z <- stats::qnorm(1 - (1 - x$level) / 2)

  est  <- c(x$tau_IV, x$tau_comp, x$Delta)
  lo   <- c(x$tau_IV - z * x$se_tau_IV, x$tau_comp - z * x$se_tau_comp, x$ci[1])
  hi   <- c(x$tau_IV + z * x$se_tau_IV, x$tau_comp + z * x$se_tau_comp, x$ci[2])
  labs <- c("tau_IV", "tau_comp", "Delta = IV - comp")
  y    <- c(3, 2, 1)

  if (!all(is.finite(c(lo, hi)))) {
    warning("Some confidence limits are not finite; plotting available values only.",
            call. = FALSE)
  }
  xr  <- range(c(lo, hi, est, 0), na.rm = TRUE)
  pad <- diff(xr) * 0.25 + 1e-6
  if (!is.finite(pad) || pad <= 0) pad <- 1
  xlim <- c(xr[1] - pad, xr[2] + pad)

  op <- graphics::par(mar = c(5.5, 9, 3.5, 2))
  on.exit(graphics::par(op))

  graphics::plot(NA, xlim = xlim, ylim = c(0.3, 3.9), yaxt = "n", xaxt = "n",
                 xlab = "", ylab = "", main = "IV-Compass diagnostic", ...)
  graphics::axis(1)
  graphics::mtext("treatment-effect estimate", side = 1, line = 2.5)
  graphics::axis(2, at = y, labels = labs, las = 1, tick = FALSE, line = -0.5)

  if (isTRUE(x$reject)) {
    graphics::rect(xlim[1], 0.55, xlim[2], 1.45,
                   col = grDevices::adjustcolor("red", alpha.f = 0.08), border = NA)
  }
  graphics::abline(h = 1.5, col = "grey85", lty = 1)
  graphics::abline(v = 0, col = "grey60", lty = 2)

  col_pt <- c("steelblue4", "darkorange3", if (isTRUE(x$reject)) "red3" else "grey30")

  ok <- is.finite(lo) & is.finite(hi)
  graphics::segments(lo[ok], y[ok], hi[ok], y[ok], lwd = 2.4, col = col_pt[ok])
  graphics::segments(lo[ok], y[ok] - 0.08, lo[ok], y[ok] + 0.08, col = col_pt[ok])
  graphics::segments(hi[ok], y[ok] - 0.08, hi[ok], y[ok] + 0.08, col = col_pt[ok])
  graphics::points(est, y, pch = c(16, 16, 18), cex = c(1.5, 1.5, 1.9), col = col_pt)

  lbl <- ifelse(ok, sprintf("%.3f  [%.3f, %.3f]", est, lo, hi), sprintf("%.3f  [NA]", est))
  graphics::text(ifelse(ok, hi, est), y, labels = lbl, pos = 4, cex = 0.82,
                 xpd = NA, col = col_pt)

  concl <- if (isTRUE(x$reject))
    sprintf("%.0f%% CI for Delta excludes 0 -> evidence against exogeneity", 100 * x$level)
  else
    sprintf("%.0f%% CI for Delta contains 0 -> no evidence against exogeneity", 100 * x$level)
  graphics::mtext(concl, side = 1, line = 4.2, cex = 0.85,
                  col = if (isTRUE(x$reject)) "red3" else "grey30")

  invisible(x)
}

#' @rdname plot.ivc
#' @export
ivc_plot <- function(x, ...) {
  if (!inherits(x, "ivc")) stop("ivc_plot() expects an 'ivc' object.", call. = FALSE)
  plot(x, ...)
}
