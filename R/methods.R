#' @export
print.ivc <- function(x, ...) {
  cat("IV-Compass (IVC) exogeneity diagnostic\n")
  cat(sprintf("  tau_IV   = %.4f\n", x$tau_IV))
  cat(sprintf("  tau_comp = %.4f   (delta_hat = %.4f)\n", x$tau_comp, x$delta_hat))
  cat(sprintf("  Delta    = %.4f   (SE = %s)\n", x$Delta,
              ifelse(is.finite(x$se), sprintf("%.4f", x$se), "NA")))
  cat(sprintf("  %.0f%% CI = [%.4f, %.4f]   p = %s\n",
              100 * x$level, x$ci[1], x$ci[2],
              ifelse(is.finite(x$p_value), sprintf("%.4f", x$p_value), "NA")))
  cat(sprintf("  method = %s, n = %d\n", x$method, x$n))
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
#' Draws the two treatment-effect estimates on a common scale and displays the
#' diagnostic contrast \eqn{\Delta} with its confidence interval. A base-graphics
#' plot is used so the package has no plotting dependencies.
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
  ests <- c(x$tau_comp, x$tau_IV)
  labs <- c("tau_comp", "tau_IV")
  xr <- range(c(ests, x$ci, 0), na.rm = TRUE)
  pad <- diff(xr) * 0.15 + 1e-6
  graphics::plot(NA, xlim = c(xr[1] - pad, xr[2] + pad), ylim = c(0.5, 3.5),
                 yaxt = "n", xlab = "treatment-effect estimate", ylab = "",
                 main = "IV-Compass diagnostic", ...)
  graphics::axis(2, at = c(1, 2, 3),
                 labels = c("tau_comp", "tau_IV", "Delta = IV - comp"), las = 1)
  graphics::abline(v = 0, col = "grey70", lty = 2)
  graphics::points(ests, c(1, 2), pch = 19, cex = 1.3)
  graphics::text(ests, c(1, 2), labels = sprintf("%.3f", ests), pos = 3)
  # Delta with CI on row 3
  if (all(is.finite(x$ci))) {
    graphics::segments(x$ci[1], 3, x$ci[2], 3, lwd = 2)
    graphics::points(x$Delta, 3, pch = 18, cex = 1.6,
                     col = if (isTRUE(x$reject)) "red3" else "black")
    graphics::text(x$Delta, 3, labels = sprintf("%.3f", x$Delta), pos = 3)
  }
  invisible(x)
}

#' @rdname plot.ivc
#' @export
ivc_plot <- function(x, ...) {
  if (!inherits(x, "ivc")) stop("ivc_plot() expects an 'ivc' object.", call. = FALSE)
  plot(x, ...)
}
