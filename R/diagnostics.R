# Conditional local independence (tetrad) test for compass-variable validity.
# Implements Lee & Kim (2026), JEEV 39(1), 55-75: with two candidate compass
# variables C1, C2, conditional local independence implies that within each
# treatment level g the tetrad
#   tau(g) = Cov(P,C1|g) Cov(Y,C2|g) - Cov(P,C2|g) Cov(Y,C1|g)
# vanishes. The sample tetrad is tested against zero with a multivariate
# delta-method SE; p-values are adjusted across treatment levels.

#' Conditional local independence (tetrad) test for candidate compass variables
#'
#' Tests the conditional local independence condition required for a valid
#' compass (detection) variable in adjusted difference-in-differences and in
#' the compass-corrected estimate of [estimate_cv()]. With two or more
#' candidate compass variables, every pair is screened: for each pair and each
#' treatment level, the tetrad statistic is computed and tested against zero
#' (Lee & Kim, 2026, JEEV 39(1)). A pair for which the null is *not* rejected
#' at any treatment level shows no empirical evidence against conditional
#' local independence.
#'
#' The test is global over the four variables involved: a rejection indicates
#' a structural inconsistency somewhere in `(P, Y, C1, C2)` but cannot single
#' out which candidate violates the assumption (identification limit noted in
#' the source paper).
#'
#' @param data A data frame.
#' @param pretest Name of the pretest column (\eqn{P}).
#' @param posttest Name of the posttest column (\eqn{Y}).
#' @param treat Name of the treatment column; treated as a discrete factor
#'   whose levels define the strata \eqn{g}.
#' @param compass Character vector (length >= 2) of candidate compass variable
#'   columns. All pairs are tested.
#' @param adjust Multiplicity adjustment across the tests performed *within a
#'   pair* (i.e., across treatment levels), passed to [stats::p.adjust()].
#'   Default `"bonferroni"`, matching the source paper. Adjustment across
#'   pairs is intentionally not applied: the screening use recommended by the
#'   paper selects the pair with the weakest evidence of violation, and
#'   cross-pair adjustment would not change that ordering.
#' @param min_n Minimum observations per treatment level required to attempt
#'   the test. Default `30`.
#'
#' @return A data frame of class `"ivc_cli"` with one row per (pair,
#'   treatment level): columns `c1`, `c2`, `level`, `n`, `tetrad`, `se`, `z`,
#'   `p_value`, `p_adj`. The pair-level decision is "no evidence of violation"
#'   when all `p_adj` for that pair exceed the chosen alpha.
#' @references Lee, S., & Kim, Y. (2026). Validation of adjusted
#'   difference-in-differences under common trend assumption violations:
#'   Application of the conditional local independence test.
#'   *Journal of Educational Evaluation, 39*(1), 55-75.
#' @seealso [estimate_cv()], [compare_iv_cv()]
#' @examples
#' d <- ivc_simulate(n = 1000, gamma = 1, kappa = 0.2, seed = 1)
#' # binary treatment strata and two noisy proxies of U as candidates
#' d$A_bin <- as.integer(d$A > stats::median(d$A))
#' d$C1 <- d$U + rnorm(nrow(d)); d$C2 <- d$U + rnorm(nrow(d))
#' ivc_cli_test(d, pretest = "P", posttest = "Y", treat = "A_bin",
#'              compass = c("C1", "C2"))
#' @export
ivc_cli_test <- function(data, pretest, posttest, treat, compass,
                         adjust = "bonferroni", min_n = 30) {
  req <- c(pretest, posttest, treat, compass)
  miss <- setdiff(req, names(data))
  if (length(miss) > 0L) {
    stop("variable(s) not found in data: ", paste(miss, collapse = ", "), call. = FALSE)
  }
  if (length(compass) < 2L) {
    stop("the tetrad test needs at least two candidate compass variables.", call. = FALSE)
  }
  g <- data[[treat]]
  levels_g <- sort(unique(g[!is.na(g)]))
  if (length(levels_g) < 2L) {
    warning("treatment has a single level; testing within that level only.", call. = FALSE)
  }
  pairs <- utils::combn(compass, 2, simplify = FALSE)

  rows <- list()
  for (pr in pairs) {
    c1 <- pr[1]; c2 <- pr[2]
    pvals <- numeric(0)
    tmp <- list()
    for (lv in levels_g) {
      sub <- data[!is.na(g) & g == lv,
                  c(pretest, posttest, c1, c2), drop = FALSE]
      sub <- sub[stats::complete.cases(sub), , drop = FALSE]
      n_g <- nrow(sub)
      res <- .ivc_tetrad_one(as.numeric(sub[[pretest]]),
                             as.numeric(sub[[posttest]]),
                             as.numeric(sub[[c1]]),
                             as.numeric(sub[[c2]]),
                             min_n = min_n)
      tmp[[length(tmp) + 1L]] <- data.frame(
        c1 = c1, c2 = c2, level = as.character(lv), n = n_g,
        tetrad = res["tetrad"], se = res["se"], z = res["z"],
        p_value = res["p"], row.names = NULL
      )
      pvals <- c(pvals, res["p"])
    }
    block <- do.call(rbind, tmp)
    block$p_adj <- stats::p.adjust(block$p_value, method = adjust)
    rows[[length(rows) + 1L]] <- block
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  attr(out, "adjust") <- adjust
  class(out) <- c("ivc_cli", class(out))
  out
}

# Tetrad statistic, delta-method SE, z, and p for one treatment stratum.
# Sigma (the asymptotic covariance of the four sample covariances) is built
# from the influence functions psi_i(a,b) = (a_i - abar)(b_i - bbar) - s_ab,
# which is distribution-free (no normality of the data required; cf. Bollen,
# 1990, on distribution-free vanishing-tetrad tests).
.ivc_tetrad_one <- function(P, Y, C1, C2, min_n = 30) {
  na4 <- c(tetrad = NA_real_, se = NA_real_, z = NA_real_, p = NA_real_)
  n <- length(P)
  if (n < min_n) return(na4)
  cP <- P - mean(P); cY <- Y - mean(Y)
  c1 <- C1 - mean(C1); c2 <- C2 - mean(C2)

  # sample covariances (divide by n for consistency with the IF construction)
  s_pc1 <- mean(cP * c1); s_yc2 <- mean(cY * c2)
  s_pc2 <- mean(cP * c2); s_yc1 <- mean(cY * c1)
  tet <- s_pc1 * s_yc2 - s_pc2 * s_yc1

  # influence functions of the four covariances, order:
  # (s_pc1, s_yc2, s_pc2, s_yc1)
  psi <- cbind(cP * c1 - s_pc1, cY * c2 - s_yc2,
               cP * c2 - s_pc2, cY * c1 - s_yc1)
  Sigma <- crossprod(psi) / n           # asymptotic cov of sqrt(n) * s-vector
  grad <- c(s_yc2, s_pc1, -s_yc1, -s_pc2)
  v <- as.numeric(t(grad) %*% Sigma %*% grad) / n
  if (!is.finite(v) || v <= 0) return(na4)
  se <- sqrt(v)
  z <- tet / se
  c(tetrad = tet, se = se, z = z, p = 2 * stats::pnorm(-abs(z)))
}

#' @export
print.ivc_cli <- function(x, digits = 4, alpha = 0.05, ...) {
  cat("Conditional local independence (tetrad) test\n")
  cat(sprintf("  p-value adjustment within pair: %s\n\n", attr(x, "adjust")))
  print.data.frame(x, digits = digits, row.names = FALSE)
  ok <- stats::aggregate(p_adj ~ c1 + c2, data = as.data.frame(x),
                         FUN = function(p) all(is.finite(p)) && all(p > alpha))
  cat("\nPair-level decision (all adjusted p >", alpha, "=> no evidence of violation):\n")
  for (i in seq_len(nrow(ok))) {
    cat(sprintf("  (%s, %s): %s\n", ok$c1[i], ok$c2[i],
                if (isTRUE(ok$p_adj[i])) "no evidence of violation"
                else "REJECTED (or not testable)"))
  }
  invisible(x)
}
