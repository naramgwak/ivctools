#' Conventional instrumental-variable (Wald) estimate
#'
#' Computes the just-identified instrumental-variable estimate of a treatment
#' effect in a pretest-posttest design. With no covariates this is the Wald
#' ratio \eqn{\widehat{\tau}_{IV} = \mathrm{cov}(Z, Y) / \mathrm{cov}(Z, A)};
#' when covariates are supplied they are partialled out by weighted least
#' squares (Frisch-Waugh-Lovell) before the ratio is formed.
#'
#' @param data A data frame containing the analysis variables.
#' @param outcome Name of the posttest outcome column (\eqn{Y}).
#' @param pretest Name of the pretest column (\eqn{P}). Accepted for a common
#'   interface with [estimate_cv()]; not used by the IV estimate itself.
#' @param treat Name of the treatment column (\eqn{A}).
#' @param instrument Name of the single instrument column (\eqn{Z}).
#' @param covariates Optional character vector of covariate column names to
#'   partial out. Default `NULL` reproduces the covariate-free Wald estimate.
#' @param weights Optional survey weights: either a numeric vector of length
#'   `nrow(data)` or the name of a column in `data`.
#'
#' @return A numeric scalar, the IV estimate \eqn{\widehat{\tau}_{IV}}.
#' @seealso [estimate_cv()], [compare_iv_cv()]
#' @examples
#' d <- ivc_simulate(n = 1000, gamma = 1, kappa = 0, seed = 1)
#' estimate_iv(d, outcome = "Y", pretest = "P", treat = "A", instrument = "Z")
#' @export
estimate_iv <- function(data, outcome, pretest, treat, instrument,
                        covariates = NULL, weights = NULL) {
  prep <- .ivc_prepare(data, outcome, pretest, treat, instrument, covariates, weights)
  est <- .ivc_estimates(prep)
  unname(est["tau_IV"])
}

#' Compass-corrected treatment-effect estimate
#'
#' Computes the compass (detection-variable corrected) estimate
#' \eqn{\widehat{\tau}_{comp}} for a pretest-posttest design. The instrument
#' \eqn{Z} is used as a detection variable: residualizing \eqn{Z}, \eqn{Y}, and
#' \eqn{P} on the treatment \eqn{A} yields a loading ratio
#' \eqn{\widehat{\delta}}, the posttest is corrected to
#' \eqn{Y^* = Y - \widehat{\delta} P}, and \eqn{\widehat{\tau}_{comp}} is the
#' coefficient of \eqn{A} in a regression of \eqn{Y^*} on \eqn{A}. Unlike
#' [estimate_iv()], this estimate targets the causal effect even when \eqn{Z}
#' violates exogeneity, provided it is a valid detection variable.
#'
#' @inheritParams estimate_iv
#' @param return_delta If `TRUE`, returns a named vector with both
#'   `tau_comp` and the loading ratio `delta_hat`. Default `FALSE`.
#'
#' @return A numeric scalar \eqn{\widehat{\tau}_{comp}}, or a named numeric
#'   vector when `return_delta = TRUE`.
#' @seealso [estimate_iv()], [compare_iv_cv()]
#' @examples
#' d <- ivc_simulate(n = 1000, gamma = 1, kappa = 0.2, seed = 1)
#' estimate_cv(d, outcome = "Y", pretest = "P", treat = "A",
#'             instrument = "Z", return_delta = TRUE)
#' @export
estimate_cv <- function(data, outcome, pretest, treat, instrument,
                        covariates = NULL, weights = NULL, return_delta = FALSE) {
  prep <- .ivc_prepare(data, outcome, pretest, treat, instrument, covariates, weights)
  est <- .ivc_estimates(prep)
  if (return_delta) return(est[c("tau_comp", "delta_hat")])
  unname(est["tau_comp"])
}

#' Compare IV and compass estimates: the IV-Compass (IVC) exogeneity diagnostic
#'
#' Implements the IVC diagnostic by contrasting the conventional IV estimate
#' \eqn{\widehat{\tau}_{IV}} with the compass-corrected estimate
#' \eqn{\widehat{\tau}_{comp}}. Under a valid instrument the two estimands
#' coincide, so the contrast \eqn{\Delta = \widehat{\tau}_{IV} -
#' \widehat{\tau}_{comp}} has expectation zero. A confidence interval for
#' \eqn{\Delta} that excludes zero is read as evidence against single-instrument
#' exogeneity. Uncertainty is quantified either by an adaptive percentile
#' bootstrap or by an analytic method-of-moments (GMM-style) standard error.
#'
#' @inheritParams estimate_iv
#' @param se Method for quantifying uncertainty in \eqn{\Delta}: `"bootstrap"`
#'   (adaptive percentile bootstrap, the default) or `"mom"` (analytic
#'   method-of-moments standard error with a normal interval).
#' @param level Confidence level for the interval. Default `0.95`.
#' @param n_boot Target number of valid bootstrap replicates (used when
#'   `se = "bootstrap"`). Default `2000`.
#' @param seed Optional integer seed for reproducible bootstrap resampling.
#'
#' @return An object of class `"ivc"`: a list with elements `tau_IV`,
#'   `tau_comp`, `delta_hat`, `Delta`, `se`, `ci` (length-2 vector),
#'   `p_value`, `reject` (logical, whether the interval excludes zero),
#'   `method`, `level`, `n`, and `call`. Has `print`, `summary`, and `plot`
#'   methods.
#' @seealso [estimate_iv()], [estimate_cv()], [ivc_power()]
#' @examples
#' d <- ivc_simulate(n = 1000, gamma = 1, kappa = 0.3, seed = 1)
#' fit <- compare_iv_cv(d, outcome = "Y", pretest = "P", treat = "A",
#'                      instrument = "Z", se = "mom")
#' fit
#' @export
compare_iv_cv <- function(data, outcome, pretest, treat, instrument,
                          covariates = NULL, weights = NULL,
                          se = c("bootstrap", "mom"),
                          level = 0.95, n_boot = 2000, seed = NULL) {
  se <- match.arg(se)
  cl <- match.call()
  prep <- .ivc_prepare(data, outcome, pretest, treat, instrument, covariates, weights)
  est <- .ivc_estimates(prep)
  if (any(!is.finite(est))) {
    stop("Estimates are not finite (instrument too weak or relevance collapses).",
         call. = FALSE)
  }
  Delta <- unname(est["Delta"])
  z_crit <- stats::qnorm(1 - (1 - level) / 2)

  if (se == "mom") {
    se_val <- .ivc_mom_se(prep, est)
    ci <- c(Delta - z_crit * se_val, Delta + z_crit * se_val)
    p_value <- 2 * stats::pnorm(-abs(Delta / se_val))
    reps <- NA_integer_
  } else {
    if (!is.null(seed)) set.seed(seed)
    bt <- .ivc_boot(prep, level = level, n_boot = n_boot)
    se_val <- bt$se
    ci <- bt$ci
    # bootstrap percentile two-sided p-value (proportion on the null side x2)
    if (!is.null(bt$draws) && length(bt$draws) > 0) {
      pp <- mean(bt$draws <= 0)
      p_value <- 2 * min(pp, 1 - pp)
    } else p_value <- NA_real_
    reps <- bt$reps
  }
  reject <- is.finite(ci[1]) && is.finite(ci[2]) && !(ci[1] <= 0 && 0 <= ci[2])

  out <- list(
    tau_IV = unname(est["tau_IV"]),
    tau_comp = unname(est["tau_comp"]),
    delta_hat = unname(est["delta_hat"]),
    Delta = Delta, se = se_val, ci = ci, p_value = p_value,
    reject = reject, method = se, level = level,
    n = prep$n, boot_reps = reps, call = cl
  )
  class(out) <- "ivc"
  out
}
