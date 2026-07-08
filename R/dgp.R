#' Simulate a pretest-posttest dataset with a shared confounder
#'
#' Generates data under the single-factor confounding model used in the IVC
#' simulation study. A latent confounder \eqn{U} loads on the pretest \eqn{P},
#' the posttest \eqn{Y}, and the treatment \eqn{A}; the instrument \eqn{Z}
#' shifts the treatment with strength `gamma`, and violates exogeneity through
#' a path into \eqn{U} of strength `kappa` (`kappa = 0` is a valid instrument).
#'
#' The structural model is
#' \deqn{U = \kappa Z + e_U,\quad A = \gamma Z + \beta_A U + e_A,}
#' \deqn{P = \beta_P U + e_P,\quad Y = \tau A + \beta_Y U + e_Y,}
#' with all error terms standard normal (scaled by the corresponding `sd_*`).
#'
#' @param n Sample size.
#' @param gamma Instrument strength (path \eqn{Z \to A}).
#' @param kappa Exogeneity-violation strength (path \eqn{Z \to U});
#'   `kappa = 0` yields a valid instrument.
#' @param tau True treatment effect. Default `0.5`.
#' @param beta_A,beta_P,beta_Y Loadings of \eqn{U} on \eqn{A}, \eqn{P}, \eqn{Y}.
#' @param sd_eU,sd_eP,sd_eY,sd_eA Error standard deviations.
#' @param seed Optional integer seed.
#'
#' @return A data frame with columns `Z`, `U`, `A`, `P`, `Y` (the latent `U`
#'   is returned for reference and is not used by the estimators).
#' @examples
#' d <- ivc_simulate(n = 500, gamma = 1, kappa = 0.2, seed = 42)
#' head(d)
#' @export
ivc_simulate <- function(n, gamma, kappa, tau = 0.5,
                         beta_A = 1, beta_P = 1, beta_Y = 1,
                         sd_eU = 1, sd_eP = 1, sd_eY = 1, sd_eA = 1,
                         seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  Z  <- stats::rnorm(n)
  eU <- stats::rnorm(n, sd = sd_eU)
  eP <- stats::rnorm(n, sd = sd_eP)
  eY <- stats::rnorm(n, sd = sd_eY)
  eA <- stats::rnorm(n, sd = sd_eA)
  U <- kappa * Z + eU
  A <- gamma * Z + beta_A * U + eA
  P <- beta_P * U + eP
  Y <- tau * A + beta_Y * U + eY
  data.frame(Z = Z, U = U, A = A, P = P, Y = Y)
}

#' Simulation-based power analysis for the IVC diagnostic
#'
#' Estimates the detection (rejection) rate of the IVC diagnostic by repeatedly
#' simulating data from [ivc_simulate()] under a given design and recording how
#' often the diagnostic interval for \eqn{\Delta} excludes zero. When
#' `kappa = 0` the returned rate is an empirical Type I error rate; when
#' `kappa > 0` it is empirical power.
#'
#' @param n Sample size per replicate.
#' @param gamma Instrument strength.
#' @param kappa Exogeneity-violation strength (`0` = valid instrument).
#' @param n_rep Number of Monte Carlo replicates. Default `500`.
#' @param se Uncertainty method passed to [compare_iv_cv()]: `"mom"` (default,
#'   fast) or `"bootstrap"`.
#' @param level Confidence level. Default `0.95`.
#' @param n_boot Bootstrap replicates when `se = "bootstrap"`.
#' @param tau,beta_A,beta_P,beta_Y,sd_eU,sd_eP,sd_eY,sd_eA Passed to
#'   [ivc_simulate()].
#' @param seed Optional integer seed.
#'
#' @return A list of class `"ivc_power"` with the detection `rate`, the number
#'   of valid replicates `n_valid`, and the design settings.
#' @examples
#' \donttest{
#' ivc_power(n = 1000, gamma = 1, kappa = 0.2, n_rep = 100, se = "mom")
#' }
#' @export
ivc_power <- function(n, gamma, kappa, n_rep = 500,
                      se = c("mom", "bootstrap"), level = 0.95, n_boot = 1000,
                      tau = 0.5, beta_A = 1, beta_P = 1, beta_Y = 1,
                      sd_eU = 1, sd_eP = 1, sd_eY = 1, sd_eA = 1, seed = NULL) {
  se <- match.arg(se)
  if (!is.null(seed)) set.seed(seed)
  rej <- logical(n_rep); ok <- logical(n_rep)
  for (r in seq_len(n_rep)) {
    d <- ivc_simulate(n, gamma, kappa, tau = tau, beta_A = beta_A,
                      beta_P = beta_P, beta_Y = beta_Y, sd_eU = sd_eU,
                      sd_eP = sd_eP, sd_eY = sd_eY, sd_eA = sd_eA)
    fit <- tryCatch(
      # weak_loading_warn = 0: a power run intentionally explores weak designs;
      # repeating the per-fit warning hundreds of times would only be noise.
      compare_iv_cv(d, "Y", "P", "A", "Z", se = se, level = level, n_boot = n_boot,
                    weak_loading_warn = 0),
      error = function(e) NULL)
    if (is.null(fit) || !is.finite(fit$ci[1]) || !is.finite(fit$ci[2])) {
      ok[r] <- FALSE; next
    }
    ok[r] <- TRUE; rej[r] <- fit$reject
  }
  out <- list(
    rate = mean(rej[ok]), n_valid = sum(ok), n_rep = n_rep,
    n = n, gamma = gamma, kappa = kappa, se = se, level = level
  )
  class(out) <- "ivc_power"
  out
}
