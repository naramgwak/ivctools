# Internal helpers for the IV-Compass (IVC) diagnostic.
# These are not exported. The public interface is estimate_iv(),
# estimate_cv(), and compare_iv_cv().
#
# Faithfulness note:
#   With no covariates and unit weights, .ivc_residualize() returns the
#   variables unchanged, so .ivc_estimates() reduces *exactly* to the study 1
#   estimators:
#     tau_IV   = cov(Z, Y) / cov(Z, A)                         (Wald ratio)
#     delta    = sum(Yr * Zr) / sum(Pr * Zr), residualized on A (compass loading)
#     tau_comp = coef(lm((Y - delta * P) ~ A))["A"]
#   When covariates X (and/or weights) are supplied, Z, A, P, Y are first
#   partialled on X by weighted least squares (Frisch-Waugh-Lovell). This is a
#   standard extension; it is documented as such and is NOT part of the
#   simulation-validated core.

# Weighted residualization of a vector on a model matrix (FWL partialling).
# Returns y unchanged when xmat is NULL.
.ivc_resid_on_X <- function(y, xmat, w) {
  if (is.null(xmat)) return(y)
  fit <- stats::lm.wfit(x = xmat, y = y, w = w)
  stats::residuals(fit)
}

# Build a covariate model matrix (with intercept) from a data frame of controls,
# or return NULL when there are no covariates.
.ivc_build_xmat <- function(data, covariates) {
  if (is.null(covariates) || length(covariates) == 0L) return(NULL)
  miss <- setdiff(covariates, names(data))
  if (length(miss) > 0L) {
    stop("covariates not found in data: ", paste(miss, collapse = ", "), call. = FALSE)
  }
  xdf <- data[, covariates, drop = FALSE]
  stats::model.matrix(~ ., data = xdf)
}

# Extract and validate the core analysis variables into a clean numeric frame.
.ivc_prepare <- function(data, outcome, pretest, treat, instrument,
                         covariates = NULL, weights = NULL) {
  req <- c(outcome, pretest, treat, instrument)
  miss <- setdiff(req, names(data))
  if (length(miss) > 0L) {
    stop("variable(s) not found in data: ", paste(miss, collapse = ", "), call. = FALSE)
  }
  Y <- as.numeric(data[[outcome]])
  P <- as.numeric(data[[pretest]])
  A <- as.numeric(data[[treat]])
  Z <- as.numeric(data[[instrument]])

  if (is.null(weights)) {
    w <- rep(1, nrow(data))
  } else if (length(weights) == 1L && is.character(weights)) {
    if (!weights %in% names(data)) stop("weights column not found: ", weights, call. = FALSE)
    w <- as.numeric(data[[weights]])
  } else {
    w <- as.numeric(weights)
    if (length(w) != nrow(data)) stop("weights length must equal nrow(data).", call. = FALSE)
  }

  xmat <- .ivc_build_xmat(data, covariates)

  # complete-case across everything used
  parts <- list(Y, P, A, Z, w)
  cc <- Reduce(`&`, lapply(parts, is.finite))
  if (!is.null(xmat)) cc <- cc & stats::complete.cases(xmat)
  list(Y = Y[cc], P = P[cc], A = A[cc], Z = Z[cc], w = w[cc],
       xmat = if (is.null(xmat)) NULL else xmat[cc, , drop = FALSE],
       n = sum(cc))
}

# Core point estimates on an already-prepared (cc) list.
# Returns named numeric: tau_IV, tau_comp, Delta, delta_hat.
.ivc_estimates <- function(prep, tol = 1e-12) {
  Y <- prep$Y; P <- prep$P; A <- prep$A; Z <- prep$Z; w <- prep$w; xmat <- prep$xmat

  # FWL: partial out covariates X from all variables (no-op when xmat is NULL).
  Yx <- .ivc_resid_on_X(Y, xmat, w)
  Px <- .ivc_resid_on_X(P, xmat, w)
  Ax <- .ivc_resid_on_X(A, xmat, w)
  Zx <- .ivc_resid_on_X(Z, xmat, w)

  na_out <- c(tau_IV = NA_real_, tau_comp = NA_real_,
              Delta = NA_real_, delta_hat = NA_real_)

  # --- tau_IV: weighted Wald ratio (just-identified IV after FWL) ---
  sw <- sum(w)
  wcov <- function(a, b) sum(w * (a - sum(w * a) / sw) * (b - sum(w * b) / sw)) / sw
  den_iv <- wcov(Zx, Ax)
  if (!is.finite(den_iv) || abs(den_iv) < tol) return(na_out)
  tau_iv <- wcov(Zx, Yx) / den_iv

  # --- delta_hat: compass loading, residualize on treatment A then ratio ---
  # resid of Z, Y, P on A (FWL-adjusted, weighted)
  Aonly <- cbind(`(Intercept)` = 1, A = Ax)
  Zr <- stats::residuals(stats::lm.wfit(Aonly, Zx, w))
  Yr <- stats::residuals(stats::lm.wfit(Aonly, Yx, w))
  Pr <- stats::residuals(stats::lm.wfit(Aonly, Px, w))
  num <- sum(w * Yr * Zr)
  den <- sum(w * Pr * Zr)
  if (!is.finite(num) || !is.finite(den) || abs(den) < tol) return(na_out)
  delta_hat <- num / den

  # --- tau_comp: regress compass-adjusted outcome Y* = Y - delta*P on A ---
  Dstar <- Yx - delta_hat * Px
  cm <- stats::lm.wfit(Aonly, Dstar, w)
  tau_comp <- unname(stats::coef(cm)["A"])
  if (!is.finite(tau_comp)) return(na_out)

  c(tau_IV = as.numeric(tau_iv),
    tau_comp = as.numeric(tau_comp),
    Delta = as.numeric(tau_iv - tau_comp),
    delta_hat = as.numeric(delta_hat))
}

# Analytic method-of-moments (GMM-style) SE for Delta.
# Faithful to study_1_gamma_with_F.R::mom_se_delta (no-covariate case).
# When covariates are present, moments are formed on FWL-residualized variables.
.ivc_mom_se <- function(prep, est, tol = 1e-12) {
  Y <- prep$Y; P <- prep$P; A <- prep$A; Z <- prep$Z; w <- prep$w; xmat <- prep$xmat
  Yx <- .ivc_resid_on_X(Y, xmat, w); Px <- .ivc_resid_on_X(P, xmat, w)
  Ax <- .ivc_resid_on_X(A, xmat, w); Zx <- .ivc_resid_on_X(Z, xmat, w)
  n <- length(Yx)

  delta  <- est["delta_hat"]; tau_c <- est["tau_comp"]; tiv <- est["tau_IV"]
  Aonly <- cbind(1, Ax)
  Zr <- stats::residuals(stats::lm.wfit(Aonly, Zx, w))
  Yr <- stats::residuals(stats::lm.wfit(Aonly, Yx, w))
  Pr <- stats::residuals(stats::lm.wfit(Aonly, Px, w))

  g1 <- Zr * (Yr - delta * Pr)
  g2 <- Ax * (Yx - delta * Px - tau_c * Ax)
  g3 <- Zx * (Yx - tiv * Ax)

  G <- matrix(0, 3, 3)
  G[1, 1] <- -mean(Zr * Pr)
  G[2, 1] <- -mean(Ax * Px)
  G[2, 2] <- -mean(Ax^2)
  G[3, 3] <- -mean(Zx * Ax)

  gmat <- cbind(g1, g2, g3)
  S <- crossprod(gmat) / n
  Ginv <- tryCatch(solve(G), error = function(e) NULL)
  if (is.null(Ginv) || any(!is.finite(Ginv))) return(NA_real_)
  V <- (1 / n) * Ginv %*% S %*% t(Ginv)
  a <- c(0, -1, 1)              # Delta = tau_IV - tau_comp
  var_delta <- as.numeric(t(a) %*% V %*% a)
  if (!is.finite(var_delta) || var_delta < 0) return(NA_real_)
  sqrt(var_delta)
}

# Adaptive percentile bootstrap for Delta.
# Faithful to study_1::boot_ci_delta (resamples rows, recomputes Delta).
.ivc_boot <- function(prep, data_cc, est_args, level = 0.95,
                      n_boot = 2000, min_ok = 1000, chunk = 500, max_draw = NULL) {
  n <- prep$n
  if (is.null(max_draw)) max_draw <- max(n_boot, max_draw <- n_boot * 3)
  alpha <- c((1 - level) / 2, 1 - (1 - level) / 2)
  deltas <- numeric(0); drawn <- 0L
  target <- max(n_boot, min_ok)
  while (length(deltas) < target && drawn < max_draw) {
    b <- min(chunk, max_draw - drawn)
    vals <- vapply(seq_len(b), function(i) {
      idx <- sample.int(n, n, replace = TRUE)
      p2 <- list(Y = prep$Y[idx], P = prep$P[idx], A = prep$A[idx],
                 Z = prep$Z[idx], w = prep$w[idx],
                 xmat = if (is.null(prep$xmat)) NULL else prep$xmat[idx, , drop = FALSE],
                 n = n)
      e <- .ivc_estimates(p2)
      if (any(!is.finite(e))) NA_real_ else unname(e["Delta"])
    }, numeric(1))
    deltas <- c(deltas, vals[is.finite(vals)])
    drawn <- drawn + b
  }
  if (length(deltas) < min_ok) {
    return(list(ci = c(NA_real_, NA_real_), se = NA_real_, reps = length(deltas)))
  }
  ci <- stats::quantile(deltas, probs = alpha, names = FALSE, type = 7)
  list(ci = ci, se = stats::sd(deltas), reps = length(deltas), draws = deltas)
}
