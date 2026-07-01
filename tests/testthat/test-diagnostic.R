test_that("compare_iv_cv returns a well-formed ivc object", {
  d <- ivc_simulate(n = 2000, gamma = 1, kappa = 0.3, seed = 21)
  fit <- compare_iv_cv(d, "Y", "P", "A", "Z", se = "mom")
  expect_s3_class(fit, "ivc")
  expect_named(fit, c("tau_IV","tau_comp","delta_hat",
                      "se_tau_IV","se_tau_comp","se_delta_hat",
                      "Delta","se","ci","p_value","reject","method",
                      "level","n","boot_reps","call"))
  expect_equal(fit$Delta, fit$tau_IV - fit$tau_comp, tolerance = 1e-12)
  expect_length(fit$ci, 2)
})

test_that("Delta has mean near zero when exogeneity holds (kappa = 0)", {
  set.seed(99)
  deltas <- replicate(200, {
    d <- ivc_simulate(n = 1500, gamma = 1, kappa = 0)
    fit <- compare_iv_cv(d, "Y", "P", "A", "Z", se = "mom")
    fit$Delta
  })
  expect_lt(abs(mean(deltas)), 0.03)
})

test_that("Type I error is controlled (conservative) at kappa = 0", {
  pw <- ivc_power(n = 1500, gamma = 1, kappa = 0, n_rep = 300, se = "mom", seed = 5)
  expect_lt(pw$rate, 0.10)
})

test_that("power increases with violation strength", {
  p_lo <- ivc_power(n = 1500, gamma = 1, kappa = 0.1, n_rep = 200, se = "mom", seed = 6)
  p_hi <- ivc_power(n = 1500, gamma = 1, kappa = 0.3, n_rep = 200, se = "mom", seed = 6)
  expect_gt(p_hi$rate, p_lo$rate)
})

test_that("bootstrap and mom SEs agree up to Monte Carlo error", {
  d <- ivc_simulate(n = 2000, gamma = 1, kappa = 0.3, seed = 31)
  f_mom  <- compare_iv_cv(d, "Y", "P", "A", "Z", se = "mom")
  f_boot <- compare_iv_cv(d, "Y", "P", "A", "Z", se = "bootstrap",
                          n_boot = 2000, seed = 1)
  expect_equal(f_boot$se / f_mom$se, 1, tolerance = 0.2)
})

test_that("individual SEs are exposed, consistent, and stable", {
  d <- ivc_simulate(n = 2000, gamma = 1, kappa = 0.3, seed = 31)
  f <- compare_iv_cv(d, "Y", "P", "A", "Z", se = "mom")

  expect_true(all(is.finite(c(f$se_tau_IV, f$se_tau_comp, f$se_delta_hat))))

  expect_equal(f$Delta,        0.1986041184, tolerance = 1e-6)
  expect_equal(f$se,           0.0364280005, tolerance = 1e-6)
  expect_equal(f$se_tau_IV,    0.0219808395, tolerance = 1e-6)
  expect_equal(f$se_tau_comp,  0.0431913709, tolerance = 1e-6)
  expect_equal(f$se_delta_hat, 0.1117883589, tolerance = 1e-6)

  ei <- estimate_iv(d, "Y", "P", "A", "Z", se = TRUE)
  expect_equal(unname(ei["se"]), f$se_tau_IV, tolerance = 1e-10)

  Zc <- d$Z - mean(d$Z); Ac <- d$A - mean(d$A); Yc <- d$Y - mean(d$Y)
  tiv <- sum(Zc * Yc) / sum(Zc * Ac); u <- Yc - tiv * Ac
  hc0 <- sqrt(sum(Zc^2 * u^2) / (sum(Zc * Ac))^2)
  expect_equal(f$se_tau_IV, hc0, tolerance = 1e-3)
})
