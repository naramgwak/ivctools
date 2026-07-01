test_that("estimate_iv reduces to the Wald ratio with no covariates", {
  d <- ivc_simulate(n = 3000, gamma = 1, kappa = 0.3, seed = 11)
  wald <- cov(d$Z, d$Y) / cov(d$Z, d$A)
  expect_equal(estimate_iv(d, "Y", "P", "A", "Z"), wald, tolerance = 1e-10)
})

test_that("estimate_cv reproduces the residualized compass estimator", {
  d <- ivc_simulate(n = 3000, gamma = 1, kappa = 0.3, seed = 12)
  Zr <- resid(lm(Z ~ A, data = d)); Yr <- resid(lm(Y ~ A, data = d)); Pr <- resid(lm(P ~ A, data = d))
  delta <- sum(Yr * Zr) / sum(Pr * Zr)
  tc <- unname(coef(lm(I(Y - delta * P) ~ A, data = d))["A"])
  out <- estimate_cv(d, "Y", "P", "A", "Z", return_delta = TRUE)
  expect_equal(unname(out["delta_hat"]), delta, tolerance = 1e-10)
  expect_equal(unname(out["tau_comp"]), tc, tolerance = 1e-10)
})

test_that("compass estimate recovers tau even under violation; IV is biased", {
  d <- ivc_simulate(n = 20000, gamma = 1, kappa = 0.4, tau = 0.5, seed = 13)
  expect_equal(estimate_cv(d, "Y", "P", "A", "Z"), 0.5, tolerance = 0.05)
  expect_gt(estimate_iv(d, "Y", "P", "A", "Z"), 0.6)  # biased upward
})

test_that("missing variables raise informative errors", {
  d <- ivc_simulate(n = 100, gamma = 1, kappa = 0, seed = 1)
  expect_error(estimate_iv(d, "nope", "P", "A", "Z"), "not found")
})
