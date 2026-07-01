test_that("compare_iv_cv returns a well-formed ivc object", {
  d <- ivc_simulate(n = 2000, gamma = 1, kappa = 0.3, seed = 21)
  fit <- compare_iv_cv(d, "Y", "P", "A", "Z", se = "mom")
  expect_s3_class(fit, "ivc")
  expect_named(fit, c("tau_IV","tau_comp","delta_hat","Delta","se","ci",
                      "p_value","reject","method","level","n","boot_reps","call"))
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
  expect_lt(pw$rate, 0.10)  # at/below nominal 5%, allowing MC noise
})

test_that("power increases with violation strength", {
  p_lo <- ivc_power(n = 1500, gamma = 1, kappa = 0.1, n_rep = 200, se = "mom", seed = 6)
  p_hi <- ivc_power(n = 1500, gamma = 1, kappa = 0.3, n_rep = 200, se = "mom", seed = 6)
  expect_gt(p_hi$rate, p_lo$rate)
})

test_that("estimators are numerically stable (seeded regression test)", {
  d <- ivc_simulate(n = 2000, gamma = 1, kappa = 0.3, seed = 31)

  # 해석적 값(MoM): RNG와 무관 -> 엄격 비교
  f_mom <- compare_iv_cv(d, "Y", "P", "A", "Z", se = "mom")
  expect_equal(f_mom$Delta, 0.1986041184, tolerance = 1e-6)
  expect_equal(f_mom$se,    0.0364280005, tolerance = 1e-6)

  # 부트스트랩 값: 시드 고정 -> 재현되되 RNG 구현 차이 대비 약한 허용오차
  f_boot <- compare_iv_cv(d, "Y", "P", "A", "Z", se = "bootstrap",
                          n_boot = 2000, seed = 1)
  expect_equal(f_boot$se, 0.0379913265, tolerance = 1e-3)
})
