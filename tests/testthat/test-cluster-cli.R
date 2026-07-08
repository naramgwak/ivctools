# Tests for the 2026-07 extensions: cluster-robust SEs (Lee & Kim, 2026,
# JEEV 39(2)) and the conditional local independence tetrad test (Lee & Kim,
# 2026, JEEV 39(1)).

.sim_school <- function(n_school = 50, m = 30, sigma_s = 2, tau = 0, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  sid <- rep(seq_len(n_school), each = m)
  Gs <- rbinom(n_school, 1, 0.5); A <- Gs[sid]
  us <- rnorm(n_school); S <- sigma_s * rnorm(n_school)
  U <- us[sid] + rnorm(n_school * m)
  Z <- U + rnorm(n_school * m)
  P <- 0.8 * U + S[sid] + rnorm(n_school * m)
  Y <- tau * A + 1.2 * U + S[sid] + rnorm(n_school * m)
  data.frame(A, Z, P, Y, sid)
}

test_that("cluster argument never changes the point estimate", {
  d <- .sim_school(seed = 101)
  a <- estimate_cv(d, "Y", "P", "A", "Z", weak_loading_warn = 0)
  b <- estimate_cv(d, "Y", "P", "A", "Z", cluster = "sid", weak_loading_warn = 0)
  expect_equal(a, b, tolerance = 1e-12)
})

test_that("cluster-robust SE exceeds the iid SE under school-level treatment", {
  d <- .sim_school(seed = 102)
  si <- estimate_cv(d, "Y", "P", "A", "Z", se = TRUE, weak_loading_warn = 0)
  sc <- estimate_cv(d, "Y", "P", "A", "Z", se = TRUE, cluster = "sid",
                    weak_loading_warn = 0)
  expect_gt(unname(sc["se_tau_comp"]), unname(si["se_tau_comp"]))
})

test_that("compare_iv_cv accepts cluster for both mom and bootstrap", {
  d <- .sim_school(n_school = 30, m = 20, seed = 103)
  fm <- compare_iv_cv(d, "Y", "P", "A", "Z", se = "mom", cluster = "sid")
  expect_true(is.finite(fm$se))
  expect_equal(fm$n_clusters, 30L)
  fb <- compare_iv_cv(d, "Y", "P", "A", "Z", se = "bootstrap",
                      cluster = "sid", n_boot = 1000, seed = 1)
  expect_true(is.finite(fb$se))
})

test_that("tetrad test controls false rejection under a valid compass pair", {
  set.seed(104)
  n <- 4000; U <- rnorm(n); G <- rbinom(n, 1, plogis(U))
  d <- data.frame(G = G,
                  P = 0.5 * U + rnorm(n), Y = 0.5 * U + rnorm(n),
                  C1 = 0.5 * U + rnorm(n), C2 = 0.5 * U + rnorm(n))
  tt <- ivc_cli_test(d, "P", "Y", "G", c("C1", "C2"))
  expect_true(all(is.finite(tt$p_adj)))
  # not a hard guarantee, but at n = 4000 with no violation this should hold
  expect_true(all(tt$p_adj > 0.01))
})

test_that("tetrad test detects a strong violation (residual corr = 0.5)", {
  set.seed(105)
  n <- 4000; U <- rnorm(n); G <- rbinom(n, 1, plogis(U))
  eP <- rnorm(n); eC1 <- 0.5 * eP + sqrt(1 - 0.25) * rnorm(n)
  d <- data.frame(G = G,
                  P = 0.5 * U + eP, Y = 0.5 * U + rnorm(n),
                  C1 = 0.5 * U + eC1, C2 = 0.5 * U + rnorm(n))
  tt <- ivc_cli_test(d, "P", "Y", "G", c("C1", "C2"))
  expect_true(any(tt$p_adj < 0.05))
})

test_that("weak compass loading triggers a warning", {
  set.seed(107)
  n <- 500
  d <- data.frame(A = rbinom(n, 1, 0.5), Z = rnorm(n),
                  P = rnorm(n), Y = rnorm(n))  # Z unrelated to P given A
  expect_warning(estimate_cv(d, "Y", "P", "A", "Z"), "weak compass loading")
})
