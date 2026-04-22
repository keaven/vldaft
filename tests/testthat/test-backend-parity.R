library(survival)

test_that("C and Rust back ends agree to ~ double precision", {
  d <- sim_aft_data()

  f_c <- vldaft(Surv(time, status) ~ x1 + x2, data = d,
                dist = "weibull", theta = 1, backend = "c")
  f_r <- vldaft(Surv(time, status) ~ x1 + x2, data = d,
                dist = "weibull", theta = 1, backend = "rust")

  expect_equal(unname(coef(f_c)), unname(coef(f_r)), tolerance = 1e-6)
  expect_equal(f_c$loglik, f_r$loglik, tolerance = 1e-8)
  expect_equal(unname(diag(vcov(f_c))), unname(diag(vcov(f_r))),
               tolerance = 1e-6)
})

test_that("Both back ends agree for all four log-time distributions (theta = 0)", {
  d <- sim_aft_data()
  ## Tolerance budget reflects the upstream CDF implementation:
  ##   - Weibull / logistic: closed form, machine precision agreement.
  ##   - Normal / Cauchy   : R's pnorm/pcauchy use log.p-aware tail
  ##                         expansions; libm::erfc / atan in Rust are
  ##                         accurate to ~1e-15 in the bulk but lose
  ##                         precision in the deep tails the optimizer
  ##                         visits during step-halving. They still pin
  ##                         down the same local maximum to ~1e-3.
  tol <- list(
    weibull  = list(coef = 1e-6, loglik = 1e-7),
    logistic = list(coef = 1e-6, loglik = 1e-7),
    normal   = list(coef = 1e-3, loglik = 1e-2),
    cauchy   = list(coef = 1e-3, loglik = 1e-2)
  )
  for (dist in names(tol)) {
    f_c <- suppressWarnings(vldaft(
      Surv(time, status) ~ x1 + x2, data = d,
      dist = dist, theta = 0, backend = "c"))
    f_r <- suppressWarnings(vldaft(
      Surv(time, status) ~ x1 + x2, data = d,
      dist = dist, theta = 0, backend = "rust"))
    expect_equal(unname(coef(f_c)), unname(coef(f_r)),
                 tolerance = tol[[dist]]$coef,
                 info = paste("distribution =", dist))
    expect_equal(f_c$loglik, f_r$loglik,
                 tolerance = tol[[dist]]$loglik,
                 info = paste("distribution =", dist))
  }
})
