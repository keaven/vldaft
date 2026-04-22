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

test_that("Both back ends agree for all five error distributions (theta = 0)", {
  d <- sim_aft_data()
  for (dist in c("weibull", "logistic", "normal", "cauchy")) {
    f_c <- vldaft(Surv(time, status) ~ x1 + x2, data = d,
                  dist = dist, theta = 0, backend = "c")
    f_r <- vldaft(Surv(time, status) ~ x1 + x2, data = d,
                  dist = dist, theta = 0, backend = "rust")
    expect_equal(unname(coef(f_c)), unname(coef(f_r)),
                 tolerance = 1e-5,
                 info = paste("distribution =", dist))
    expect_equal(f_c$loglik, f_r$loglik,
                 tolerance = 1e-6,
                 info = paste("distribution =", dist))
  }
})
