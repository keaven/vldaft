library(survival)

test_that("parse_vldaft_formula splits on | correctly", {
  d <- sim_aft_data()
  parts <- vldaft:::parse_vldaft_formula(
    Surv(time, status) ~ x1 + x2 | x1, d)
  expect_s3_class(parts$response, "Surv")
  expect_equal(attr(terms(parts$loc_formula), "term.labels"), c("x1", "x2"))
  expect_equal(attr(terms(parts$scale_formula), "term.labels"), "x1")
})

test_that("a formula without | gives a NULL scale component", {
  d <- sim_aft_data()
  parts <- vldaft:::parse_vldaft_formula(
    Surv(time, status) ~ x1 + x2, d)
  expect_null(parts$scale_formula)
})

test_that("separate scale covariates are actually used by the fit", {
  d <- sim_aft_data()
  fit <- vldaft(Surv(time, status) ~ x1 + x2 | x1, data = d,
                dist = "weibull", theta = 0)
  expect_true("gamma:x1" %in% names(fit$coefficients))
})
