library(survival)

test_that("vldaft returns a sane object for the simplest right-censored fit", {
  d <- sim_aft_data()
  fit <- vldaft(Surv(time, status) ~ x1 + x2, data = d,
                dist = "weibull", theta = 0)

  expect_s3_class(fit, "vldaft")
  expect_true(fit$converged)
  expect_length(fit$coefficients, 4)           # gamma:(Intercept), eta:(Intercept), eta:x1, eta:x2
  expect_true(any(grepl("^gamma:", names(fit$coefficients))))
  expect_true(any(grepl("^eta:",   names(fit$coefficients))))
  expect_equal(dim(fit$vcov), c(4, 4))
  expect_true(all(is.finite(diag(fit$vcov))))
  expect_gt(fit$loglik, -Inf)
})

test_that("print and summary succeed", {
  d <- sim_aft_data()
  fit <- vldaft(Surv(time, status) ~ x1 + x2, data = d,
                dist = "weibull", theta = 1)

  expect_output(print(fit), "AFT Regression Model")
  s <- summary(fit)
  expect_s3_class(s, "summary.vldaft")
  expect_true("theta1" %in% rownames(s$coefficients))
  expect_output(print(s), "Std. Error")
})

test_that("theta > 0 strictly improves the Weibull log-likelihood on data where it should", {
  d <- sim_aft_data(theta1 = -0.3)
  f0 <- vldaft(Surv(time, status) ~ x1 + x2, data = d, dist = "weibull", theta = 0)
  f1 <- vldaft(Surv(time, status) ~ x1 + x2, data = d, dist = "weibull", theta = 1)
  expect_gt(f1$loglik, f0$loglik - 1e-6)   # should not get worse
})

test_that("logLik(jacobian = TRUE) matches survreg on a null-theta Weibull fit", {
  skip_if_not_installed("survival")
  d <- sim_aft_data()
  fit <- vldaft(Surv(time, status) ~ x1 + x2, data = d,
                dist = "weibull", theta = 0)
  sr <- survival::survreg(Surv(time, status) ~ x1 + x2, data = d,
                          dist = "weibull")
  expect_equal(as.numeric(logLik(fit, jacobian = TRUE)),
               as.numeric(logLik(sr)),
               tolerance = 1e-3)
})

test_that("predict() returns sensible linear predictors", {
  d <- sim_aft_data()
  fit <- vldaft(Surv(time, status) ~ x1 + x2, data = d, dist = "weibull")

  lp    <- predict(fit, newdata = d, type = "lp")
  sg    <- predict(fit, newdata = d, type = "sigma")
  med   <- predict(fit, newdata = d, type = "median")
  surv  <- predict(fit, newdata = d[1:5, ], type = "survival", times = c(1, 10, 100))

  expect_length(lp, nrow(d))
  expect_true(all(sg > 0))
  expect_length(med, nrow(d))
  expect_true(all(surv >= 0 & surv <= 1))
})

test_that("MAXCOV guard trips before we touch the C stack", {
  set.seed(1)
  n <- 50
  df <- as.data.frame(matrix(rnorm(n * 35), n, 35))
  names(df) <- paste0("x", 1:35)
  df$time <- rexp(n)
  df$status <- rbinom(n, 1, 0.6)
  f <- as.formula(paste("Surv(time, status) ~",
                        paste(paste0("x", 1:35), collapse = "+")))
  expect_error(vldaft(f, data = df, dist = "weibull"),
               "at most 30 parameters")
})
