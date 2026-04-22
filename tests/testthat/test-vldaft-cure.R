library(survival)

test_that("compiled cure evaluator matches the R cure-family log-likelihood", {
  d <- sim_cure_data(n = 40, dist = "weibull", seed = 101)
  y <- Surv(d$time, d$status)
  fam <- Cure_Weibull()

  Xmu <- model.matrix(~ x1 + x2, d)
  Xsc <- model.matrix(~ 1, d)
  Xcu <- model.matrix(~ x1 + x2, d)

  beta <- c(log(mean(d$sigma_true)),
            unname(coef(lm(mu_true ~ x1 + x2, d))),
            qlogis(0.2), 0, 0)

  dm <- cbind(Xsc, Xmu, Xcu, time = d$time, status = d$status, start = 0, t0 = 0)
  c_val <- .Call(
    "vldaft_cure_eval",
    dm,
    ncol(dm) - 4L,
    ncol(dm) - 3L,
    -1L,
    -1L,
    as.integer(ncol(Xsc) + (0:(ncol(Xmu) - 1L))),
    as.integer(0:(ncol(Xsc) - 1L)),
    as.integer(ncol(Xsc) + ncol(Xmu) + (0:(ncol(Xcu) - 1L))),
    0L,
    -1L,
    1L,
    1.0,
    as.double(beta),
    0L
  )

  mu <- as.numeric(Xmu %*% beta[2:4])
  log_sigma <- rep(beta[1], nrow(d))
  cure_p <- stats::plogis(as.numeric(Xcu %*% beta[5:7]))
  r_ll <- sum(fam$pdf(y, list(mu = mu, sigma = exp(log_sigma), nu = cure_p), log = TRUE))

  expect_equal(c_val$loglik, r_ll, tolerance = 1e-8)
})

test_that("vldaft fits a standalone cure model and predicts cure probabilities", {
  d <- sim_cure_data(n = 250, dist = "weibull", seed = 102)
  fit <- vldaft(
    Surv(time, status) ~ x1 + x2 | 1 | x1 + x2,
    data = d,
    dist = "weibull",
    cure = TRUE,
    backend = "c"
  )

  expect_s3_class(fit, "vldaft")
  expect_true(fit$converged)
  expect_true(isTRUE(fit$cure))
  expect_true(any(grepl("^cure:", names(fit$coefficients))))
  expect_equal(length(coef(fit)), 7)

  p_hat <- predict(fit, newdata = d[1:10, ], type = "cure")
  s_hat <- predict(fit, newdata = d[1:5, ], type = "survival", times = c(1, 5, 10))

  expect_true(all(p_hat > 0 & p_hat < 1))
  expect_true(all(s_hat >= 0 & s_hat <= 1))
})

test_that("cure models are only implemented for the C backend", {
  d <- sim_cure_data(n = 30, dist = "weibull", seed = 103)
  expect_error(
    vldaft(
      Surv(time, status) ~ x1 + x2 | 1 | x1 + x2,
      data = d,
      dist = "weibull",
      cure = TRUE,
      backend = "rust"
    ),
    "currently implemented only for backend = 'c'"
  )
})
