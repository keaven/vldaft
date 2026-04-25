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

test_that("compiled cure evaluator handles left-censored observations", {
  d <- sim_cure_data(n = 40, dist = "weibull", seed = 106)
  detection_limit <- stats::quantile(d$time, 0.3)
  d$left_time <- pmax(d$time, detection_limit)
  d$left_status <- as.integer(d$time > detection_limit)
  y <- Surv(d$left_time, d$left_status, type = "left")
  fam <- Cure_Weibull()

  Xmu <- model.matrix(~ x1 + x2, d)
  Xsc <- model.matrix(~ 1, d)
  Xcu <- model.matrix(~ x1 + x2, d)
  beta <- c(log(mean(d$sigma_true)),
            unname(coef(lm(mu_true ~ x1 + x2, d))),
            qlogis(0.2), 0, 0)

  dm <- cbind(Xsc, Xmu, Xcu, time = d$left_time,
              status = ifelse(d$left_status == 1L, 1L, -1L),
              start = 0, t0 = 0, time2 = 0)
  c_val <- .Call(
    "vldaft_cure_eval",
    dm,
    ncol(dm) - 5L,
    ncol(dm) - 4L,
    -1L,
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

test_that("cure convergence follows optim convergence", {
  d <- sim_cure_data(n = 80, dist = "weibull", seed = 105)
  expect_warning(
    fit <- vldaft(
      Surv(time, status) ~ x1 + x2 | 1 | x1 + x2,
      data = d,
      dist = "weibull",
      cure = TRUE,
      backend = "c",
      control = list(maxiter = 1)
    ),
    "optimizer did not converge"
  )

  expect_false(fit$converged)
  expect_equal(fit$optim_convergence, 1)
})

test_that("standalone cure models accept interval-censored data", {
  d <- sim_cure_data(n = 80, dist = "weibull", seed = 107)
  d$lower <- pmax(d$time * 0.8, .Machine$double.eps)
  d$upper <- d$time * 1.2

  fit <- vldaft(
    Surv(lower, upper, type = "interval2") ~ x1 + x2 | 1 | x1 + x2,
    data = d,
    dist = "weibull",
    cure = TRUE,
    backend = "c",
    control = list(maxiter = 40)
  )

  expect_s3_class(fit, "vldaft")
  expect_equal(fit$backend, "c")
  expect_true(all(fit$status == 2L))
  expect_true(is.finite(fit$loglik))
  expect_equal(fit$hessian_method, "optimHess")
})

test_that("reported cure intercept is on the raw covariate scale", {
  d <- sim_cure_data(n = 200, dist = "weibull", seed = 104)
  fit <- vldaft(
    Surv(time, status) ~ x1 + x2 | 1 | x1 + x2,
    data = d,
    dist = "weibull",
    cure = TRUE,
    backend = "c"
  )

  X_cure <- model.matrix(~ x1 + x2, d)
  cure_cf <- coef(fit)[grepl("^cure:", names(coef(fit)))]
  cure_cf <- cure_cf[paste0("cure:", colnames(X_cure))]

  expect_equal(
    stats::plogis(as.numeric(X_cure %*% unname(cure_cf))),
    as.numeric(predict(fit, newdata = d, type = "cure")),
    tolerance = 1e-8
  )
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
