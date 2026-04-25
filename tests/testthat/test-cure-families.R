library(survival)

shift_link_parameter <- function(par, target, delta) {
  out <- par
  if (target == "sigma") {
    out$sigma <- par$sigma * exp(delta)
  } else if (target == "nu") {
    out$nu <- stats::plogis(stats::qlogis(par$nu) + delta)
  } else {
    out[[target]] <- par[[target]] + delta
  }
  out
}

test_that("cure-family constructors expose the expected parameters", {
  cure_weibull <- Cure_Weibull(theta = 1)
  expect_equal(cure_weibull$names, c("mu", "sigma", "nu", "theta1"))
  expect_equal(unname(cure_weibull$links),
               c("identity", "log", "logit", "identity"))

  cure_exp <- Cure_Exponential()
  expect_equal(cure_exp$names, c("mu", "nu"))
  expect_equal(unname(cure_exp$links), c("identity", "logit"))

  par <- list(mu = c(0, 0), sigma = c(1, 1), nu = c(0.2, 0.4), theta1 = c(0, 0))
  expect_equal(cure_weibull$cure_probability(par), c(0.2, 0.4))
  surv <- cure_weibull$survival(c(1, 2, 4), par)
  expect_equal(dim(surv), c(2, 3))
  expect_true(all(diff(surv[1, ]) <= 0))
})

test_that("Cure_Weibull analytic scores match finite differences when theta = 0", {
  d <- sim_cure_data(n = 12, dist = "weibull", theta1 = 0, seed = 11)
  y <- Surv(d$time, d$status)
  fam <- Cure_Weibull(theta = 0)
  par <- list(
    mu = d$mu_true + 0.1,
    sigma = d$sigma_true * 1.05,
    nu = pmin(pmax(d$cure_true * 0.95, 0.05), 0.95)
  )

  total_loglik <- function(pars) {
    sum(fam$pdf(y, pars, log = TRUE))
  }

  for (target in c("mu", "sigma", "nu")) {
    eps <- 1e-6
    ll_plus <- total_loglik(shift_link_parameter(par, target, eps))
    ll_minus <- total_loglik(shift_link_parameter(par, target, -eps))
    num_grad <- (ll_plus - ll_minus) / (2 * eps)
    ana_grad <- sum(fam$score[[target]](y, par))
    expect_equal(ana_grad, num_grad, tolerance = 1e-4)
  }
})

test_that("AFT and cure families handle left and interval censoring", {
  y_left <- Surv(c(1, 2, 3, 4), c(1, 0, 1, 0), type = "left")
  y_interval <- Surv(c(1, 2, 3, 4, 4), c(2, Inf, 3, 5, NA),
                     type = "interval2")

  aft <- AFT_Weibull(theta = 0)
  cure <- Cure_Weibull(theta = 0)
  par_aft <- list(mu = rep(0.1, 5), sigma = rep(1.1, 5))
  par_cure <- list(mu = rep(0.1, 5), sigma = rep(1.1, 5), nu = rep(0.25, 5))

  expect_true(all(is.finite(aft$pdf(y_left, list(mu = rep(0.1, 4),
                                                  sigma = rep(1.1, 4)),
                                      log = TRUE))))
  expect_true(all(is.finite(aft$pdf(y_interval, par_aft, log = TRUE))))
  expect_true(all(is.finite(cure$pdf(y_interval, par_cure, log = TRUE))))
  expect_true(all(is.finite(aft$score$mu(y_interval, par_aft))))
  expect_true(all(is.finite(cure$score$nu(y_interval, par_cure))))
})

test_that("interval-censored cure-family scores match finite differences", {
  y <- Surv(c(1, 2, 3, 4, 4), c(2, Inf, 3, 5, NA), type = "interval2")
  fam <- Cure_Weibull(theta = 0)
  par <- list(
    mu = c(0.1, 0.2, -0.1, 0.3, 0.4),
    sigma = rep(1.1, 5),
    nu = rep(0.25, 5)
  )
  total_loglik <- function(pars) {
    sum(fam$pdf(y, pars, log = TRUE))
  }

  for (target in c("mu", "sigma", "nu")) {
    eps <- 1e-6
    ll_plus <- total_loglik(shift_link_parameter(par, target, eps))
    ll_minus <- total_loglik(shift_link_parameter(par, target, -eps))
    num_grad <- (ll_plus - ll_minus) / (2 * eps)
    ana_grad <- sum(fam$score[[target]](y, par))
    expect_equal(ana_grad, num_grad, tolerance = 1e-4)
  }
})

test_that("Cure_Exponential likelihood recovers intercept-only parameters", {
  set.seed(21)
  n <- 400
  mu_true <- 1.3
  p_true <- 0.25
  tau <- -log(p_true)
  u <- stats::runif(n)
  has_event <- u <= (1 - p_true)
  t_event <- rep(Inf, n)
  f0 <- -log1p(-u[has_event]) / tau
  t_event[has_event] <- exp(mu_true + log(-log(1 - f0)))
  c_time <- stats::rexp(n, rate = 0.02)
  time <- pmin(t_event, c_time)
  status <- as.integer(t_event <= c_time)

  fam <- Cure_Exponential()
  y <- Surv(time, status)
  objective <- function(eta) {
    par <- list(
      mu = rep(eta[1], n),
      nu = rep(stats::plogis(eta[2]), n)
    )
    -sum(fam$pdf(y, par, log = TRUE))
  }

  opt <- stats::optim(c(log(median(time)), stats::qlogis(0.15)),
                      objective, method = "BFGS")
  expect_equal(opt$convergence, 0)
  expect_equal(opt$par[1], mu_true, tolerance = 0.2)
  expect_equal(stats::plogis(opt$par[2]), p_true, tolerance = 0.08)
})

test_that("Cure_Exponential handles left truncation", {
  y <- Surv(c(1, 1.5, 2), c(4, 3, 3), c(1, 0, 1), type = "counting")
  fam <- Cure_Exponential()
  par <- list(mu = c(0, 0, 0), nu = c(0.3, 0.3, 0.3))

  ll <- fam$pdf(y, par, log = TRUE)
  sc_mu <- fam$score$mu(y, par)
  sc_nu <- fam$score$nu(y, par)

  expect_true(all(is.finite(ll)))
  expect_true(all(is.finite(sc_mu)))
  expect_true(all(is.finite(sc_nu)))
})
