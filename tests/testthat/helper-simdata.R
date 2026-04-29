## Generate a small Weibull AFT dataset for tests.
##
## The ground truth is a varying-dispersion AFT:
##   log(T) = beta0 + beta1 * x1 + beta2 * x2 + sigma(x1) * eps,
##   log(sigma) = gamma0 + theta1 * (mu - mean(mu)),  eps ~ Gumbel(0,1).
##
## Uses a deterministic seed so tests are reproducible.
sim_aft_data <- function(n = 300, theta1 = -0.2, seed = 20240421) {
  set.seed(seed)
  x1 <- rnorm(n)
  x2 <- rbinom(n, 1, 0.4)
  beta <- c(b0 = 2.0, b1 = 0.5, b2 = -0.3)
  mu <- beta["b0"] + beta["b1"] * x1 + beta["b2"] * x2
  log_sigma <- -0.3 + theta1 * (mu - mean(mu))
  sigma <- exp(log_sigma)
  eps <- log(-log(runif(n)))   # standard Gumbel / extreme-value
  t_event <- exp(mu + sigma * eps)
  c_event <- rexp(n, rate = 0.03)
  time <- pmin(t_event, c_event)
  status <- as.integer(t_event <= c_event)
  data.frame(time = time, status = status, x1 = x1, x2 = x2)
}

## Quantile for the standardized baseline error used by the AFT families.
sim_aft_quantile <- function(p, dist = c("weibull", "logistic", "normal",
                                         "cauchy", "gamma", "exponential"),
                             dist_nu = 1) {
  dist <- match.arg(dist)
  if (dist == "exponential") {
    dist <- "weibull"
  }
  switch(dist,
    weibull = log(-log(1 - p)),
    logistic = log(p / (1 - p)),
    normal = stats::qnorm(p),
    cauchy = stats::qcauchy(p),
    gamma = {
      enu <- exp(dist_nu)
      (log(stats::qgamma(p, shape = enu, scale = 1)) - dist_nu) / enu
    }
  )
}

## Simulate from a Poisson-mixture cure model whose baseline is one of the
## package's standardized AFT families.
sim_cure_data <- function(n = 400,
                          dist = c("weibull", "logistic", "normal",
                                   "cauchy", "gamma", "exponential"),
                          theta1 = 0,
                          dist_nu = 1,
                          seed = 20260422,
                          censor_rate = 0.03) {
  dist <- match.arg(dist)
  set.seed(seed)

  x1 <- rnorm(n)
  x2 <- rbinom(n, 1, 0.4)
  mu <- 1.8 + 0.5 * x1 - 0.35 * x2
  if (dist == "exponential") {
    sigma <- rep(1, n)
  } else {
    log_sigma <- -0.25 + theta1 * (mu - mean(mu))
    sigma <- exp(log_sigma)
  }

  nu_lp <- qlogis(0.22) - 0.35 * x1 + 0.55 * x2
  cure_p <- stats::plogis(nu_lp)
  tau <- -log(cure_p)

  u <- stats::runif(n)
  has_event <- u <= (1 - cure_p)
  t_event <- rep(Inf, n)
  f0 <- -log1p(-u[has_event]) / tau[has_event]
  eps <- sim_aft_quantile(f0, dist = dist, dist_nu = dist_nu)
  t_event[has_event] <- exp(mu[has_event] + sigma[has_event] * eps)

  c_time <- stats::rexp(n, rate = censor_rate)
  time <- pmin(t_event, c_time)
  status <- as.integer(t_event <= c_time)

  data.frame(
    time = time,
    status = status,
    x1 = x1,
    x2 = x2,
    mu_true = mu,
    sigma_true = sigma,
    cure_true = cure_p
  )
}
