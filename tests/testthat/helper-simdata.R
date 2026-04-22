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
