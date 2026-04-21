#' Fit an Accelerated Failure Time Model with Non-Proportional Hazards
#'
#' Fits the Anderson (1991) AFT regression model where both location (mu) and
#' scale (sigma) are functions of covariates, and sigma can depend on mu
#' through a theta polynomial.
#'
#' @param formula A formula of the form \code{Surv(time, status) ~ loc_vars | scale_vars}.
#'   Left of \code{|} specifies location (mu) covariates, right specifies
#'   scale (gamma) covariates. If no \code{|} is present, all covariates enter
#'   the location model and scale has intercept only.
#' @param data A data frame.
#' @param dist Character string specifying the error distribution:
#'   \code{"weibull"} (default), \code{"logistic"}, \code{"normal"},
#'   \code{"cauchy"}, or \code{"gamma"}.
#' @param theta Integer, order of the theta polynomial linking log(sigma) to mu.
#'   0 = no coupling (default), 1 = linear in mu, 2 = quadratic, etc.
#' @param theta_vars A one-sided formula specifying which location covariates
#'   form mu* (the mu that feeds into the theta polynomial for sigma).
#'   Default \code{NULL} means all location covariates.
#' @param nu Numeric, shape parameter for the gamma distribution (only used
#'   when \code{dist = "gamma"}).
#' @param init Numeric vector of initial parameter values, or \code{NULL}
#'   for zeros.
#' @param control A list of control parameters. See Details.
#' @param adjust Logical, whether to center covariates (default \code{TRUE}).
#' @param backend Character, which computation backend to use:
#'   \code{"c"} (default) or \code{"rust"}.
#'
#' @details
#' The model assumes
#' \deqn{\log(T) = \mu + \sigma \cdot \varepsilon}{log(T) = mu + sigma * epsilon}
#' where the location is a linear function of covariates,
#' \deqn{\mu = \beta_0 + \beta_1 x_1 + \cdots + \beta_p x_p}{mu = b0 + b1*x1 + ... + bp*xp}
#' and the log-dispersion can depend on covariates and on the location:
#' \deqn{\log(\sigma) = \gamma_0 + \gamma_1 z_1 + \cdots + \gamma_q z_q + \theta_1 \mu^{*} + \theta_2 \mu^{*2} + \cdots}{log(sigma) = g0 + g1*z1 + ... + gq*zq + theta1*mu* + theta2*mu*^2 + ...}
#' Here \eqn{\mu^{*} = \mu - \bar{\mu}}{mu* = mu - mean(mu)} is the centered location.
#' When \eqn{\theta_1 = 0}{theta1 = 0}, this reduces to the standard linear location AFT
#' model with proportional hazards under a Weibull distribution.
#' The error distribution \eqn{F} can be Weibull (extreme value), logistic,
#' normal, Cauchy, or gamma.
#'
#' Control parameters:
#' \describe{
#'   \item{acc}{Convergence accuracy (default 0.0001)}
#'   \item{maxiter}{Maximum Newton-Raphson iterations (default 50)}
#'   \item{maxhalv}{Maximum step halvings (default 20)}
#' }
#'
#' @return An object of class \code{"vldaft"}.
#'
#' @references
#' Anderson, K.M. (1991). A nonproportional hazards Weibull accelerated failure
#' time regression model. \emph{Biometrics}, 47, 281-288.
#'
#' @importFrom survival Surv
#' @importFrom stats model.frame model.matrix model.response update
#'   as.formula pnorm dnorm dgamma pgamma printCoefmat
#' @useDynLib vldaft, .registration = TRUE
#' @export
vldaft <- function(formula, data, dist = c("weibull", "logistic", "normal",
                                            "cauchy", "gamma"),
                   theta = 0L, theta_vars = NULL, nu = 1.0,
                   init = NULL, control = list(), adjust = TRUE,
                   backend = c("c", "rust")) {

  dist <- match.arg(dist)
  backend <- match.arg(backend)
  cl <- match.call()

  ## Parse formula: Surv(time, status) ~ loc | scale
  ## or: Surv(time, status) ~ loc (scale = intercept only)
  f_parts <- parse_vldaft_formula(formula, data)

  ## Build model matrices
  loc_mm <- model.matrix(f_parts$loc_formula, data)
  if (!is.null(f_parts$scale_formula)) {
    scale_mm <- model.matrix(f_parts$scale_formula, data)
  } else {
    scale_mm <- model.matrix(~ 1, data)
  }

  ## Get Surv response
  surv_obj <- f_parts$response
  if (!inherits(surv_obj, "Surv"))
    stop("Response must be a Surv() object")

  surv_type <- attr(surv_obj, "type")
  n <- nrow(surv_obj)

  ## Extract time, event, optional start time
  if (surv_type == "right") {
    time_var <- surv_obj[, "time"]
    status_var <- surv_obj[, "status"]
    start_var <- rep(0, n)
    t0_var <- rep(0, n)
    has_start <- FALSE
  } else if (surv_type == "counting") {
    t0_var <- surv_obj[, "start"]
    time_var <- surv_obj[, "stop"]
    status_var <- surv_obj[, "status"]
    start_var <- as.integer(t0_var > 0)
    has_start <- TRUE
  } else {
    stop("Only 'right' and 'counting' Surv types are supported")
  }

  ## Determine theta_vars (which location covariates form mu*)
  loc_names <- colnames(loc_mm)
  nloc <- ncol(loc_mm)
  nsc <- ncol(scale_mm)

  if (theta > 0 && !is.null(theta_vars)) {
    tv_mm <- model.matrix(theta_vars, data)
    tv_names <- colnames(tv_mm)
    ## theta_vars must be a subset of location vars
    if (!all(tv_names %in% loc_names))
      stop("theta_vars must be a subset of location covariates")
    ## Reorder location columns: non-theta-vars first (spmlo1), then theta-vars
    non_tv <- setdiff(loc_names, tv_names)
    loc_mm <- loc_mm[, c(non_tv, tv_names), drop = FALSE]
    loc_names <- colnames(loc_mm)
    mlo1_val <- length(non_tv)
  } else if (theta > 0) {
    ## All location vars couple to sigma via theta
    mlo1_val <- 0L
  } else {
    mlo1_val <- -1L  # signal: no override needed
  }

  ## Build combined data matrix: [scale_covs, loc_covs, time, event, start, t0]
  ncol_total <- nsc + nloc + 4  # +4 for time, event, start, t0
  data_matrix <- matrix(0, nrow = n, ncol = ncol_total)

  ## Scale covariates: columns 0..(nsc-1)
  data_matrix[, 1:nsc] <- scale_mm

  ## Location covariates: columns nsc..(nsc+nloc-1)
  data_matrix[, (nsc + 1):(nsc + nloc)] <- loc_mm

  ## Time: column nsc+nloc
  data_matrix[, nsc + nloc + 1] <- time_var

  ## Event: column nsc+nloc+1
  data_matrix[, nsc + nloc + 2] <- status_var

  ## Start indicator: column nsc+nloc+2
  data_matrix[, nsc + nloc + 3] <- start_var

  ## Start time: column nsc+nloc+3
  data_matrix[, nsc + nloc + 4] <- t0_var

  ## Column indices (0-based for C)
  scale_cols <- as.integer(0:(nsc - 1))
  loc_cols <- as.integer(nsc:(nsc + nloc - 1))
  time_col <- as.integer(nsc + nloc)
  event_col <- as.integer(nsc + nloc + 1)
  start_col_idx <- if (has_start) as.integer(nsc + nloc + 2) else -1L
  t0_col_idx <- if (has_start) as.integer(nsc + nloc + 3) else -1L

  ## Distribution code
  dist_code <- switch(dist,
    weibull = 1L, logistic = 2L, normal = 3L, cauchy = 4L, gamma = 5L)

  ## Control parameters
  ctrl <- list(acc = 0.0001, maxiter = 50L, maxhalv = 20L)
  ctrl[names(control)] <- control

  ## Call C or Rust backend
  call_name <- if (backend == "rust") "wrap__vldaft_fit_rust" else "vldaft_fit"
  fit <- .Call(call_name,
    data_matrix,
    time_col,
    event_col,
    as.integer(1L),       # event value = 1
    as.integer(-1L),      # right censoring value (not used; status=0 means censored)
    start_col_idx,
    t0_col_idx,
    loc_cols,
    scale_cols,
    as.integer(theta),
    as.integer(mlo1_val),
    dist_code,
    as.double(nu),
    if (!is.null(init)) as.double(init) else NULL,
    as.double(ctrl$acc),
    as.integer(ctrl$maxiter),
    as.integer(ctrl$maxhalv),
    as.integer(adjust)
  )

  ## Check convergence
  if (fit$iter < 0) {
    conv_msg <- switch(as.character(fit$iter),
      "-1" = "Unable to compute partials (large exponent)",
      "-2" = "Unable to compute partials",
      "-3" = "Unable to compute partials at solution",
      "-4" = paste("Failed to converge after", ctrl$maxiter, "iterations"),
      "-5" = "Information matrix not positive definite at solution",
      "-6" = "Unable to increase likelihood",
      "-7" = "Unable to compute likelihood",
      "Unknown error"
    )
    warning("vldaft did not converge: ", conv_msg)
  }

  ## Name the coefficients
  gamma_names <- paste0("gamma:", colnames(scale_mm))
  eta_names <- paste0("eta:", colnames(loc_mm))
  theta_names <- if (theta > 0) paste0("theta", 1:theta) else character(0)
  coef_names <- c(gamma_names, eta_names, theta_names)
  names(fit$coefficients) <- coef_names
  rownames(fit$vcov) <- coef_names
  colnames(fit$vcov) <- coef_names

  ## Build result object
  result <- list(
    coefficients = fit$coefficients,
    vcov = fit$vcov,
    loglik = fit$loglik,
    loglik_init = fit$loglik_init,
    score = fit$score,
    iter = fit$iter,
    npar = fit$npar,
    nobs = fit$nobs,
    dist = dist,
    theta = theta,
    theta_vars = theta_vars,
    formula = formula,
    call = cl,
    scale_names = colnames(scale_mm),
    loc_names = colnames(loc_mm),
    converged = fit$iter > 0
  )
  class(result) <- "vldaft"
  result
}


#' Parse vldaft formula
#'
#' Splits a formula like Surv(time, status) ~ x1 + x2 | x3 + x4
#' into location and scale parts.
#'
#' @param formula The formula.
#' @param data The data frame.
#' @return A list with response, loc_formula, scale_formula.
#' @keywords internal
parse_vldaft_formula <- function(formula, data) {
  ## Get response
  mf <- model.frame(update(formula, . ~ 1), data = data)
  response <- model.response(mf)

  ## Get RHS as character
  rhs <- as.character(formula)[3]

  ## Split on |
  if (grepl("\\|", rhs)) {
    parts <- strsplit(rhs, "\\|")[[1]]
    loc_rhs <- trimws(parts[1])
    scale_rhs <- trimws(parts[2])
    loc_formula <- as.formula(paste("~", loc_rhs))
    scale_formula <- as.formula(paste("~", scale_rhs))
  } else {
    loc_formula <- as.formula(paste("~", rhs))
    scale_formula <- NULL
  }

  list(response = response, loc_formula = loc_formula, scale_formula = scale_formula)
}
