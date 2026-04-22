#' Fit an Accelerated Failure Time Model with Non-Proportional Hazards
#'
#' Fits the Anderson (1991) AFT regression model where both location (mu) and
#' scale (sigma) are functions of covariates, and sigma can depend on mu
#' through a theta polynomial.
#'
#' @param formula A formula of the form
#'   \code{Surv(time, status) ~ loc_vars | scale_vars}.
#'   Left of \code{|} specifies location (mu) covariates, right specifies
#'   scale (gamma) covariates. If no \code{|} is present, all covariates
#'   enter the location model and the scale has an intercept only.
#' @param data A data frame.
#' @param dist Character: error distribution. One of \code{"weibull"}
#'   (default), \code{"logistic"}, \code{"normal"}, \code{"cauchy"},
#'   \code{"gamma"}.
#' @param theta Integer, order of the theta polynomial linking
#'   \eqn{\log(\sigma)} to \eqn{\mu}. \code{0} = no coupling (default),
#'   \code{1} = linear, \code{2} = quadratic, etc.
#' @param theta_vars A one-sided formula naming the location covariates that
#'   form \eqn{\mu^*} (the \eqn{\mu} that feeds into the theta polynomial
#'   for \eqn{\sigma}). Default \code{NULL} means all location covariates.
#' @param nu Numeric, shape parameter for the gamma distribution (only used
#'   when \code{dist = "gamma"}).
#' @param init Numeric vector of initial parameter values, or \code{NULL}
#'   for zeros.
#' @param control A list of control parameters: \code{acc} (convergence
#'   tolerance, default \code{1e-4}), \code{maxiter} (max Newton-Raphson
#'   iterations, default 50), \code{maxhalv} (max step halvings, default 20).
#' @param adjust Logical; center covariates (default \code{TRUE}). Centering
#'   does not change location coefficients but does shift the interpretation
#'   of \eqn{\mu^*} used by the theta polynomial.
#' @param backend Character, \code{"c"} (default) or \code{"rust"}.
#'
#' @details
#' The model assumes
#' \deqn{\log(T) = \mu + \sigma \cdot \varepsilon,}{log(T) = mu + sigma * epsilon,}
#' with
#' \deqn{\mu = \beta_0 + \beta_1 x_1 + \cdots + \beta_p x_p,}{mu = b0 + b1 x1 + ... + bp xp,}
#' \deqn{\log(\sigma) = \gamma_0 + \gamma_1 z_1 + \cdots + \gamma_q z_q + \theta_1 \mu^{*} + \theta_2 \mu^{*2} + \cdots,}{log(sigma) = g0 + g1 z1 + ... + gq zq + theta1 mustar + theta2 mustar^2 + ...,}
#' where \eqn{\mu^{*} = \mu - \bar{\mu}}{mu* = mu - mean(mu)} is the centered
#' location. Setting \eqn{\theta = 0}{theta = 0} recovers the standard linear
#' location AFT model (proportional hazards under the Weibull error). The
#' estimator is Newton-Raphson with modified Cholesky safeguards.
#'
#' Coefficient names are prefixed to mark their role:
#' \itemize{
#'   \item \code{gamma:} -- scale (dispersion) regression coefficients
#'   \item \code{eta:} -- location regression coefficients
#'   \item \code{theta1}, \code{theta2}, ... -- polynomial coupling terms
#' }
#'
#' The returned log-likelihood is on the log-time scale (no Jacobian). To
#' compare with \code{survival::survreg}, add
#' \eqn{\sum_i d_i \log(t_i)}{sum_i d_i log(t_i)}; \code{logLik()} exposes this
#' via its \code{jacobian} argument.
#'
#' The compiled back end has a compile-time cap of \code{MAXCOV = 30}
#' parameters (scale + location + theta). The R wrapper checks this before
#' dispatching.
#'
#' @return An object of class \code{"vldaft"}.
#'
#' @references
#' Anderson, K.M. (1991). A nonproportional hazards Weibull accelerated failure
#' time regression model. \emph{Biometrics}, 47, 281-288.
#'
#' @importFrom survival Surv
#' @importFrom stats model.frame model.matrix model.response update
#'   as.formula pnorm dnorm dgamma pgamma printCoefmat terms delete.response
#' @useDynLib vldaft, .registration = TRUE
#' @export
vldaft <- function(formula, data,
                   dist = c("weibull", "logistic", "normal", "cauchy", "gamma"),
                   theta = 0L, theta_vars = NULL, nu = 1.0,
                   init = NULL, control = list(), adjust = TRUE,
                   backend = c("c", "rust")) {

  dist <- match.arg(dist)
  backend <- match.arg(backend)
  cl <- match.call()

  if (missing(data))
    stop("`data` is required", call. = FALSE)
  if (!is.numeric(theta) || length(theta) != 1 || theta < 0 || theta != as.integer(theta))
    stop("`theta` must be a single non-negative integer", call. = FALSE)
  theta <- as.integer(theta)

  f_parts <- parse_vldaft_formula(formula, data)

  mf_loc <- model.frame(f_parts$loc_formula, data, na.action = stats::na.pass)
  if (anyNA(mf_loc))
    stop("NA values in location covariates; please handle NAs before calling vldaft()",
         call. = FALSE)
  loc_mm <- model.matrix(f_parts$loc_formula, mf_loc)

  if (!is.null(f_parts$scale_formula)) {
    mf_sc <- model.frame(f_parts$scale_formula, data, na.action = stats::na.pass)
    if (anyNA(mf_sc))
      stop("NA values in scale covariates; please handle NAs before calling vldaft()",
           call. = FALSE)
    scale_mm <- model.matrix(f_parts$scale_formula, mf_sc)
  } else {
    scale_mm <- model.matrix(~ 1, data.frame(.y = rep(1, nrow(data))))
  }

  surv_obj <- f_parts$response
  if (!inherits(surv_obj, "Surv"))
    stop("Response must be a Surv() object", call. = FALSE)

  surv_type <- attr(surv_obj, "type")
  n <- nrow(surv_obj)

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
    stop("Only 'right' and 'counting' Surv types are supported", call. = FALSE)
  }

  if (any(time_var <= 0))
    stop("Survival times must be strictly positive", call. = FALSE)

  loc_names <- colnames(loc_mm)
  nloc <- ncol(loc_mm)
  nsc <- ncol(scale_mm)

  if (theta > 0 && !is.null(theta_vars)) {
    tv_mm <- model.matrix(theta_vars, data)
    tv_names <- colnames(tv_mm)
    if (!all(tv_names %in% loc_names))
      stop("`theta_vars` must be a subset of location covariates", call. = FALSE)
    non_tv <- setdiff(loc_names, tv_names)
    loc_mm <- loc_mm[, c(non_tv, tv_names), drop = FALSE]
    loc_names <- colnames(loc_mm)
    mlo1_val <- length(non_tv)
  } else if (theta > 0) {
    mlo1_val <- 0L
  } else {
    mlo1_val <- -1L
  }

  total_par <- nsc + nloc + theta
  if (total_par > 30L)
    stop(sprintf(
      "vldaft supports at most 30 parameters (scale + location + theta); got %d",
      total_par), call. = FALSE)

  ncol_total <- nsc + nloc + 4L
  data_matrix <- matrix(0, nrow = n, ncol = ncol_total)
  data_matrix[, seq_len(nsc)]                   <- scale_mm
  data_matrix[, nsc + seq_len(nloc)]            <- loc_mm
  data_matrix[, nsc + nloc + 1L]                <- time_var
  data_matrix[, nsc + nloc + 2L]                <- status_var
  data_matrix[, nsc + nloc + 3L]                <- start_var
  data_matrix[, nsc + nloc + 4L]                <- t0_var

  scale_cols <- as.integer(0:(nsc - 1))
  loc_cols   <- as.integer(nsc:(nsc + nloc - 1))
  time_col   <- as.integer(nsc + nloc)
  event_col  <- as.integer(nsc + nloc + 1)
  start_col_idx <- if (has_start) as.integer(nsc + nloc + 2) else -1L
  t0_col_idx    <- if (has_start) as.integer(nsc + nloc + 3) else -1L

  dist_code <- switch(dist,
    weibull = 1L, logistic = 2L, normal = 3L, cauchy = 4L, gamma = 5L)

  ctrl <- list(acc = 1e-4, maxiter = 50L, maxhalv = 20L)
  ctrl[names(control)] <- control

  call_name <- if (backend == "rust") "wrap__vldaft_fit_rust" else "vldaft_fit"
  fit <- .Call(call_name,
    data_matrix,
    time_col,
    event_col,
    1L,
    -1L,
    start_col_idx,
    t0_col_idx,
    loc_cols,
    scale_cols,
    theta,
    as.integer(mlo1_val),
    dist_code,
    as.double(nu),
    if (!is.null(init)) as.double(init) else NULL,
    as.double(ctrl$acc),
    as.integer(ctrl$maxiter),
    as.integer(ctrl$maxhalv),
    as.integer(isTRUE(adjust))
  )

  if (fit$iter < 0) {
    conv_msg <- switch(as.character(fit$iter),
      "-1" = "Unable to compute partials (large exponent)",
      "-2" = "Unable to compute partials",
      "-3" = "Unable to compute partials at solution",
      "-4" = paste("Failed to converge after", ctrl$maxiter, "iterations"),
      "-5" = "Information matrix not positive-definite at solution",
      "-6" = "Unable to increase likelihood",
      "-7" = "Unable to compute likelihood",
      "Unknown error"
    )
    warning("vldaft did not converge: ", conv_msg, call. = FALSE)
  }

  gamma_names <- paste0("gamma:", colnames(scale_mm))
  eta_names   <- paste0("eta:", colnames(loc_mm))
  theta_names <- if (theta > 0) paste0("theta", seq_len(theta)) else character(0)
  coef_names  <- c(gamma_names, eta_names, theta_names)
  names(fit$coefficients) <- coef_names
  rownames(fit$vcov) <- coef_names
  colnames(fit$vcov) <- coef_names

  nevent <- sum(status_var == 1)

  result <- list(
    coefficients = fit$coefficients,
    vcov         = fit$vcov,
    loglik       = fit$loglik,
    loglik_init  = fit$loglik_init,
    score        = fit$score,
    iter         = fit$iter,
    npar         = fit$npar,
    nobs         = fit$nobs,
    nevent       = nevent,
    dist         = dist,
    theta        = theta,
    theta_vars   = theta_vars,
    nu           = nu,
    adjust       = isTRUE(adjust),
    formula      = formula,
    call         = cl,
    backend      = backend,
    scale_names  = colnames(scale_mm),
    loc_names    = colnames(loc_mm),
    time         = time_var,
    status       = status_var,
    start        = t0_var,
    converged    = fit$iter > 0
  )
  class(result) <- "vldaft"
  result
}


#' Parse a vldaft model formula
#'
#' Splits a formula like \code{Surv(time, status) ~ x1 + x2 | x3 + x4}
#' into its location and scale components. Uses a term-based split so it
#' is robust to multiline formulas and to expressions on either side.
#'
#' @param formula The model formula.
#' @param data The data frame (used only to evaluate the response).
#' @return A list with elements \code{response}, \code{loc_formula},
#'   and \code{scale_formula} (possibly \code{NULL}).
#' @keywords internal
parse_vldaft_formula <- function(formula, data) {
  if (!inherits(formula, "formula") || length(formula) != 3L)
    stop("`formula` must be a two-sided formula", call. = FALSE)

  mf <- model.frame(update(formula, . ~ 1), data = data)
  response <- model.response(mf)

  rhs <- formula[[3L]]
  if (length(rhs) >= 2L && identical(rhs[[1L]], as.name("|"))) {
    loc_formula <- stats::reformulate(deparse(rhs[[2L]], width.cutoff = 500L),
                                      response = NULL)
    scale_formula <- stats::reformulate(deparse(rhs[[3L]], width.cutoff = 500L),
                                        response = NULL)
  } else {
    loc_formula <- stats::reformulate(deparse(rhs, width.cutoff = 500L),
                                      response = NULL)
    scale_formula <- NULL
  }

  list(response = response,
       loc_formula = loc_formula,
       scale_formula = scale_formula)
}
