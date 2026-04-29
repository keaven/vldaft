#' Fit an Accelerated Failure Time Model with Non-Proportional Hazards
#'
#' Fits the Anderson (1991) AFT regression model where both location (mu) and
#' scale (sigma) are functions of covariates, and sigma can depend on mu
#' through a theta polynomial.
#'
#' @param formula A formula of the form
#'   \code{Surv(time, status) ~ loc_vars | scale_vars}. The response may use
#'   right, left, counting-process, or interval-censored \code{Surv()} forms.
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
#' @param cure Logical; if \code{TRUE}, fit the Poisson-mixture cure-model
#'   extension with an additional cure-fraction regression block. When the
#'   formula contains a third \code{|}-separated block, \code{cure = TRUE}
#'   is implied.
#' @param init Numeric vector of initial parameter values, or \code{NULL}
#'   for zeros.
#' @param control A list of control parameters: \code{acc} (convergence
#'   tolerance, default \code{1e-4}), \code{maxiter} (max Newton-Raphson
#'   iterations, default 50), \code{maxhalv} (max step halvings, default 20).
#' @param adjust Logical; center covariates (default \code{TRUE}). Centering
#'   does not change location coefficients but does shift the interpretation
#'   of \eqn{\mu^*} used by the theta polynomial. The fitted object stores the
#'   training centering constants, reports coefficients back on the raw
#'   covariate scale, and reuses the original training centering inside
#'   \code{predict()}.
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
#' compare with \code{survival::survreg}, subtract
#' \eqn{\sum_i d_i \log(t_i)}{sum_i d_i log(t_i)}; \code{logLik()} exposes this
#' via its \code{jacobian} argument.
#'
#' The compiled back end has a compile-time cap of \code{MAXCOV = 30}
#' parameters (scale + location + theta). The R wrapper checks this before
#' dispatching. The C and Rust back ends support right, left, counting-process,
#' and interval-censored responses. Cure-model fits use BFGS with a numerical
#' Hessian around the compiled cure likelihood/score evaluator.
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
#'   setNames
#' @useDynLib vldaft, .registration = TRUE
#' @export
vldaft <- function(formula, data,
                   dist = c("weibull", "logistic", "normal", "cauchy", "gamma"),
                   theta = 0L, theta_vars = NULL, nu = 1.0, cure = FALSE,
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

  cure <- isTRUE(cure) || !is.null(f_parts$cure_formula)

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

  if (cure) {
    if (!is.null(f_parts$cure_formula)) {
      mf_cu <- model.frame(f_parts$cure_formula, data, na.action = stats::na.pass)
      if (anyNA(mf_cu))
        stop("NA values in cure covariates; please handle NAs before calling vldaft()",
             call. = FALSE)
      cure_mm <- model.matrix(f_parts$cure_formula, mf_cu)
    } else {
      cure_mm <- model.matrix(~ 1, data.frame(.y = rep(1, nrow(data))))
    }
  } else {
    cure_mm <- NULL
  }

  surv_obj <- f_parts$response
  if (!inherits(surv_obj, "Surv"))
    stop("Response must be a Surv() object", call. = FALSE)

  surv_type <- attr(surv_obj, "type")
  n <- nrow(surv_obj)

  if (surv_type == "right") {
    time_var <- surv_obj[, "time"]
    status_var <- surv_obj[, "status"]
    time2_var <- rep(NA_real_, n)
    start_var <- rep(0, n)
    t0_var <- rep(0, n)
    has_start <- FALSE
  } else if (surv_type == "left") {
    time_var <- surv_obj[, "time"]
    status_var <- ifelse(surv_obj[, "status"] == 1, 1L, -1L)
    time2_var <- rep(NA_real_, n)
    start_var <- rep(0, n)
    t0_var <- rep(0, n)
    has_start <- FALSE
  } else if (surv_type == "counting") {
    t0_var <- surv_obj[, "start"]
    time_var <- surv_obj[, "stop"]
    status_var <- surv_obj[, "status"]
    time2_var <- rep(NA_real_, n)
    start_var <- as.integer(t0_var > 0)
    has_start <- TRUE
  } else if (surv_type == "interval") {
    raw_status <- surv_obj[, "status"]
    time_var <- surv_obj[, "time1"]
    time2_var <- ifelse(raw_status == 3, surv_obj[, "time2"], NA_real_)
    status_var <- ifelse(raw_status == 1, 1L,
                         ifelse(raw_status == 2, -1L,
                                ifelse(raw_status == 3, 2L, 0L)))
    start_var <- rep(0, n)
    t0_var <- rep(0, n)
    has_start <- FALSE
  } else {
    stop("Only 'right', 'left', 'counting', and 'interval' Surv types are supported",
         call. = FALSE)
  }

  if (any(!is.finite(time_var) | time_var <= 0))
    stop("Survival times must be strictly positive", call. = FALSE)
  is_interval_row <- status_var == 2L
  if (any(is_interval_row & (!is.finite(time2_var) | time2_var <= time_var)))
    stop("Interval-censored upper bounds must be finite and greater than lower bounds",
         call. = FALSE)

  if (cure && backend != "c") {
    stop("Cure-model fitting is currently implemented only for backend = 'c'",
         call. = FALSE)
  }

  loc_names <- colnames(loc_mm)
  nloc <- ncol(loc_mm)
  nsc <- ncol(scale_mm)
  ncu <- if (cure) ncol(cure_mm) else 0L

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

  scale_means <- if (isTRUE(adjust)) .vldaft_block_means(scale_mm) else .vldaft_zero_block_means(scale_mm)
  loc_means   <- if (isTRUE(adjust)) .vldaft_block_means(loc_mm) else .vldaft_zero_block_means(loc_mm)
  cure_means  <- if (cure && isTRUE(adjust)) .vldaft_block_means(cure_mm) else .vldaft_zero_block_means(cure_mm)

  if (cure && isTRUE(adjust)) {
    loc_mm <- .vldaft_apply_training_centering(loc_mm, loc_means)
    scale_mm <- .vldaft_apply_training_centering(scale_mm, scale_means)
    if (cure) cure_mm <- .vldaft_apply_training_centering(cure_mm, cure_means)
  }

  total_par <- nsc + nloc + theta + ncu
  if (total_par > 30L)
    stop(sprintf(
      "vldaft supports at most 30 parameters (scale + location + cure + theta); got %d",
      total_par), call. = FALSE)

  if (cure && max(nloc, nsc, ncu, theta) > 30L) {
    stop("Each cure-model parameter block must have at most 30 coefficients",
         call. = FALSE)
  }

  ncol_total <- nsc + nloc + ncu + 5L
  data_matrix <- matrix(0, nrow = n, ncol = ncol_total)
  data_matrix[, seq_len(nsc)]                   <- scale_mm
  data_matrix[, nsc + seq_len(nloc)]            <- loc_mm
  if (cure) {
    data_matrix[, nsc + nloc + seq_len(ncu)]    <- cure_mm
  }
  data_matrix[, nsc + nloc + ncu + 1L]          <- time_var
  data_matrix[, nsc + nloc + ncu + 2L]          <- status_var
  data_matrix[, nsc + nloc + ncu + 3L]          <- start_var
  data_matrix[, nsc + nloc + ncu + 4L]          <- t0_var
  data_matrix[, nsc + nloc + ncu + 5L]          <- ifelse(is.na(time2_var), 0, time2_var)

  scale_cols <- as.integer(0:(nsc - 1))
  loc_cols   <- as.integer(nsc:(nsc + nloc - 1))
  cure_cols  <- if (cure) as.integer((nsc + nloc):(nsc + nloc + ncu - 1L)) else integer(0)
  time_col   <- as.integer(nsc + nloc + ncu)
  event_col  <- as.integer(nsc + nloc + ncu + 1L)
  start_col_idx <- if (has_start) as.integer(nsc + nloc + ncu + 2) else -1L
  t0_col_idx    <- if (has_start) as.integer(nsc + nloc + ncu + 3) else -1L
  time2_col      <- if (surv_type == "interval") as.integer(nsc + nloc + ncu + 4) else -1L

  dist_code <- switch(dist,
    weibull = 1L, logistic = 2L, normal = 3L, cauchy = 4L, gamma = 5L)

  ctrl <- list(acc = 1e-4, maxiter = 50L, maxhalv = 20L)
  ctrl[names(control)] <- control

  backend_used <- backend
  if (cure) {
    fit <- .fit_vldaft_cure(
      data_matrix = data_matrix,
      time_col = time_col,
      event_col = event_col,
      time2_col = time2_col,
      start_col_idx = start_col_idx,
      t0_col_idx = t0_col_idx,
      loc_cols = loc_cols,
      scale_cols = scale_cols,
      cure_cols = cure_cols,
      theta = theta,
      mlo1_val = mlo1_val,
      dist = dist,
      base_nu = nu,
      init = init,
      control = ctrl
    )
  } else {
    if (backend == "rust") {
      fit <- .Call("wrap__vldaft_fit_rust",
        data_matrix,
        time_col,
        time2_col,
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
    } else {
      fit <- .Call("vldaft_fit",
        data_matrix,
        time_col,
        time2_col,
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
    }
  }

  if (!is.null(fit$optim_convergence) && fit$optim_convergence != 0L) {
    conv_msg <- if (!is.null(fit$optim_message) && nzchar(fit$optim_message)) {
      fit$optim_message
    } else {
      paste("optimizer convergence code", fit$optim_convergence)
    }
    warning("vldaft optimizer did not converge: ", conv_msg, call. = FALSE)
  } else if (fit$iter < 0) {
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
  cure_names  <- if (cure) paste0("cure:", colnames(cure_mm)) else character(0)
  theta_names <- if (theta > 0) paste0("theta", seq_len(theta)) else character(0)
  coef_names  <- c(gamma_names, eta_names, cure_names, theta_names)
  names(fit$coefficients) <- coef_names
  rownames(fit$vcov) <- coef_names
  colnames(fit$vcov) <- coef_names

  fit$coefficients_internal <- fit$coefficients
  fit$vcov_internal <- fit$vcov
  transformed <- .vldaft_transform_parameters(
    coefficients = fit$coefficients_internal,
    vcov = fit$vcov_internal,
    scale_names = colnames(scale_mm),
    loc_names = colnames(loc_mm),
    cure_names = if (cure) colnames(cure_mm) else character(0),
    scale_means = scale_means,
    loc_means = loc_means,
    cure_means = cure_means,
    adjust = isTRUE(adjust)
  )
  fit$coefficients <- transformed$coefficients
  fit$vcov <- transformed$vcov

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
    mlo1         = mlo1_val,
    theta_vars   = theta_vars,
    nu           = nu,
    cure         = cure,
    adjust       = isTRUE(adjust),
    formula      = formula,
    call         = cl,
    backend      = backend_used,
    scale_names  = colnames(scale_mm),
    scale_means  = scale_means,
    loc_names    = colnames(loc_mm),
    loc_means    = loc_means,
    cure_names   = if (cure) colnames(cure_mm) else character(0),
    cure_means   = cure_means,
    coefficients_internal = fit$coefficients_internal,
    vcov_internal = fit$vcov_internal,
    time         = time_var,
    time2        = time2_var,
    status       = status_var,
    start        = t0_var,
    optim_convergence = if (!is.null(fit$optim_convergence)) fit$optim_convergence else NA_integer_,
    optim_message = if (!is.null(fit$optim_message)) fit$optim_message else NA_character_,
    hessian_method = if (!is.null(fit$hessian_method)) fit$hessian_method else NA_character_,
    converged    = if (!is.null(fit$converged)) isTRUE(fit$converged) else fit$iter > 0
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

  c(list(response = response), .vldaft_rhs_formulas(formula))
}

.vldaft_rhs_formulas <- function(formula) {
  if (!inherits(formula, "formula") || length(formula) != 3L)
    stop("`formula` must be a two-sided formula", call. = FALSE)

  rhs <- formula[[3L]]
  rhs_parts <- .flatten_vldaft_rhs(rhs)

  loc_formula <- stats::reformulate(deparse(rhs_parts[[1L]], width.cutoff = 500L),
                                    response = NULL)
  scale_formula <- if (length(rhs_parts) >= 2L) {
    stats::reformulate(deparse(rhs_parts[[2L]], width.cutoff = 500L),
                       response = NULL)
  } else NULL
  cure_formula <- if (length(rhs_parts) >= 3L) {
    stats::reformulate(deparse(rhs_parts[[3L]], width.cutoff = 500L),
                       response = NULL)
  } else NULL
  if (length(rhs_parts) > 3L) {
    stop("`formula` supports at most three RHS blocks: location | scale | cure",
         call. = FALSE)
  }

  list(loc_formula = loc_formula,
       scale_formula = scale_formula,
       cure_formula = cure_formula)
}

.flatten_vldaft_rhs <- function(rhs) {
  if (is.call(rhs) && identical(rhs[[1L]], as.name("|"))) {
    c(.flatten_vldaft_rhs(rhs[[2L]]), list(rhs[[3L]]))
  } else {
    list(rhs)
  }
}

.vldaft_zero_block_means <- function(mm) {
  if (is.null(mm)) return(stats::setNames(numeric(0), character(0)))
  stats::setNames(rep(0, ncol(mm)), colnames(mm))
}

.vldaft_block_means <- function(mm) {
  if (is.null(mm)) return(stats::setNames(numeric(0), character(0)))
  mm <- as.matrix(mm)
  means <- .vldaft_zero_block_means(mm)
  if (!ncol(mm)) return(means)
  idx <- seq_len(ncol(mm))
  intercept_col <- which(colnames(mm) == "(Intercept)")
  if (length(intercept_col) == 1L) idx <- setdiff(idx, intercept_col)
  if (length(idx)) {
    means[idx] <- colMeans(mm[, idx, drop = FALSE])
  }
  means
}

.vldaft_center_model_matrix <- function(mm) {
  mm <- .vldaft_apply_training_centering(mm, .vldaft_block_means(mm))
  mm
}

.vldaft_apply_training_centering <- function(mm, means) {
  if (is.null(mm)) return(mm)
  mm <- as.matrix(mm)
  if (!ncol(mm) || !length(means)) return(mm)
  idx_names <- setdiff(intersect(colnames(mm), names(means)), "(Intercept)")
  if (length(idx_names)) {
    mm[, idx_names] <- sweep(mm[, idx_names, drop = FALSE], 2, means[idx_names], FUN = "-")
  }
  mm
}

.vldaft_transform_parameters <- function(coefficients, vcov, scale_names, loc_names,
                                         cure_names, scale_means, loc_means,
                                         cure_means, adjust) {
  if (!isTRUE(adjust)) {
    return(list(coefficients = coefficients, vcov = vcov))
  }

  coef_names <- names(coefficients)
  Tmat <- diag(length(coefficients))
  rownames(Tmat) <- coef_names
  colnames(Tmat) <- coef_names

  add_block_transform <- function(prefix, block_names, means) {
    if (!length(block_names)) return()
    int_pos <- match("(Intercept)", block_names, nomatch = 0L)
    if (int_pos == 0L) return()
    idx <- match(paste0(prefix, block_names), coef_names)
    slope_pos <- setdiff(seq_along(block_names), int_pos)
    if (!length(slope_pos)) return()
    Tmat[idx[int_pos], idx[slope_pos]] <<- -unname(means[block_names[slope_pos]])
  }

  add_block_transform("gamma:", scale_names, scale_means)
  add_block_transform("eta:", loc_names, loc_means)
  add_block_transform("cure:", cure_names, cure_means)

  coef_raw <- as.numeric(Tmat %*% coefficients)
  names(coef_raw) <- coef_names
  vcov_raw <- Tmat %*% vcov %*% t(Tmat)
  rownames(vcov_raw) <- coef_names
  colnames(vcov_raw) <- coef_names

  list(coefficients = coef_raw, vcov = vcov_raw)
}

.vldaft_cure_dist_code <- function(dist) {
  switch(dist,
         weibull = 1L, logistic = 2L, normal = 3L, cauchy = 4L, gamma = 5L,
         stop("Unsupported cure-model distribution: ", dist, call. = FALSE))
}

.fit_vldaft_cure <- function(data_matrix, time_col, event_col, time2_col, start_col_idx,
                             t0_col_idx, loc_cols, scale_cols, cure_cols,
                             theta, mlo1_val, dist, base_nu, init, control) {
  nobs <- nrow(data_matrix)
  nloc <- length(loc_cols)
  nsc <- length(scale_cols)
  ncu <- length(cure_cols)
  dist_code <- .vldaft_cure_dist_code(dist)

  if (is.null(init)) {
    logt <- log(data_matrix[, time_col + 1L])
    p0 <- min(0.5, max(0.05, mean(data_matrix[, event_col + 1L] == 0) / 2))
    init <- c(
      rep(0, nsc),
      c(mean(logt), rep(0, max(0, nloc - 1L))),
      c(stats::qlogis(p0), rep(0, max(0, ncu - 1L))),
      rep(0, theta)
    )
    if (nsc > 0) init[1L] <- log(stats::sd(logt))
  } else {
    init <- as.double(init)
  }

  expected_len <- nsc + nloc + ncu + theta
  if (length(init) != expected_len) {
    stop(sprintf("`init` must have length %d for this cure model; got %d",
                 expected_len, length(init)),
         call. = FALSE)
  }

  cache <- new.env(parent = emptyenv())
  cache$par <- NULL
  cache$val <- NULL
  evaluate <- function(par) {
    if (!is.null(cache$par) && isTRUE(all.equal(par, cache$par, tolerance = 0))) {
      return(cache$val)
    }
    val <- .Call(
      "vldaft_cure_eval",
      data_matrix,
      as.integer(time_col),
      as.integer(event_col),
      as.integer(time2_col),
      as.integer(start_col_idx),
      as.integer(t0_col_idx),
      as.integer(loc_cols),
      as.integer(scale_cols),
      as.integer(cure_cols),
      as.integer(theta),
      as.integer(mlo1_val),
      as.integer(dist_code),
      as.double(base_nu),
      as.double(par),
      0L
    )
    cache$par <- par
    cache$val <- val
    val
  }
  safe_evaluate <- function(par) {
    tryCatch(
      evaluate(par),
      error = function(e) list(
        loglik = -1e12 - sum(par^2),
        gradient = rep(0, length(par))
      )
    )
  }

  objective <- function(par) {
    val <- safe_evaluate(par)
    -val$loglik
  }
  gradient <- function(par) {
    val <- safe_evaluate(par)
    -val$gradient
  }

  opt <- stats::optim(
    par = init,
    fn = objective,
    gr = gradient,
    method = "BFGS",
    control = list(reltol = control$acc, maxit = control$maxiter)
  )

  final_eval <- safe_evaluate(opt$par)
  init_eval <- safe_evaluate(init)
  hess <- stats::optimHess(opt$par, fn = objective, gr = gradient)
  vcov <- tryCatch(solve(hess), error = function(e) {
    warning("Unable to invert cure-model Hessian at solution", call. = FALSE)
    matrix(NA_real_, nrow = expected_len, ncol = expected_len)
  })

  list(
    coefficients = opt$par,
    vcov = vcov,
    loglik = final_eval$loglik,
    loglik_init = init_eval$loglik,
    score = sum(final_eval$gradient^2),
    iter = opt$counts[["function"]],
    npar = expected_len,
    nobs = nobs,
    optim_convergence = opt$convergence,
    optim_message = if (is.null(opt$message)) "" else opt$message,
    hessian_method = "optimHess",
    converged = opt$convergence == 0
  )
}

.vldaft_log1mexp <- function(logx) {
  safe_logx <- pmin(logx, 0)
  out <- log1p(-exp(safe_logx))
  use_expm1 <- safe_logx > -log(2)
  out[use_expm1] <- log(-expm1(safe_logx[use_expm1]))
  out[logx >= 0] <- -Inf
  out
}

.vldaft_logdiffexp <- function(log_hi, log_lo) {
  log_hi + .vldaft_log1mexp(pmin(log_lo - log_hi, 0))
}
