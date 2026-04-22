#' Methods for vldaft Objects
#'
#' Print, summary, coefficient extraction, variance-covariance matrix,
#' log-likelihood, and prediction for \code{vldaft} model fits.
#'
#' @param x,object A \code{vldaft} object.
#' @param ... Additional arguments (mostly ignored; see specific methods).
#'
#' @return
#' \code{print} returns \code{x} invisibly.
#' \code{summary} returns a \code{summary.vldaft} object.
#' \code{coef} returns the coefficient vector.
#' \code{vcov} returns the variance-covariance matrix.
#' \code{logLik} returns a \code{logLik} object with \code{df} and
#'   \code{nobs} attributes.
#' \code{predict} returns a numeric vector (or matrix) on the requested scale.
#'
#' @name vldaft-methods
NULL


#' @rdname vldaft-methods
#' @param digits Number of significant digits for printed output.
#' @export
print.vldaft <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("AFT Regression Model (Anderson 1991)\n")
  cat("Distribution:", x$dist,
      "  Observations:", x$nobs,
      "  Events:", x$nevent, "\n")
  cat("Log-likelihood:", format(x$loglik, digits = digits),
      "  Iterations:", x$iter, "\n")
  if (x$theta > 0) cat("Theta polynomial order:", x$theta, "\n")

  scale_cf <- x$coefficients[grepl("^gamma:", names(x$coefficients))]
  loc_cf   <- x$coefficients[grepl("^eta:",   names(x$coefficients))]
  theta_cf <- x$coefficients[grepl("^theta[0-9]+$", names(x$coefficients))]

  cat("\nLocation (mu) coefficients [eta]:\n")
  print(loc_cf, digits = digits)
  cat("\nLog-dispersion (log sigma) coefficients [gamma]:\n")
  print(scale_cf, digits = digits)
  if (length(theta_cf)) {
    cat("\nLocation->dispersion coupling [theta]:\n")
    print(theta_cf, digits = digits)
  }
  if (!x$converged) cat("\n*** Model did not converge ***\n")
  invisible(x)
}


#' @rdname vldaft-methods
#' @export
summary.vldaft <- function(object, ...) {
  se <- sqrt(diag(object$vcov))
  z <- object$coefficients / se
  p <- 2 * stats::pnorm(-abs(z))
  tab <- cbind(Estimate = object$coefficients,
               `Std. Error` = se,
               `z value`   = z,
               `Pr(>|z|)`  = p)
  out <- list(
    call        = object$call,
    coefficients = tab,
    loglik      = object$loglik,
    loglik_init = object$loglik_init,
    score       = object$score,
    dist        = object$dist,
    theta       = object$theta,
    nobs        = object$nobs,
    nevent      = object$nevent,
    iter        = object$iter,
    converged   = object$converged
  )
  class(out) <- "summary.vldaft"
  out
}


#' @rdname vldaft-methods
#' @export
print.summary.vldaft <- function(x,
                                 digits = max(3L, getOption("digits") - 3L),
                                 ...) {
  cat("AFT Regression Model (Anderson 1991)\n\n")
  cat("Call:\n")
  print(x$call)
  cat("\nDistribution:", x$dist, "\n")
  cat("Observations:", x$nobs, "  Events:", x$nevent, "\n")
  cat("Log-likelihood:", format(x$loglik, digits = digits),
      " (init:", format(x$loglik_init, digits = digits), ")\n")
  cat("Score statistic:", format(x$score, digits = digits), "\n")
  cat("Iterations:", x$iter, "\n\n")
  cat("Coefficients:\n")
  stats::printCoefmat(x$coefficients, digits = digits,
                      P.values = TRUE, has.Pvalue = TRUE)
  if (!x$converged) cat("\n*** Model did not converge ***\n")
  invisible(x)
}


#' @rdname vldaft-methods
#' @export
coef.vldaft <- function(object, ...) object$coefficients


#' @rdname vldaft-methods
#' @export
vcov.vldaft <- function(object, ...) object$vcov


#' @rdname vldaft-methods
#' @param jacobian Logical (default \code{FALSE}). The C engine returns the
#'   log-likelihood on the standardized residual scale -- it does not include
#'   the \eqn{-\log t_i} Jacobian term that arises when changing variables
#'   from log-time to time. Set \code{jacobian = TRUE} to subtract
#'   \eqn{\sum_i d_i \log t_i}{sum(d * log(t))} so the returned value matches
#'   the time-scale convention used by \code{survival::survreg}.
#' @export
logLik.vldaft <- function(object, jacobian = FALSE, ...) {
  ll <- object$loglik
  if (isTRUE(jacobian)) {
    d <- object$status
    t <- object$time
    if (is.null(d) || is.null(t))
      stop("logLik(jacobian = TRUE) needs time/status stored on the fit; ",
           "refit with the current version of vldaft().", call. = FALSE)
    ll <- ll - sum(d * log(t))
  }
  attr(ll, "df") <- object$npar
  attr(ll, "nobs") <- object$nobs
  class(ll) <- "logLik"
  ll
}


#' @rdname vldaft-methods
#' @param newdata Optional data frame for \code{predict}. If omitted, the
#'   training data is implied via the stored \code{time}/\code{status};
#'   supply \code{newdata} to predict on new covariates.
#' @param type One of \code{"lp"} (linear predictor for mu, the default),
#'   \code{"logsigma"}, \code{"sigma"}, \code{"median"} (predicted log-time
#'   plus \code{sigma * q0.5} where \code{q0.5} is the median of the
#'   standardized error), or \code{"survival"} (requires \code{times}).
#' @param times Numeric vector of times at which to evaluate survival when
#'   \code{type = "survival"}.
#' @export
predict.vldaft <- function(object, newdata,
                           type = c("lp", "logsigma", "sigma", "median", "survival"),
                           times = NULL, ...) {
  type <- match.arg(type)
  if (missing(newdata))
    stop("`newdata` is required for predict.vldaft()", call. = FALSE)

  f <- parse_vldaft_formula(object$formula, newdata)
  loc_mm   <- stats::model.matrix(f$loc_formula, newdata)
  if (!is.null(f$scale_formula)) {
    scale_mm <- stats::model.matrix(f$scale_formula, newdata)
  } else {
    scale_mm <- stats::model.matrix(~ 1, data.frame(.y = rep(1, nrow(newdata))))
  }

  eta_cf   <- object$coefficients[grep("^eta:", names(object$coefficients))]
  gamma_cf <- object$coefficients[grep("^gamma:", names(object$coefficients))]
  theta_cf <- object$coefficients[grep("^theta[0-9]+$", names(object$coefficients))]

  loc_mm   <- loc_mm[, sub("^eta:", "", names(eta_cf)), drop = FALSE]
  scale_mm <- scale_mm[, sub("^gamma:", "", names(gamma_cf)), drop = FALSE]

  mu <- as.numeric(loc_mm %*% eta_cf)
  log_sigma <- as.numeric(scale_mm %*% gamma_cf)
  if (length(theta_cf)) {
    mu_c <- mu - mean(mu)
    for (k in seq_along(theta_cf)) {
      log_sigma <- log_sigma + theta_cf[k] * mu_c^k
    }
  }

  if (type == "lp") return(mu)
  if (type == "logsigma") return(log_sigma)
  if (type == "sigma") return(exp(log_sigma))
  if (type == "median") {
    q05 <- .vldaft_quantile(0.5, object$dist, object$nu)
    return(mu + exp(log_sigma) * q05)
  }

  if (is.null(times))
    stop("`times` is required when type = 'survival'", call. = FALSE)
  sigma <- exp(log_sigma)
  S <- vapply(times, function(t) {
    w <- (log(t) - mu) / sigma
    .vldaft_survival(w, object$dist, object$nu)
  }, numeric(length(mu)))
  if (length(mu) == 1) as.numeric(S) else S
}


## Internal: standardized residual quantile for the built-in distributions.
.vldaft_quantile <- function(p, dist, nu) {
  switch(dist,
    weibull  = log(-log(1 - p)),
    logistic = log(p / (1 - p)),
    normal   = stats::qnorm(p),
    cauchy   = stats::qcauchy(p),
    gamma    = {
      enu <- exp(nu)
      (log(stats::qgamma(p, shape = enu, scale = 1)) - nu) / enu
    },
    stop("Unsupported distribution: ", dist))
}

## Internal: survival function S(w) = P(epsilon > w).
.vldaft_survival <- function(w, dist, nu) {
  switch(dist,
    weibull  = exp(-exp(w)),
    logistic = 1 - stats::plogis(w),
    normal   = stats::pnorm(w, lower.tail = FALSE),
    cauchy   = stats::pcauchy(w, lower.tail = FALSE),
    gamma    = {
      enu <- exp(nu)
      stats::pgamma(exp(enu * w + nu), shape = enu, scale = 1,
                    lower.tail = FALSE)
    },
    stop("Unsupported distribution: ", dist))
}
