#' Methods for vldaft Objects
#'
#' Print, summary, coefficient extraction, variance-covariance matrix,
#' and log-likelihood for \code{vldaft} model fits.
#'
#' @param x,object A \code{vldaft} object.
#' @param ... Additional arguments (ignored).
#'
#' @return
#' \code{print} returns \code{x} invisibly.
#' \code{summary} returns a \code{summary.vldaft} object.
#' \code{coef} returns the coefficient vector.
#' \code{vcov} returns the variance-covariance matrix.
#' \code{logLik} returns a \code{logLik} object with \code{df} and \code{nobs} attributes.
#'
#' @name vldaft-methods
NULL

#' @rdname vldaft-methods
#' @export
print.vldaft <- function(x, ...) {
  cat("AFT Regression Model (Anderson 1991)\n")
  cat("Distribution:", x$dist, "\n")
  cat("Observations:", x$nobs, "\n")
  cat("Log-likelihood:", format(x$loglik, digits = 6), "\n")
  if (x$theta > 0) cat("Theta order:", x$theta, "\n")
  cat("\nCoefficients:\n")
  print(x$coefficients)
  if (!x$converged) cat("\n*** Model did not converge ***\n")
  invisible(x)
}

#' @rdname vldaft-methods
#' @export
summary.vldaft <- function(object, ...) {
  se <- sqrt(diag(object$vcov))
  z <- object$coefficients / se
  p <- 2 * pnorm(-abs(z))
  tab <- cbind(Estimate = object$coefficients, `Std. Error` = se,
               `z value` = z, `Pr(>|z|)` = p)
  out <- list(
    call = object$call,
    coefficients = tab,
    loglik = object$loglik,
    loglik_init = object$loglik_init,
    score = object$score,
    dist = object$dist,
    theta = object$theta,
    nobs = object$nobs,
    iter = object$iter,
    converged = object$converged
  )
  class(out) <- "summary.vldaft"
  out
}

#' @rdname vldaft-methods
#' @export
print.summary.vldaft <- function(x, ...) {
  cat("AFT Regression Model (Anderson 1991)\n\n")
  cat("Call:\n")
  print(x$call)
  cat("\nDistribution:", x$dist, "\n")
  cat("Observations:", x$nobs, "\n")
  cat("Log-likelihood:", format(x$loglik, digits = 6),
      " (init:", format(x$loglik_init, digits = 6), ")\n")
  cat("Score statistic:", format(x$score, digits = 6), "\n")
  cat("Iterations:", x$iter, "\n\n")
  cat("Coefficients:\n")
  printCoefmat(x$coefficients, P.values = TRUE, has.Pvalue = TRUE)
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
#' @export
logLik.vldaft <- function(object, ...) {
  ll <- object$loglik
  attr(ll, "df") <- object$npar
  attr(ll, "nobs") <- object$nobs
  class(ll) <- "logLik"
  ll
}
