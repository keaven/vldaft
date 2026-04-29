#' AFT and Cure Family Constructors for gamlss2
#'
#' These functions create \code{gamlss2} family objects for accelerated
#' failure time models and their Poisson-mixture cure-model extensions.
#' They handle \code{Surv()} responses with right censoring, left censoring,
#' interval censoring, and left truncation (counting-process form).
#'
#' @param theta Integer, order of the theta polynomial linking
#'   \eqn{\log(\sigma)} to \eqn{\mu}. Default 0 (no coupling).
#' @param nu Numeric, shape parameter for the gamma distribution
#'   (only used by \code{AFT_Gamma}; the parameterization mirrors
#'   the C/Rust engines, which treat \code{enu = exp(nu)} as the
#'   gamma shape). Default \code{1}.
#' @param base_nu Numeric, baseline gamma-shape parameter used by
#'   \code{Cure_Gamma()}. This is separate from the cure-fraction
#'   regression parameter \code{nu}. Default \code{1}.
#' @importFrom stats sd
#'
#' @details
#' The plain AFT model is
#' \deqn{\log(T) = \mu + \sigma_{eff} \cdot \varepsilon}
#' where
#' \deqn{\log(\sigma_{eff}) = \log(\sigma) + \sum_k \theta_k (\mu - \bar\mu)^k.}
#'
#' The cure-model extension wraps the baseline survival distribution
#' \eqn{S_0(t)} as
#' \deqn{S(t) = \exp\{-\tau (1 - S_0(t))\}, \qquad \tau = -\log(p),}
#' where \eqn{p = \mathrm{logit}^{-1}(\nu)} is the cure probability.
#'
#' Parameters:
#' \describe{
#'   \item{mu}{Location parameter (identity link).}
#'   \item{sigma}{Scale parameter (log link), omitted for
#'     \code{Cure_Exponential()}.}
#'   \item{nu}{Cure probability \eqn{p} on the response scale
#'     (logit link) for the cure-model families.}
#'   \item{theta1, theta2, ...}{Coupling parameters (identity link), if
#'     \code{theta > 0}.}
#' }
#'
#' @return A \code{gamlss2.family} object.
#'
#' @name AFT_families
NULL

#' @rdname AFT_families
#' @export
AFT_Weibull <- function(theta = 0L) {
  .aft_family("weibull", theta = theta)
}

#' @rdname AFT_families
#' @export
AFT_Logistic <- function(theta = 0L) {
  .aft_family("logistic", theta = theta)
}

#' @rdname AFT_families
#' @export
AFT_Normal <- function(theta = 0L) {
  .aft_family("normal", theta = theta)
}

#' @rdname AFT_families
#' @export
AFT_Cauchy <- function(theta = 0L) {
  .aft_family("cauchy", theta = theta)
}

#' @rdname AFT_families
#' @export
AFT_Gamma <- function(theta = 0L, nu = 1) {
  .aft_family("gamma", theta = theta, dist_nu = nu)
}

#' @rdname AFT_families
#' @export
Cure_Weibull <- function(theta = 0L) {
  .aft_family("weibull", theta = theta, cure = TRUE)
}

#' @rdname AFT_families
#' @export
Cure_Logistic <- function(theta = 0L) {
  .aft_family("logistic", theta = theta, cure = TRUE)
}

#' @rdname AFT_families
#' @export
Cure_Normal <- function(theta = 0L) {
  .aft_family("normal", theta = theta, cure = TRUE)
}

#' @rdname AFT_families
#' @export
Cure_Cauchy <- function(theta = 0L) {
  .aft_family("cauchy", theta = theta, cure = TRUE)
}

#' @rdname AFT_families
#' @export
Cure_Gamma <- function(theta = 0L, base_nu = 1) {
  .aft_family("gamma", theta = theta, cure = TRUE, dist_nu = base_nu)
}

#' @rdname AFT_families
#' @export
Cure_Exponential <- function() {
  .aft_family("weibull", theta = 0L, cure = TRUE, fixed_sigma = TRUE)
}


## Baseline helpers ---------------------------------------------------------

.aft_family <- function(dist, theta = 0L, dist_nu = 1, cure = FALSE,
                        fixed_sigma = FALSE) {
  theta <- as.integer(theta)
  if (theta < 0L) {
    stop("`theta` must be non-negative", call. = FALSE)
  }
  if (fixed_sigma && theta > 0L) {
    stop("`theta` is not supported when sigma is fixed", call. = FALSE)
  }

  par_names <- "mu"
  par_links <- c(mu = "identity")
  if (!fixed_sigma) {
    par_names <- c(par_names, "sigma")
    par_links <- c(par_links, sigma = "log")
  }
  if (cure) {
    par_names <- c(par_names, "nu")
    par_links <- c(par_links, nu = "logit")
  }
  if (theta > 0L) {
    th_names <- paste0("theta", seq_len(theta))
    par_names <- c(par_names, th_names)
    th_links <- rep("identity", theta)
    names(th_links) <- th_names
    par_links <- c(par_links, th_links)
  }

  pdf <- function(y, par, log = FALSE) {
    state <- .aft_family_state(y, par, dist = dist, theta = theta,
                               dist_nu = dist_nu, cure = cure,
                               fixed_sigma = fixed_sigma)
    ll <- state$ll
    if (log) ll else exp(ll)
  }

  score <- list()
  score$mu <- function(y, par, ...) {
    .aft_score_parameter(y, par, target = "mu", dist = dist, theta = theta,
                         dist_nu = dist_nu, cure = cure,
                         fixed_sigma = fixed_sigma)
  }

  if (!fixed_sigma) {
    score$sigma <- function(y, par, ...) {
      .aft_score_parameter(y, par, target = "sigma", dist = dist,
                           theta = theta, dist_nu = dist_nu, cure = cure,
                           fixed_sigma = fixed_sigma)
    }
  }

  if (cure) {
    score$nu <- function(y, par, ...) {
      .aft_score_parameter(y, par, target = "nu", dist = dist, theta = theta,
                           dist_nu = dist_nu, cure = TRUE,
                           fixed_sigma = fixed_sigma)
    }
  }

  if (theta > 0L) {
    for (kk in seq_len(theta)) {
      local({
        k <- kk
        name <- paste0("theta", k)
        score[[name]] <<- function(y, par, ...) {
          .aft_score_parameter(y, par, target = name, dist = dist,
                               theta = theta, dist_nu = dist_nu,
                               cure = cure, fixed_sigma = fixed_sigma)
        }
      })
    }
  }

  fam <- list(
    family = .aft_family_name(dist, cure = cure, fixed_sigma = fixed_sigma),
    names = par_names,
    links = par_links,
    pdf = pdf,
    score = score
  )

  fam$initialize <- list(
    mu = function(y, ...) {
      resp <- .aft_response_parts(y)
      rep(mean(resp$logt), length(resp$logt))
    }
  )
  if (!fixed_sigma) {
    fam$initialize$sigma <- function(y, ...) {
      resp <- .aft_response_parts(y)
      rep(stats::sd(resp$logt), length(resp$logt))
    }
  }
  if (cure) {
    fam$initialize$nu <- function(y, ...) {
      resp <- .aft_response_parts(y)
      censor_frac <- mean(resp$status == 0)
      p0 <- min(0.5, max(0.05, censor_frac / 2))
      rep(p0, length(resp$logt))
    }
  }
  if (theta > 0L) {
    for (k in seq_len(theta)) {
      fam$initialize[[paste0("theta", k)]] <- function(y, ...) {
        rep(0, NROW(y))
      }
    }
  }

  fam$valid.response <- function(x) {
    inherits(x, "Surv") || (is.numeric(x) && all(x > 0))
  }

  fam$cure_probability <- function(par) {
    if (!cure) {
      stop("Cure probability is only defined for cure-model families",
           call. = FALSE)
    }
    .aft_cure_probability(par$nu)
  }

  fam$survival <- function(times, par) {
    par_grid <- .aft_expand_parameters(par, length(times))
    state <- .aft_family_state(rep(times, each = length(par$mu)), par_grid,
                               dist = dist, theta = theta, dist_nu = dist_nu,
                               cure = cure, fixed_sigma = fixed_sigma,
                               y_is_time = TRUE)
    S <- matrix(state$survival, nrow = length(par$mu), byrow = FALSE)
    if (length(par$mu) == 1L) as.numeric(S) else S
  }

  class(fam) <- "gamlss2.family"
  fam
}

.aft_family_name <- function(dist, cure = FALSE, fixed_sigma = FALSE) {
  if (fixed_sigma) {
    return(if (cure) "Cure_Exponential" else "AFT_Exponential")
  }
  dist_label <- switch(
    dist,
    weibull = "Weibull",
    logistic = "Logistic",
    normal = "Normal",
    cauchy = "Cauchy",
    gamma = "Gamma",
    stop("Unsupported distribution: ", dist, call. = FALSE)
  )
  paste0(if (cure) "Cure_" else "AFT_", dist_label)
}

.aft_response_parts <- function(y) {
  if (inherits(y, "Surv")) {
    stype <- attr(y, "type")
    if (stype == "right") {
      list(
        logt = log(y[, "time"]),
        logt2 = rep(NA_real_, nrow(y)),
        status = y[, "status"],
        t0 = rep(0, nrow(y)),
        has_start = rep(FALSE, nrow(y)),
        is_counting = FALSE
      )
    } else if (stype == "left") {
      list(
        logt = log(y[, "time"]),
        logt2 = rep(NA_real_, nrow(y)),
        status = ifelse(y[, "status"] == 1, 1L, -1L),
        t0 = rep(0, nrow(y)),
        has_start = rep(FALSE, nrow(y)),
        is_counting = FALSE
      )
    } else if (stype == "counting") {
      t0 <- y[, "start"]
      list(
        logt = log(y[, "stop"]),
        logt2 = rep(NA_real_, nrow(y)),
        status = y[, "status"],
        t0 = t0,
        has_start = t0 > 0,
        is_counting = TRUE
      )
    } else if (stype == "interval") {
      raw_status <- y[, "status"]
      list(
        logt = log(y[, "time1"]),
        logt2 = ifelse(raw_status == 3, log(y[, "time2"]), NA_real_),
        status = ifelse(raw_status == 1, 1L,
                        ifelse(raw_status == 2, -1L,
                               ifelse(raw_status == 3, 2L, 0L))),
        t0 = rep(0, nrow(y)),
        has_start = rep(FALSE, nrow(y)),
        is_counting = FALSE
      )
    } else {
      stop("Only 'right', 'left', 'counting', and 'interval' Surv types are supported",
           call. = FALSE)
    }
  } else {
    list(
      logt = log(y),
      logt2 = rep(NA_real_, length(y)),
      status = rep(1, length(y)),
      t0 = rep(0, length(y)),
      has_start = rep(FALSE, length(y)),
      is_counting = FALSE
    )
  }
}

.aft_effective_scale <- function(mu, par, theta, fixed_sigma = FALSE) {
  n <- length(mu)
  if (fixed_sigma) {
    return(list(
      log_sigma_eff = rep(0, n),
      sigma_eff = rep(1, n),
      mu_c = rep(0, n)
    ))
  }

  log_sigma_eff <- log(par$sigma)
  mu_c <- mu - mean(mu)
  if (theta > 0L) {
    for (k in seq_len(theta)) {
      log_sigma_eff <- log_sigma_eff + par[[paste0("theta", k)]] * mu_c^k
    }
  }

  list(
    log_sigma_eff = log_sigma_eff,
    sigma_eff = exp(log_sigma_eff),
    mu_c = mu_c
  )
}

.aft_baseline_quantities <- function(w, dist, dist_nu) {
  out <- switch(dist,
    weibull = {
      ew <- exp(w)
      list(
        logf = w - ew,
        logS = -ew,
        dlogf = 1 - ew,
        dlogS = -ew
      )
    },
    logistic = {
      ew <- exp(w)
      p <- 1 / (1 + ew)
      list(
        logf = w - 2 * log1p(ew),
        logS = -log1p(ew),
        dlogf = 1 - 2 * ew * p,
        dlogS = -ew * p
      )
    },
    normal = {
      list(
        logf = stats::dnorm(w, log = TRUE),
        logS = stats::pnorm(w, lower.tail = FALSE, log.p = TRUE),
        dlogf = -w,
        dlogS = -stats::dnorm(w) / stats::pnorm(w, lower.tail = FALSE)
      )
    },
    cauchy = {
      surv <- pmax(0.5 - atan(w) / pi, .Machine$double.eps)
      list(
        logf = -log(pi) - log1p(w^2),
        logS = log(surv),
        dlogf = -2 * w / (1 + w^2),
        dlogS = -1 / (pi * (1 + w^2) * surv)
      )
    },
    gamma = {
      enu <- exp(dist_nu)
      u <- exp(enu * w + dist_nu)
      surv <- stats::pgamma(u, enu, 1, lower.tail = FALSE)
      list(
        logf = dist_nu + enu * (dist_nu + w) - u - lgamma(enu),
        logS = stats::pgamma(u, enu, 1, lower.tail = FALSE, log.p = TRUE),
        dlogf = enu * (1 - u),
        dlogS = -stats::dgamma(u, enu, 1) * enu * u / surv
      )
    },
    stop("Unsupported distribution: ", dist, call. = FALSE)
  )

  out$S <- exp(out$logS)
  out$F <- 1 - out$S
  out
}

.aft_parameter_derivatives <- function(target, mu, sigma_eff, w, mu_c, par,
                                       theta, fixed_sigma = FALSE) {
  if (identical(target, "mu")) {
    dlogse <- rep(0, length(mu))
    if (!fixed_sigma && theta > 0L) {
      for (k in seq_len(theta)) {
        dlogse <- dlogse + par[[paste0("theta", k)]] * k * mu_c^(k - 1)
      }
    }
    return(list(dw = -1 / sigma_eff - w * dlogse, dlogse = dlogse))
  }

  if (identical(target, "sigma")) {
    if (fixed_sigma) {
      stop("sigma is fixed for this family", call. = FALSE)
    }
    return(list(dw = -w, dlogse = rep(1, length(mu))))
  }

  if (startsWith(target, "theta")) {
    k <- as.integer(sub("^theta", "", target))
    mupow <- mu_c^k
    return(list(dw = -w * mupow, dlogse = mupow))
  }

  stop("Unsupported target parameter: ", target, call. = FALSE)
}

.aft_cure_probability <- function(p) {
  pmin(pmax(p, .Machine$double.eps), 1 - .Machine$double.eps)
}

.aft_cure_tau <- function(p) {
  -log(.aft_cure_probability(p))
}

.aft_index_parameter <- function(x, idx) {
  n <- if (is.logical(idx)) sum(idx) else length(idx)
  if (length(x) == 1L) {
    rep(x, n)
  } else {
    x[idx]
  }
}

.aft_expand_parameters <- function(par, times) {
  lapply(par, rep, times = times)
}

.aft_family_state <- function(y, par, dist, theta, dist_nu, cure,
                              fixed_sigma = FALSE, y_is_time = FALSE) {
  mu <- par$mu
  scale_state <- .aft_effective_scale(mu, par, theta, fixed_sigma = fixed_sigma)

  if (y_is_time) {
    resp <- list(
      logt = log(y),
      logt2 = rep(NA_real_, length(y)),
      status = rep(0, length(y)),
      t0 = rep(0, length(y)),
      has_start = rep(FALSE, length(y)),
      is_counting = FALSE
    )
  } else {
    resp <- .aft_response_parts(y)
  }

  w <- (resp$logt - mu) / scale_state$sigma_eff
  base <- .aft_baseline_quantities(w, dist = dist, dist_nu = dist_nu)
  base_event_ll <- base$logf - resp$logt - scale_state$log_sigma_eff
  exact <- resp$status == 1L
  right <- resp$status == 0L
  left <- resp$status == -1L
  interval <- resp$status == 2L
  base2 <- NULL

  if (cure) {
    p <- .aft_cure_probability(par$nu)
    tau <- .aft_cure_tau(p)
    log_surv <- -tau * base$F
    ll <- numeric(length(resp$status))
    ll[exact] <- log(tau[exact]) + base_event_ll[exact] - tau[exact] * base$F[exact]
    ll[right] <- log_surv[right]
    ll[left] <- .vldaft_log1mexp(log_surv[left])
    if (any(interval)) {
      w2 <- (resp$logt2[interval] - mu[interval]) / scale_state$sigma_eff[interval]
      base2 <- .aft_baseline_quantities(w2, dist = dist, dist_nu = dist_nu)
      log_surv2 <- -tau[interval] * base2$F
      ll[interval] <- .vldaft_logdiffexp(log_surv[interval], log_surv2)
    }
    survival <- exp(log_surv)
  } else {
    ll <- numeric(length(resp$status))
    ll[exact] <- base_event_ll[exact]
    ll[right] <- base$logS[right]
    ll[left] <- .vldaft_log1mexp(base$logS[left])
    if (any(interval)) {
      w2 <- (resp$logt2[interval] - mu[interval]) / scale_state$sigma_eff[interval]
      base2 <- .aft_baseline_quantities(w2, dist = dist, dist_nu = dist_nu)
      ll[interval] <- .vldaft_logdiffexp(base$logS[interval], base2$logS)
    }
    survival <- base$S
  }

  if (resp$is_counting && any(resp$has_start)) {
    w0 <- (log(resp$t0[resp$has_start]) - mu[resp$has_start]) /
      scale_state$sigma_eff[resp$has_start]
    base0 <- .aft_baseline_quantities(w0, dist = dist, dist_nu = dist_nu)
    if (cure) {
      tau <- .aft_cure_tau(par$nu[resp$has_start])
      ll[resp$has_start] <- ll[resp$has_start] + tau * base0$F
    } else {
      ll[resp$has_start] <- ll[resp$has_start] - base0$logS
    }
  }

  list(
    ll = ll,
    survival = survival,
    mu = mu,
    sigma_eff = scale_state$sigma_eff,
    log_sigma_eff = scale_state$log_sigma_eff,
    mu_c = scale_state$mu_c,
    resp = resp,
    base = base,
    base2 = base2
  )
}

.aft_score_parameter <- function(y, par, target, dist, theta, dist_nu, cure,
                                 fixed_sigma = FALSE) {
  state <- .aft_family_state(y, par, dist = dist, theta = theta,
                             dist_nu = dist_nu, cure = cure,
                             fixed_sigma = fixed_sigma)

  if (identical(target, "nu")) {
    p <- .aft_cure_probability(par$nu)
    tau <- .aft_cure_tau(p)
    exact <- state$resp$status == 1L
    right <- state$resp$status == 0L
    left <- state$resp$status == -1L
    interval <- state$resp$status == 2L
    log_surv <- -tau * state$base$F
    cure_surv <- exp(log_surv)
    s <- numeric(length(state$resp$status))
    s[exact] <- -(1 - p[exact]) / tau[exact] + (1 - p[exact]) * state$base$F[exact]
    s[right] <- (1 - p[right]) * state$base$F[right]
    s[left] <- -cure_surv[left] * (1 - p[left]) * state$base$F[left] /
      pmax(1 - cure_surv[left], .Machine$double.eps)
    if (any(interval)) {
      p_i <- p[interval]
      tau_i <- tau[interval]
      log_surv1 <- log_surv[interval]
      cure_surv1 <- exp(log_surv1)
      cure_surv2 <- exp(-tau_i * state$base2$F)
      denom <- pmax(cure_surv1 - cure_surv2, .Machine$double.eps)
      s[interval] <- (cure_surv1 * (1 - p_i) * state$base$F[interval] -
                        cure_surv2 * (1 - p_i) * state$base2$F) / denom
    }
    if (state$resp$is_counting && any(state$resp$has_start)) {
      has_start <- state$resp$has_start
      w0 <- (log(state$resp$t0[has_start]) - state$mu[has_start]) /
        state$sigma_eff[has_start]
      base0 <- .aft_baseline_quantities(w0, dist = dist, dist_nu = dist_nu)
      s[has_start] <- s[has_start] - (1 - p[has_start]) * base0$F
    }
    return(s)
  }

  deriv <- .aft_parameter_derivatives(
    target = target,
    mu = state$mu,
    sigma_eff = state$sigma_eff,
    w = (state$resp$logt - state$mu) / state$sigma_eff,
    mu_c = state$mu_c,
    par = par,
    theta = theta,
    fixed_sigma = fixed_sigma
  )

  event_score <- state$base$dlogf * deriv$dw - deriv$dlogse
  surv_score <- state$base$dlogS * deriv$dw
  surv_deriv <- state$base$S * surv_score
  exact <- state$resp$status == 1L
  right <- state$resp$status == 0L
  left <- state$resp$status == -1L
  interval <- state$resp$status == 2L

  if (cure) {
    tau <- .aft_cure_tau(par$nu)
    log_surv <- -tau * state$base$F
    cure_surv <- exp(log_surv)
    s <- numeric(length(state$resp$status))
    s[exact] <- event_score[exact] + tau[exact] * surv_deriv[exact]
    s[right] <- tau[right] * surv_deriv[right]
    s[left] <- -cure_surv[left] * tau[left] * surv_deriv[left] /
      pmax(1 - cure_surv[left], .Machine$double.eps)
    if (any(interval)) {
      w2 <- (state$resp$logt2[interval] - state$mu[interval]) /
        state$sigma_eff[interval]
      deriv2 <- .aft_parameter_derivatives(
        target = target,
        mu = state$mu[interval],
        sigma_eff = state$sigma_eff[interval],
        w = w2,
        mu_c = state$mu_c[interval],
        par = lapply(par, .aft_index_parameter, idx = interval),
        theta = theta,
        fixed_sigma = fixed_sigma
      )
      surv_score2 <- state$base2$dlogS * deriv2$dw
      surv_deriv2 <- state$base2$S * surv_score2
      cure_surv1 <- cure_surv[interval]
      cure_surv2 <- exp(-tau[interval] * state$base2$F)
      denom <- pmax(cure_surv1 - cure_surv2, .Machine$double.eps)
      s[interval] <- (cure_surv1 * tau[interval] * surv_deriv[interval] -
                        cure_surv2 * tau[interval] * surv_deriv2) / denom
    }
  } else {
    s <- numeric(length(state$resp$status))
    s[exact] <- event_score[exact]
    s[right] <- surv_score[right]
    s[left] <- -surv_deriv[left] / pmax(state$base$F[left], .Machine$double.eps)
    if (any(interval)) {
      w2 <- (state$resp$logt2[interval] - state$mu[interval]) /
        state$sigma_eff[interval]
      deriv2 <- .aft_parameter_derivatives(
        target = target,
        mu = state$mu[interval],
        sigma_eff = state$sigma_eff[interval],
        w = w2,
        mu_c = state$mu_c[interval],
        par = lapply(par, .aft_index_parameter, idx = interval),
        theta = theta,
        fixed_sigma = fixed_sigma
      )
      surv_score2 <- state$base2$dlogS * deriv2$dw
      surv_deriv2 <- state$base2$S * surv_score2
      denom <- pmax(state$base$S[interval] - state$base2$S, .Machine$double.eps)
      s[interval] <- (surv_deriv[interval] - surv_deriv2) / denom
    }
  }

  if (state$resp$is_counting && any(state$resp$has_start)) {
    has_start <- state$resp$has_start
    w0 <- (log(state$resp$t0[has_start]) - state$mu[has_start]) /
      state$sigma_eff[has_start]
    base0 <- .aft_baseline_quantities(w0, dist = dist, dist_nu = dist_nu)
    deriv0 <- .aft_parameter_derivatives(
      target = target,
      mu = state$mu[has_start],
      sigma_eff = state$sigma_eff[has_start],
      w = w0,
      mu_c = state$mu_c[has_start],
      par = lapply(par, .aft_index_parameter, idx = has_start),
      theta = theta,
      fixed_sigma = fixed_sigma
    )
    surv_score0 <- base0$dlogS * deriv0$dw
    if (cure) {
      tau0 <- .aft_cure_tau(par$nu[has_start])
      s[has_start] <- s[has_start] - tau0 * base0$S * surv_score0
    } else {
      s[has_start] <- s[has_start] - surv_score0
    }
  }

  s
}
