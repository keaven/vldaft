#' AFT Family Constructors for gamlss2
#'
#' These functions create gamlss2 family objects for accelerated failure time
#' models with five error distributions. Each family handles \code{Surv()}
#' responses with right censoring, left censoring, and left truncation.
#'
#' @param theta Integer, order of the theta polynomial linking log(sigma) to mu.
#'   Default 0 (no coupling).
#'
#' @details
#' The model is:
#' \deqn{\log(T) = \mu + \sigma_{eff} \cdot \varepsilon}
#' where \eqn{\log(\sigma_{eff}) = \log(\sigma) + \sum_k \theta_k (\mu - \bar\mu)^k}
#'
#' Parameters:
#' \describe{
#'   \item{mu}{Location parameter (identity link)}
#'   \item{sigma}{Scale parameter (log link)}
#'   \item{theta1, theta2, ...}{Coupling parameters (identity link), if theta > 0}
#' }
#'
#' @return A \code{gamlss2.family} object.
#'
#' @name AFT_families
NULL

#' @rdname AFT_families
#' @export
AFT_Weibull <- function(theta = 0L) {
  .aft_family("weibull", theta)
}

#' @rdname AFT_families
#' @export
AFT_Logistic <- function(theta = 0L) {
  .aft_family("logistic", theta)
}

#' @rdname AFT_families
#' @export
AFT_Normal <- function(theta = 0L) {
  .aft_family("normal", theta)
}

#' @rdname AFT_families
#' @export
AFT_Cauchy <- function(theta = 0L) {
  .aft_family("cauchy", theta)
}

#' @rdname AFT_families
#' @export
AFT_Gamma <- function(theta = 0L, nu = 1) {
  .aft_family("gamma", theta, nu = nu)
}


## Internal: build the family object
.aft_family <- function(dist, theta = 0L, nu = 1) {

  ## Parameter names and links
  par_names <- c("mu", "sigma")
  par_links <- c(mu = "identity", sigma = "log")
  if (theta > 0) {
    th_names <- paste0("theta", seq_len(theta))
    par_names <- c(par_names, th_names)
    lnks <- rep("identity", theta)
    names(lnks) <- th_names
    par_links <- c(par_links, lnks)
  }

  ## Density/survival/hazard helpers (on standardized scale w = (log(t) - mu) / sigma_eff)
  ## Returns list(logf, logS, logS_lower) for a single observation
  .log_densities <- function(w, dist, nu) {
    switch(dist,
      weibull = {
        logf <- w - exp(w)
        logS <- -exp(w)
      },
      logistic = {
        logf <- w - 2 * log(1 + exp(w))
        logS <- -log(1 + exp(w))
      },
      normal = {
        logf <- dnorm(w, log = TRUE)
        logS <- pnorm(w, lower.tail = FALSE, log.p = TRUE)
      },
      cauchy = {
        logf <- -log(pi) - log(1 + w^2)
        logS <- -atan(w) / pi + 0.5
        logS <- log(pmax(logS, .Machine$double.eps))  # logS is actually log(S)
      },
      gamma = {
        enu <- exp(nu)
        logf <- nu + enu * (nu + w) - exp(enu * w + nu) - lgamma(enu)
        ## S = 1 - pgamma(exp(enu*w+nu), enu, 1)
        logS <- pgamma(exp(enu * w + nu), enu, 1, lower.tail = FALSE, log.p = TRUE)
      }
    )
    list(logf = logf, logS = logS)
  }

  ## Score components: d(log f)/dw and d(log S)/dw
  .score_w <- function(w, dist, nu) {
    switch(dist,
      weibull = {
        ew <- exp(w)
        dlogf <- 1 - ew
        dlogS <- -ew
      },
      logistic = {
        ew <- exp(w)
        p <- 1 / (1 + ew)
        dlogf <- 1 - 2 * ew * p
        dlogS <- -ew * p
      },
      normal = {
        dlogf <- -w
        ## d(logS)/dw = -phi(w)/Phi(-w)
        dlogS <- -dnorm(w) / pnorm(w, lower.tail = FALSE)
      },
      cauchy = {
        dlogf <- -2 * w / (1 + w^2)
        dlogS <- -1 / (pi * (1 + w^2)) / pmax(0.5 - atan(w)/pi, .Machine$double.eps)
      },
      gamma = {
        enu <- exp(nu)
        u <- exp(enu * w + nu)
        dlogf <- enu * (1 - u)
        ## d(logS)/dw = -f_gamma(u) * du/dw / S
        ## du/dw = enu * u
        dlogS <- -dgamma(u, enu, 1) * enu * u /
                  pgamma(u, enu, 1, lower.tail = FALSE)
      }
    )
    list(dlogf = dlogf, dlogS = dlogS)
  }

  ## PDF function
  pdf <- function(y, par, log = FALSE) {
    mu <- par$mu
    sigma <- par$sigma
    log_sigma_eff <- log(sigma)
    if (theta > 0) {
      mu_c <- mu - mean(mu)
      for (k in seq_len(theta)) {
        log_sigma_eff <- log_sigma_eff + par[[paste0("theta", k)]] * mu_c^k
      }
    }
    sigma_eff <- exp(log_sigma_eff)

    ## Handle Surv objects
    if (inherits(y, "Surv")) {
      stype <- attr(y, "type")
      if (stype == "right") {
        logt <- log(y[, "time"])
        status <- y[, "status"]
        w <- (logt - mu) / sigma_eff
        dens <- .log_densities(w, dist, nu)
        ## uncensored: f(t) = f_w(w) / (t * sigma_eff)
        ## log-likelihood = logf(w) - log(t) - log(sigma_eff) for events
        ##                = logS(w) for right-censored
        ll <- ifelse(status == 1,
                     dens$logf - logt - log_sigma_eff,
                     dens$logS)
      } else if (stype == "counting") {
        t0 <- y[, "start"]
        logt <- log(y[, "stop"])
        status <- y[, "status"]
        w <- (logt - mu) / sigma_eff
        dens <- .log_densities(w, dist, nu)
        ll <- ifelse(status == 1,
                     dens$logf - logt - log_sigma_eff,
                     dens$logS)
        ## Subtract left-truncation contribution
        has_start <- t0 > 0
        if (any(has_start)) {
          w0 <- (log(t0[has_start]) - mu[has_start]) / sigma_eff[has_start]
          dens0 <- .log_densities(w0, dist, nu)
          ll[has_start] <- ll[has_start] - dens0$logS
        }
      } else {
        stop("Only 'right' and 'counting' Surv types supported")
      }
    } else {
      ## Non-Surv: treat as uncensored
      logt <- log(y)
      w <- (logt - mu) / sigma_eff
      dens <- .log_densities(w, dist, nu)
      ll <- dens$logf - logt - log_sigma_eff
    }

    if (log) ll else exp(ll)
  }

  ## Score functions
  score <- list()

  score$mu <- function(y, par, ...) {
    mu <- par$mu; sigma <- par$sigma
    log_se <- log(sigma)
    if (theta > 0) {
      mu_c <- mu - mean(mu)
      for (k in seq_len(theta)) log_se <- log_se + par[[paste0("theta", k)]] * mu_c^k
    }
    se <- exp(log_se)

    if (inherits(y, "Surv")) {
      stype <- attr(y, "type")
      logt <- if (stype == "right") log(y[, "time"]) else log(y[, "stop"])
      status <- if (stype == "right") y[, "status"] else y[, "status"]
    } else {
      logt <- log(y); status <- rep(1, length(y))
    }

    w <- (logt - mu) / se
    sc <- .score_w(w, dist, nu)

    ## d(ll)/d(mu) via chain rule
    ## w = (logt - mu)/se, dw/dmu = -1/se
    ## Also d(log_se)/dmu = sum_k theta_k * k * mu_c^(k-1)
    dlogse_dmu <- 0
    if (theta > 0) for (k in seq_len(theta)) dlogse_dmu <- dlogse_dmu + par[[paste0("theta", k)]] * k * mu_c^(k - 1)

    ## For event: d/dmu [logf(w) - log(se)] = dlogf/dw * dw/dmu - dlogse/dmu
    ##   dw/dmu = -1/se - w * dlogse_dmu
    dw_dmu <- -1/se - w * dlogse_dmu

    s <- ifelse(status == 1,
                sc$dlogf * dw_dmu - dlogse_dmu,
                sc$dlogS * dw_dmu)

    ## Left truncation adjustment
    if (inherits(y, "Surv") && attr(y, "type") == "counting") {
      t0 <- y[, "start"]; has_start <- t0 > 0
      if (any(has_start)) {
        w0 <- (log(t0[has_start]) - mu[has_start]) / se[has_start]
        sc0 <- .score_w(w0, dist, nu)
        dw0_dmu <- -1/se[has_start] - w0 * dlogse_dmu[has_start]
        s[has_start] <- s[has_start] - sc0$dlogS * dw0_dmu
      }
    }
    s
  }

  score$sigma <- function(y, par, ...) {
    mu <- par$mu; sigma <- par$sigma
    log_se <- log(sigma)
    if (theta > 0) {
      mu_c <- mu - mean(mu)
      for (k in seq_len(theta)) log_se <- log_se + par[[paste0("theta", k)]] * mu_c^k
    }
    se <- exp(log_se)

    if (inherits(y, "Surv")) {
      stype <- attr(y, "type")
      logt <- if (stype == "right") log(y[, "time"]) else log(y[, "stop"])
      status <- if (stype == "right") y[, "status"] else y[, "status"]
    } else {
      logt <- log(y); status <- rep(1, length(y))
    }

    w <- (logt - mu) / se
    sc <- .score_w(w, dist, nu)

    ## d(ll)/d(sigma) -- sigma has log link, so we need d/d(log(sigma))
    ## d(log_se)/d(log(sigma)) = 1
    ## dw/d(log(sigma)) = -w
    s <- ifelse(status == 1,
                sc$dlogf * (-w) - 1,
                sc$dlogS * (-w))

    ## Left truncation
    if (inherits(y, "Surv") && attr(y, "type") == "counting") {
      t0 <- y[, "start"]; has_start <- t0 > 0
      if (any(has_start)) {
        w0 <- (log(t0[has_start]) - mu[has_start]) / se[has_start]
        sc0 <- .score_w(w0, dist, nu)
        s[has_start] <- s[has_start] - sc0$dlogS * (-w0)
      }
    }
    s
  }

  ## Theta score functions
  if (theta > 0) {
    for (kk in seq_len(theta)) {
      local({
        k <- kk
        score[[paste0("theta", k)]] <<- function(y, par, ...) {
          mu <- par$mu; sigma <- par$sigma
          log_se <- log(sigma)
          mu_c <- mu - mean(mu)
          if (theta > 0) for (j in seq_len(theta)) log_se <- log_se + par[[paste0("theta", j)]] * mu_c^j
          se <- exp(log_se)

          if (inherits(y, "Surv")) {
            stype <- attr(y, "type")
            logt <- if (stype == "right") log(y[, "time"]) else log(y[, "stop"])
            status <- if (stype == "right") y[, "status"] else y[, "status"]
          } else {
            logt <- log(y); status <- rep(1, length(y))
          }

          w <- (logt - mu) / se
          sc <- .score_w(w, dist, nu)

          ## d(log_se)/d(theta_k) = mu_c^k
          ## dw/d(theta_k) = -w * mu_c^k
          mupow <- mu_c^k

          s <- ifelse(status == 1,
                      sc$dlogf * (-w * mupow) - mupow,
                      sc$dlogS * (-w * mupow))

          if (inherits(y, "Surv") && attr(y, "type") == "counting") {
            t0 <- y[, "start"]; has_start <- t0 > 0
            if (any(has_start)) {
              w0 <- (log(t0[has_start]) - mu[has_start]) / se[has_start]
              sc0 <- .score_w(w0, dist, nu)
              s[has_start] <- s[has_start] - sc0$dlogS * (-w0 * mupow[has_start])
            }
          }
          s
        }
      })
    }
  }

  ## Build family
  fam <- list(
    family = paste0("AFT_", dist),
    names = par_names,
    links = par_links,
    pdf = pdf,
    score = score
  )

  ## Initialize
  fam$initialize <- list(
    mu = function(y, ...) {
      if (inherits(y, "Surv")) {
        logt <- if (attr(y, "type") == "right") log(y[, "time"]) else log(y[, "stop"])
      } else logt <- log(y)
      rep(mean(logt), length(logt))
    },
    sigma = function(y, ...) {
      if (inherits(y, "Surv")) {
        logt <- if (attr(y, "type") == "right") log(y[, "time"]) else log(y[, "stop"])
      } else logt <- log(y)
      rep(sd(logt), length(logt))
    }
  )
  if (theta > 0) {
    for (k in seq_len(theta)) {
      fam$initialize[[paste0("theta", k)]] <- function(y, ...) {
        rep(0, NROW(y))
      }
    }
  }

  ## Valid response
  fam$valid.response <- function(x) {
    inherits(x, "Surv") || (is.numeric(x) && all(x > 0))
  }

  class(fam) <- "gamlss2.family"
  fam
}
