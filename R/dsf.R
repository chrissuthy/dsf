#' Fit a destination selection function
#'
#' Fits a destination selection function (DSF) by maximum likelihood to
#' redistribution data in which each individual's origin, its destination,
#' and the full set of candidate destinations are known.
#'
#' @details
#' Individual \eqn{i}, observed in year \eqn{t_i}, exerts on each candidate
#' destination \eqn{j} the hazard
#' \deqn{\lambda_{ij} = \exp\left(-\frac{d(s_i, v_j)^2}{2\sigma_i^2}\right)\,\theta_{j t_i}}
#' and the destination at which it is observed is a categorical draw with
#' probabilities \eqn{\pi_{ij} = \lambda_{ij} / \sum_k \lambda_{ik}}, the sum
#' running over the destinations operational in year \eqn{t_i}. The
#' log-likelihood is \eqn{\sum_i \log \pi_{i, e_i}} where \eqn{e_i} is the
#' observed destination.
#'
#' Covariates enter through two linear predictors:
#' \eqn{\log \sigma_i = X_i \beta_\sigma} for individual-level covariates
#' (with intercept) and
#' \eqn{\log \theta_{jt} = Z_{jt} \beta_\theta + \log(\mathrm{effort}_{jt})}
#' for destination-level covariates, which may vary by year. Because the
#' normalisation is within individual and each individual belongs to one
#' year, a term constant across destinations within a year (an intercept, a
#' year effect) cancels and is refused.
#'
#' With `cost` and `landscape` supplied, Euclidean distance is replaced by the
#' least-cost (ecological) distance through a resistance surface whose cell
#' cost is \eqn{\exp(\alpha^\top c(x))}, and \eqn{\alpha} is estimated
#' alongside the other parameters (coefficients named `cost:`). Paths are
#' computed on the landscape aggregated to `cost_res`; see [make_cost_fun()].
#' \eqn{\alpha = 0} recovers (grid-approximate) Euclidean distance, so a
#' likelihood-ratio test against the Euclidean fit tests the landscape effect.
#'
#' Parameters are estimated on the log scale with [stats::optim()] (BFGS);
#' standard errors are from the inverse numerical Hessian.
#'
#' @param sigma One-sided formula for \eqn{\log \sigma_i}; terms refer to
#'   columns of the individual table. Default `~ 1`.
#' @param theta One-sided formula for \eqn{\log \theta_{jt}}; terms refer to
#'   columns of the destination table and of `trap_tv`. Must have no
#'   intercept (`~ 0 + ...`). Default `~ 1`.
#' @param data A `"dsf_frame"` from [dsf_frame()], or the individual data
#'   frame (then `trap_df`, `op_mat`, `trap_tv` and `exclude_natal` are
#'   passed to [dsf_frame()]).
#' @param trap_df,op_mat,trap_tv,exclude_natal Passed to [dsf_frame()] when
#'   `data` is a data frame; ignored when it is a frame.
#' @param cost Optional one-sided formula naming layers of `landscape` that
#'   define the resistance surface, e.g. `~ cov`. `NULL` (default) fits the
#'   Euclidean model.
#' @param landscape `SpatRaster` with the layers named in `cost`.
#' @param cost_res Cell size for least-cost path calculations; the landscape
#'   is aggregated to this resolution. Coarser is faster.
#' @param start Optional numeric vector of starting values, in the order
#'   `c(beta_sigma, beta_theta, alpha)`.
#' @param hessian Logical; compute standard errors?
#'
#' @return An object of class `"dsf"`: a list with
#' \describe{
#'   \item{`coef`, `se`}{Named estimates and standard errors on the log scale.
#'     Names are prefixed `sigma:`, `theta:` or `cost:`.}
#'   \item{`logLik`, `AIC`, `npar`, `convergence`}{Fit summaries;
#'     `convergence` is the [stats::optim()] code (0 = success).}
#'   \item{`sigma`, `theta`, `cost`}{The formulas.}
#'   \item{`Xs`, `Xt`}{Design matrices (`Xt` has one row per
#'     destination-year).}
#'   \item{`frame`}{The `"dsf_frame"` the model was fitted to.}
#'   \item{`d2`}{n-by-J matrix of squared distances at the estimate.}
#'   \item{`cost_fun`}{For the non-Euclidean model, the distance calculator
#'     from [make_cost_fun()]; otherwise `NULL`.}
#'   \item{`n`, `J`, `optim`}{Sizes and the raw optimiser output.}
#' }
#'
#' @seealso [dsf_frame()], [simulate_dsf()]
#'
#' @examples
#' \donttest{
#' sim <- simulate_dsf(~ age_c, ~ 0 + cov,
#'                     beta_sigma = c(log(5), 0.3), beta_theta = 0.8,
#'                     surface = TRUE, exclude_natal = TRUE, seed = 1)
#' f0 <- dsf(~ 1,     ~ 1,       sim$ind_df, sim$trap_df, sim$op_mat, exclude_natal = TRUE)
#' f1 <- dsf(~ age_c, ~ 0 + cov, sim$ind_df, sim$trap_df, sim$op_mat, exclude_natal = TRUE)
#' f1
#' 2 * (f1$logLik - f0$logLik)   # likelihood-ratio statistic, 2 df
#'
#' # the same through a frame
#' fr <- dsf_frame(sim$ind_df, sim$trap_df, sim$op_mat, exclude_natal = TRUE)
#' dsf(~ age_c, ~ 0 + cov, fr)
#' }
#'
#' @importFrom stats model.matrix optim median setNames
#' @export
dsf <- function(sigma = ~ 1, theta = ~ 1, data, trap_df = NULL, op_mat = NULL,
                trap_tv = NULL, exclude_natal = FALSE,
                cost = NULL, landscape = NULL, cost_res = 1,
                start = NULL, hessian = TRUE) {

  ## ------------------------------------------------------------ frame
  if (inherits(data, "dsf_frame")) {
    fr <- data
  } else {
    if (is.null(trap_df)) stop("trap_df is required when data is not a dsf_frame", call. = FALSE)
    if (is.null(op_mat))
      message("op_mat not supplied: assuming every destination available in every year")
    fr <- dsf_frame(data, trap_df, op_mat, trap_tv, exclude_natal)
  }
  if (is.null(fr$obs)) stop("ind$dest is required for fitting", call. = FALSE)
  n <- fr$n; J <- fr$J

  ## ------------------------------------------------------------ design matrices
  Xs <- stats::model.matrix(sigma, fr$ind)
  Xt <- theta_design(theta, fr)
  p <- ncol(Xs); q <- ncol(Xt)
  off_rows <- theta_rows(as.vector(fr$offset), fr)      # n x J log effort

  ## ------------------------------------------------------------ distances
  d2 <- fr$d2
  cost_fun <- NULL; k <- 0
  if (!is.null(cost)) {
    if (is.null(landscape)) stop("landscape is required when cost is supplied", call. = FALSE)
    cost_fun <- make_cost_fun(landscape, cost, cost_res,
                              fr$ind[, c("src_x", "src_y")], fr$traps[, c("x", "y")])
    k <- length(cost_fun$names)
  }
  mask <- ifelse(fr$avail, 0, -Inf)
  obs  <- fr$obs

  ## ------------------------------------------------------------ negative log-likelihood
  nll <- function(par) {
    log_sigma <- drop(Xs %*% par[seq_len(p)])
    lt <- off_rows
    if (q > 0) lt <- lt + theta_rows(drop(Xt %*% par[p + seq_len(q)]), fr)
    dd <- if (k > 0) cost_fun$fun(par[p + q + seq_len(k)]) else d2
    ll <- -dd / (2 * exp(2 * log_sigma)) + lt + mask
    mx <- apply(ll, 1, max)
    -sum(ll[obs] - (mx + log(rowSums(exp(ll - mx)))))
  }

  ## ------------------------------------------------------------ fit
  if (is.null(start))
    start <- c(log(stats::median(sqrt(d2[obs])) + 0.1), rep(0, p - 1 + q + k))
  if (length(start) != p + q + k)
    stop("start must have length ", p + q + k, call. = FALSE)

  fit <- stats::optim(start, nll, method = "BFGS", hessian = hessian,
                      control = list(maxit = 500, reltol = 1e-10))

  nms <- c(paste0("sigma:", colnames(Xs)),
           if (q > 0) paste0("theta:", colnames(Xt)),
           if (k > 0) paste0("cost:", cost_fun$names))
  if (k > 0) d2 <- cost_fun$fun(fit$par[p + q + seq_len(k)])
  est <- stats::setNames(fit$par, nms)
  se  <- if (hessian) {
    tryCatch(stats::setNames(sqrt(diag(solve(fit$hessian))), nms),
             error = function(e) stats::setNames(rep(NA_real_, length(nms)), nms))
  } else NULL

  structure(list(coef = est, se = se,
                 logLik = -fit$value, npar = length(est),
                 AIC = 2 * fit$value + 2 * length(est),
                 convergence = fit$convergence,
                 sigma = sigma, theta = theta, cost = cost, Xs = Xs, Xt = Xt,
                 n = n, J = J, frame = fr, avail = fr$avail, d2 = d2, obs = obs,
                 cost_fun = cost_fun, cost_res = cost_res,
                 ind_df = fr$ind, trap_df = fr$traps, optim = fit),
            class = "dsf")
}

#' @describeIn dsf Print a fitted DSF.
#' @param x A `"dsf"` object.
#' @param ... Ignored.
#' @export
print.dsf <- function(x, ...) {
  cat("Destination selection function\n")
  cat("  sigma:", deparse(x$sigma), "   theta:", deparse(x$theta),
      if (!is.null(x$cost)) paste("   cost:", deparse(x$cost)) else "", "\n")
  cat("  n =", x$n, "  destinations =", x$J, "  years =", x$frame$T,
      "  logLik =", round(x$logLik, 2), "  AIC =", round(x$AIC, 2),
      if (x$convergence != 0) "  [optim did not converge]" else "", "\n\n")
  tab <- cbind(Estimate = x$coef)
  if (!is.null(x$se))
    tab <- cbind(tab, SE = x$se,
                 lower = x$coef - 1.96 * x$se, upper = x$coef + 1.96 * x$se)
  print(round(tab, 4))
  i <- match("sigma:(Intercept)", names(x$coef))
  if (!is.na(i))
    cat("\n  sigma at reference level:", round(exp(x$coef[i]), 3), "\n")
  invisible(x)
}
