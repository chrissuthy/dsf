#' Fit a destination selection function
#'
#' Fits a destination selection function (DSF) by maximum likelihood to
#' redistribution data in which each individual's origin, its destination,
#' and the full set of candidate destinations are known.
#'
#' @details
#' Individual \eqn{i} with origin \eqn{s_i} exerts, on each candidate
#' destination \eqn{j}, the hazard
#' \deqn{\lambda_{ij} = \exp\left(-\frac{d(s_i, v_j)^2}{2\sigma_i^2}\right)\,\theta_j}
#' and the destination at which it is observed is a categorical draw with
#' probabilities \eqn{\pi_{ij} = \lambda_{ij} / \sum_k \lambda_{ik}}, the sum
#' running over the destinations available to individual \eqn{i}. The
#' log-likelihood is \eqn{\sum_i \log \pi_{i, e_i}} where \eqn{e_i} is the
#' observed destination.
#'
#' Covariates enter through two linear predictors:
#' \eqn{\log \sigma_i = X_i \beta_\sigma} for individual-level covariates
#' (with intercept) and \eqn{\log \theta_j = Z_j \beta_\theta + \log(\mathrm{effort}_j)}
#' for destination-level covariates (no intercept, since a multiplier common
#' to all destinations cancels in the normalisation). Availability is defined
#' by `op_mat`: the destinations available to individual \eqn{i} are the
#' traps operational in `ind_df$obs_year[i]`.
#'
#' Parameters are estimated on the log scale with [stats::optim()] (BFGS);
#' standard errors are from the inverse numerical Hessian.
#'
#' @param sigma One-sided formula for \eqn{\log \sigma_i}; terms refer to
#'   columns of `ind_df`. Default `~ 1`.
#' @param theta One-sided formula for \eqn{\log \theta_j}; terms refer to
#'   columns of `trap_df`. Must have no intercept (`~ 0 + ...`). Default
#'   `~ 1`.
#' @param ind_df Data frame with one row per redistribution event. Required
#'   columns: `src_x`, `src_y` (origin coordinates), `dest` (integer row
#'   index of the destination in `trap_df`); `obs_year` if `op_mat` is
#'   supplied; plus any covariates named in `sigma`.
#' @param trap_df Data frame with one row per candidate destination. Required
#'   columns: `x`, `y`; optional `effort` (positive, used as an offset on
#'   \eqn{\log \theta_j}); plus any covariates named in `theta`.
#' @param op_mat Optional matrix of traps by year (rows in `trap_df` order,
#'   columns named by year), positive where the trap was operational. If
#'   `NULL`, every destination is available to every individual, with a
#'   message.
#' @param exclude_natal Logical; remove each individual's nearest trap from
#'   its choice set. Use when origins are trap locations and an individual
#'   cannot be recorded at its natal trap.
#' @param start Optional numeric vector of starting values, in the order
#'   `c(beta_sigma, beta_theta)`.
#' @param hessian Logical; compute standard errors?
#'
#' @return An object of class `"dsf"`: a list with
#' \describe{
#'   \item{`coef`, `se`}{Named estimates and standard errors on the log scale.
#'     Names are prefixed `sigma:` or `theta:`.}
#'   \item{`logLik`, `AIC`, `npar`, `convergence`}{Fit summaries;
#'     `convergence` is the [stats::optim()] code (0 = success).}
#'   \item{`sigma`, `theta`}{The formulas.}
#'   \item{`Xs`, `Xt`}{Design matrices.}
#'   \item{`avail`}{Logical n-by-J availability matrix used in the likelihood.}
#'   \item{`d2`}{n-by-J matrix of squared distances.}
#'   \item{`n`, `J`, `ind_df`, `trap_df`, `optim`}{Inputs and the raw
#'     optimiser output.}
#' }
#'
#' @seealso [simulate_dsf()] to generate data of the required form.
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
#' }
#'
#' @importFrom stats model.matrix optim median setNames
#' @export
dsf <- function(sigma = ~ 1, theta = ~ 1, ind_df, trap_df, op_mat = NULL,
                exclude_natal = FALSE, start = NULL, hessian = TRUE) {

  ## ------------------------------------------------------------ checks
  need_i <- c("src_x", "src_y", "dest")
  if (!all(need_i %in% names(ind_df)))
    stop("ind_df must contain: ", paste(need_i, collapse = ", "), call. = FALSE)
  if (!all(c("x", "y") %in% names(trap_df)))
    stop("trap_df must contain x and y", call. = FALSE)
  n <- nrow(ind_df); J <- nrow(trap_df)
  if (!all(ind_df$dest %in% seq_len(J)))
    stop("ind_df$dest must index rows of trap_df", call. = FALSE)

  ## ------------------------------------------------------------ design matrices
  Xs <- stats::model.matrix(sigma, ind_df)
  if (is_null_formula(theta)) {
    Xt <- matrix(0, J, 0)
  } else {
    Xt <- stats::model.matrix(theta, trap_df)
    if ("(Intercept)" %in% colnames(Xt))
      stop("theta must have no intercept (not identifiable): use ~ 0 + ...",
           call. = FALSE)
  }
  p <- ncol(Xs); q <- ncol(Xt)
  offset <- if (!is.null(trap_df$effort)) {
    if (any(trap_df$effort <= 0)) stop("effort must be positive", call. = FALSE)
    log(trap_df$effort)
  } else rep(0, J)

  ## ------------------------------------------------------------ distances & availability
  d2 <- outer(ind_df$src_x, trap_df$x, "-")^2 + outer(ind_df$src_y, trap_df$y, "-")^2

  if (is.null(op_mat)) {
    message("op_mat not supplied: assuming every destination available to every individual")
    avail <- matrix(TRUE, n, J)
  } else {
    if (is.null(ind_df$obs_year))
      stop("ind_df$obs_year is required when op_mat is supplied", call. = FALSE)
    yr <- as.character(ind_df$obs_year)
    if (nrow(op_mat) != J)
      stop("op_mat must have one row per trap", call. = FALSE)
    if (!all(yr %in% colnames(op_mat)))
      stop("obs_year values missing from colnames(op_mat): ",
           paste(setdiff(unique(yr), colnames(op_mat)), collapse = ", "), call. = FALSE)
    avail <- t(op_mat[, yr, drop = FALSE]) > 0
  }
  if (exclude_natal) avail[cbind(seq_len(n), apply(d2, 1, which.min))] <- FALSE

  obs <- cbind(seq_len(n), ind_df$dest)
  if (any(rowSums(avail) == 0))
    stop("individuals with no available destination: ",
         paste(which(rowSums(avail) == 0), collapse = ", "), call. = FALSE)
  if (!all(avail[obs]))
    stop("individuals observed at a destination marked unavailable: ",
         paste(which(!avail[obs]), collapse = ", "), call. = FALSE)
  mask <- ifelse(avail, 0, -Inf)

  ## ------------------------------------------------------------ negative log-likelihood
  nll <- function(par) {
    log_sigma <- drop(Xs %*% par[seq_len(p)])
    log_theta <- offset + if (q > 0) drop(Xt %*% par[p + seq_len(q)]) else 0
    ll <- -d2 / (2 * exp(2 * log_sigma)) +
          matrix(log_theta, n, J, byrow = TRUE) + mask
    mx <- apply(ll, 1, max)
    -sum(ll[obs] - (mx + log(rowSums(exp(ll - mx)))))
  }

  ## ------------------------------------------------------------ fit
  if (is.null(start))
    start <- c(log(stats::median(sqrt(d2[obs])) + 0.1), rep(0, p - 1 + q))
  if (length(start) != p + q)
    stop("start must have length ", p + q, call. = FALSE)

  fit <- stats::optim(start, nll, method = "BFGS", hessian = hessian,
                      control = list(maxit = 500, reltol = 1e-10))

  nms <- c(paste0("sigma:", colnames(Xs)),
           if (q > 0) paste0("theta:", colnames(Xt)))
  est <- stats::setNames(fit$par, nms)
  se  <- if (hessian) {
    tryCatch(stats::setNames(sqrt(diag(solve(fit$hessian))), nms),
             error = function(e) stats::setNames(rep(NA_real_, length(nms)), nms))
  } else NULL

  structure(list(coef = est, se = se,
                 logLik = -fit$value, npar = length(est),
                 AIC = 2 * fit$value + 2 * length(est),
                 convergence = fit$convergence,
                 sigma = sigma, theta = theta, Xs = Xs, Xt = Xt,
                 n = n, J = J, avail = avail, d2 = d2, obs = obs, offset = offset,
                 ind_df = ind_df, trap_df = trap_df, optim = fit),
            class = "dsf")
}

#' @describeIn dsf Print a fitted DSF.
#' @param x A `"dsf"` object.
#' @param ... Ignored.
#' @export
print.dsf <- function(x, ...) {
  cat("Destination selection function\n")
  cat("  sigma:", deparse(x$sigma), "   theta:", deparse(x$theta), "\n")
  cat("  n =", x$n, "  destinations =", x$J,
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
