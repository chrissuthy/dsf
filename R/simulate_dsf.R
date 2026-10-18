#' Simulate redistribution data from a destination selection function
#'
#' Generates a landscape, a trap array, a set of individuals with known
#' origins, and the trap in which each individual was captured, by drawing
#' from a destination selection function (DSF). The function returns only
#' what a field programme would record (individuals, traps, and the
#' trap-by-year operational matrix), so that its output can be passed
#' directly to [dsf()].
#'
#' @details
#' Individual \eqn{i} with origin \eqn{s_i} and destination \eqn{j} has hazard
#' \deqn{\lambda_{ij} = \exp\left(-\frac{d(s_i, v_j)^2}{2\sigma_i^2}\right)\,\theta_j}
#' with \eqn{\log \sigma_i = X_i \beta_\sigma} (individual covariates, with
#' intercept) and \eqn{\log \theta_j = Z_j \beta_\theta} (destination
#' covariates, no intercept). The destination is a single categorical draw
#' from the hazards normalised over the traps operational in the individual's
#' observation year.
#'
#' The landscape covariate is a standardised Matérn Gaussian random field
#' simulated by circulant embedding ([fields::sim.rf()]); `smoothness = 0.5`
#' gives the exponential covariance. It is attached to traps as `cov` and to
#' origins as `src_cov`, so `theta_formula = ~ 0 + cov` and
#' `sigma_formula = ~ src_cov` work without further processing. The raster is
#' returned for use as a cost surface in the non-Euclidean DSF.
#'
#' If `design` is supplied, the landscape, traps and individuals are taken
#' from it and all design arguments are ignored. This allows the same design
#' to be reused with different parameters or across replicate dispersal draws.
#'
#' @param sigma_formula One-sided formula for the individual-level linear
#'   predictor of \eqn{\log \sigma_i}. Terms refer to columns of the
#'   individual table (`age`, `age_c`, `src_cov`, ...). Default `~ 1`.
#' @param theta_formula One-sided formula for the destination-level linear
#'   predictor of \eqn{\log \theta_j}. Must have no intercept (use
#'   `~ 0 + cov`); a term constant across destinations is not identifiable.
#'   Default `~ 1` (all \eqn{\theta_j = 1}).
#' @param beta_sigma Numeric vector of coefficients for `sigma_formula`, on
#'   the log scale, in the order of `colnames(model.matrix(sigma_formula))`.
#'   For `~ 1`, a single value `log(sigma)`.
#' @param beta_theta Numeric vector of coefficients for `theta_formula`, on
#'   the log scale. `NULL` when `theta_formula = ~ 1`.
#' @param ntraps Number of traps.
#' @param nyears Number of birth cohorts (years 1 to `nyears`).
#' @param nborn Number of individuals born per cohort.
#' @param prage Probabilities of age at observation 1, 2, 3 years.
#' @param L Side of the square landscape (km); origin at (0, 0).
#' @param ncell Number of raster cells per side.
#' @param buffer Traps and uniform births are placed at least this far (km)
#'   from the landscape edge.
#' @param surface Logical; simulate a landscape covariate? If `FALSE` the
#'   covariate is zero everywhere.
#' @param aRange Range parameter of the Matérn covariance, in km (the same
#'   units as `L`), so the patch size does not change with `ncell`.
#' @param smoothness Matérn smoothness; 0.5 is exponential.
#' @param births `"traps"` places each origin at a randomly chosen trap
#'   (the parent's capture location); `"uniform"` draws origins uniformly
#'   inside the buffer.
#' @param op_mat Optional trap-by-year operational matrix (0/1), with one row
#'   per trap and columns named by year. If `NULL`, every trap is open in
#'   every observation year.
#' @param exclude_natal Logical; remove each individual's nearest trap
#'   (distance zero when `births = "traps"`) from its choice set?
#' @param design Optional design list (`landscape`, `traps`, `animals`) from a
#'   previous call, returned as `$design`.
#' @param seed Optional integer passed to [set.seed()].
#'
#' @return A list with components
#' \describe{
#'   \item{`ind_df`}{Data frame, one row per individual: `id`, `src_x`,
#'     `src_y`, `birth_year`, `age`, `obs_year`, `age_c` (= age - 1),
#'     `src_cov`, and the recorded destination `dest` (row index of
#'     `trap_df`), `dest_x`, `dest_y`.}
#'   \item{`trap_df`}{Data frame, one row per trap: `trap_id`, `x`, `y`, `cov`.}
#'   \item{`op_mat`}{Integer matrix, traps by year, 1 if operational.}
#'   \item{`landscape`}{`SpatRaster` of the landscape covariate.}
#'   \item{`design`}{List of `landscape`, `traps`, `animals` for reuse.}
#' }
#'
#' @seealso [dsf()] to fit the model to the returned data.
#'
#' @examples
#' \donttest{
#' sim <- simulate_dsf(~ age_c, ~ 0 + cov,
#'                     beta_sigma = c(log(5), 0.3), beta_theta = 0.8,
#'                     exclude_natal = TRUE, seed = 1)
#' str(sim$ind_df)
#' fit <- dsf(~ age_c, ~ 0 + cov, sim$ind_df, sim$trap_df, sim$op_mat,
#'            exclude_natal = TRUE)
#' fit
#'
#' # reuse the same landscape and traps with a different parameter set
#' sim2 <- simulate_dsf(~ 1, ~ 1, beta_sigma = log(8), design = sim$design)
#' }
#'
#' @importFrom stats model.matrix runif sd
#' @importFrom fields matern.image.cov sim.rf
#' @importFrom terra rast ext crds extract
#' @export
simulate_dsf <- function(sigma_formula = ~ 1, theta_formula = ~ 1,
                         beta_sigma, beta_theta = NULL,
                         ntraps = 250, nyears = 5, nborn = 100,
                         prage = c(0.5, 0.3, 0.2),
                         L = 100, ncell = 1000, buffer = 10,
                         surface = TRUE, aRange = 10, smoothness = 0.5,
                         births = c("traps", "uniform"),
                         op_mat = NULL, exclude_natal = FALSE,
                         design = NULL, seed = NULL) {

  if (!is.null(seed)) set.seed(seed)
  births <- match.arg(births)

  ## ------------------------------------------------------------ design
  if (is.null(design)) {
    nall <- nyears * nborn

    if (surface) {
      ## aRange is in km; convert to cells. The circulant embedding needs a
      ## torus comfortably larger than the correlation range, so enlarge it
      ## (M, N) to at least 4 ranges, rounded up to a power of two.
      range_cells <- aRange / (L / ncell)
      M <- 2^ceiling(log2(max(2 * ncell, 4 * range_cells)))
      obj <- fields::matern.image.cov(grid = list(x = seq_len(ncell), y = seq_len(ncell)),
                                      aRange = range_cells, smoothness = smoothness,
                                      M = M, N = M, setup = TRUE)
      m <- fields::sim.rf(obj)
      m <- (m - mean(m)) / stats::sd(m)
    } else {
      m <- matrix(0, ncell, ncell)
    }
    landscape <- terra::rast(t(m[, ncell:1]), extent = terra::ext(0, L, 0, L), crs = "")
    names(landscape) <- "cov"

    cells  <- as.data.frame(terra::crds(landscape))
    inside <- cells$x > buffer & cells$x < L - buffer &
              cells$y > buffer & cells$y < L - buffer
    pick   <- sample(which(inside), ntraps)
    traps  <- data.frame(trap_id = seq_len(ntraps), x = cells$x[pick], y = cells$y[pick])
    traps$cov <- terra::extract(landscape, traps[, c("x", "y")])$cov

    if (births == "traps") {
      b  <- sample(ntraps, nall, replace = TRUE)
      sx <- traps$x[b]; sy <- traps$y[b]
    } else {
      sx <- stats::runif(nall, buffer, L - buffer)
      sy <- stats::runif(nall, buffer, L - buffer)
    }
    animals <- data.frame(id = seq_len(nall), src_x = sx, src_y = sy,
                          birth_year = rep(seq_len(nyears), each = nborn),
                          age = sample(1:3, nall, replace = TRUE, prob = prage))
    animals$obs_year <- animals$birth_year + animals$age
    animals$age_c    <- animals$age - 1
    animals$src_cov  <- terra::extract(landscape, animals[, c("src_x", "src_y")])$cov

    design <- list(landscape = landscape, traps = traps, animals = animals)
  } else {
    landscape <- design$landscape
    traps     <- design$traps
    animals   <- design$animals
  }
  n <- nrow(animals); J <- nrow(traps)

  ## ------------------------------------------------------------ sigma_i
  Xs <- stats::model.matrix(sigma_formula, animals)
  if (ncol(Xs) != length(beta_sigma))
    stop("beta_sigma needs ", ncol(Xs), " value(s): ",
         paste(colnames(Xs), collapse = ", "), call. = FALSE)
  log_sigma <- drop(Xs %*% beta_sigma)

  ## ------------------------------------------------------------ theta_j
  if (is_null_formula(theta_formula)) {
    if (!is.null(beta_theta))
      stop("beta_theta supplied but theta_formula is ~ 1", call. = FALSE)
    log_theta <- rep(0, J)
  } else {
    Xt <- stats::model.matrix(theta_formula, traps)
    if ("(Intercept)" %in% colnames(Xt))
      stop("theta_formula must have no intercept (not identifiable): use ~ 0 + ...",
           call. = FALSE)
    if (ncol(Xt) != length(beta_theta))
      stop("beta_theta needs ", ncol(Xt), " value(s): ",
           paste(colnames(Xt), collapse = ", "), call. = FALSE)
    log_theta <- drop(Xt %*% beta_theta)
  }
  if (!is.null(traps$effort)) {
    if (any(traps$effort <= 0)) stop("effort must be positive", call. = FALSE)
    log_theta <- log_theta + log(traps$effort)
  }

  ## ------------------------------------------------------------ op_mat
  if (is.null(op_mat)) {
    yrs    <- sort(unique(animals$obs_year))
    op_mat <- matrix(1L, J, length(yrs), dimnames = list(NULL, as.character(yrs)))
  }
  yr <- as.character(animals$obs_year)
  if (nrow(op_mat) != J)
    stop("op_mat must have one row per trap", call. = FALSE)
  if (!all(yr %in% colnames(op_mat)))
    stop("obs_year values missing from colnames(op_mat): ",
         paste(setdiff(unique(yr), colnames(op_mat)), collapse = ", "), call. = FALSE)

  ## ------------------------------------------------------------ hazard & draw
  d2 <- outer(animals$src_x, traps$x, "-")^2 + outer(animals$src_y, traps$y, "-")^2
  avail <- t(op_mat[, yr, drop = FALSE]) > 0
  if (exclude_natal) avail[cbind(seq_len(n), apply(d2, 1, which.min))] <- FALSE
  if (any(rowSums(avail) == 0))
    stop("individuals with no available destination: ",
         paste(which(rowSums(avail) == 0), collapse = ", "), call. = FALSE)

  log_lambda <- -d2 / (2 * exp(2 * log_sigma)) +
                matrix(log_theta, n, J, byrow = TRUE) +
                ifelse(avail, 0, -Inf)
  mx   <- apply(log_lambda, 1, max)
  pi   <- exp(log_lambda - mx) / rowSums(exp(log_lambda - mx))
  dest <- apply(pi, 1, function(p) sample.int(J, 1, prob = p))

  ## ------------------------------------------------------------ recorded data
  ind_df <- animals
  ind_df$dest   <- dest
  ind_df$dest_x <- traps$x[dest]
  ind_df$dest_y <- traps$y[dest]

  list(ind_df = ind_df, trap_df = traps, op_mat = op_mat,
       landscape = landscape, design = design)
}
