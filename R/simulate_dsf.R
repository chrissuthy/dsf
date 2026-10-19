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
#' To hold the design fixed across replicates while redrawing dispersal, give
#' the same `landscape_seed` on every call (or pass back `$design`). To let the
#' design vary as well, leave `landscape_seed = NULL`.
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
#' @param surface Logical; simulate a landscape covariate (a Gaussian random
#'   field)? Only needed for destination covariates or the non-Euclidean DSF.
#'   If `FALSE` (default) the covariate is zero everywhere and the raster is
#'   constant.
#' @param aRange Range parameter of the Matérn covariance, in km (the same
#'   units as `L`), so the patch size does not change with `ncell`.
#' @param smoothness Matérn smoothness; 0.5 is exponential.
#' @param trap_layout `"random"` (default) places traps uniformly at random
#'   inside the buffer; `"grid"` places them on a regular lattice, so that a
#'   given `ntraps` is a density rather than a draw. The lattice is `nx` by
#'   `ny` with `nx * ny = ntraps`, as close to square as the factorisation of
#'   `ntraps` allows (a prime `ntraps` gives a single row; choose a composite
#'   number, and a square number for equal spacing in both directions).
#' @param jitter For `trap_layout = "grid"`, the maximum random displacement
#'   of each trap from its lattice node, as a fraction of the lattice spacing
#'   (uniform in each direction). `0` (default) gives a regular lattice.
#' @param trap_cov Destination covariates. `NULL` (default) extracts `cov`
#'   at each trap from the landscape raster. Otherwise a numeric vector (one
#'   value per trap, stored as `cov`), a data frame with one row per trap
#'   (columns added to `trap_df`, replacing any of the same name), or a
#'   function of `ntraps` returning either. A function is evaluated inside
#'   the `landscape_seed` scope, so the covariate is part of the fixed design.
#' @param births `"traps"` places each origin at a randomly chosen trap
#'   (the parent's capture location); `"uniform"` draws origins uniformly
#'   inside the buffer.
#' @param cost Optional one-sided formula naming landscape layers that define
#'   a resistance surface (e.g. `~ cov`); requires `surface = TRUE`. With it,
#'   dispersal uses least-cost distance with cell cost
#'   \eqn{\exp(\alpha^\top c(x))}; see [make_cost_fun()].
#' @param beta_cost Resistance coefficients \eqn{\alpha} for `cost`.
#' @param cost_res Cell size for least-cost paths; the landscape is aggregated
#'   to this resolution.
#' @param op_mat Optional trap-by-year operational matrix (0/1), with one row
#'   per trap and columns named by year. If `NULL`, every trap is open in
#'   every observation year.
#' @param exclude_natal Logical; remove each individual's nearest trap
#'   (distance zero when `births = "traps"`) from its choice set?
#' @param landscape_seed Optional integer. If supplied, the landscape, trap
#'   array and cohort design are drawn from a generator seeded with this value
#'   and the previous RNG state is then restored, so repeated calls with the
#'   same `landscape_seed` share one design while the dispersal draw varies.
#'   If `NULL` (default), each call draws a fresh design.
#' @param design Optional design list (`landscape`, `traps`, `animals`) from a
#'   previous call, returned as `$design`; overrides `landscape_seed`.
#' @param seed Optional integer passed to [set.seed()] for the whole call.
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
#' # destination covariate needs a landscape surface
#' sim <- simulate_dsf(~ age_c, ~ 0 + cov,
#'                     beta_sigma = c(log(5), 0.3), beta_theta = 0.8,
#'                     surface = TRUE, exclude_natal = TRUE, seed = 1)
#' str(sim$ind_df)
#' fit <- dsf(~ age_c, ~ 0 + cov, sim$ind_df, sim$trap_df, sim$op_mat,
#'            exclude_natal = TRUE)
#' fit
#'
#' # same design, different dispersal draw: fix the landscape seed
#' sim2 <- simulate_dsf(~ 1, ~ 1, beta_sigma = log(5), landscape_seed = 99)
#' sim3 <- simulate_dsf(~ 1, ~ 1, beta_sigma = log(5), landscape_seed = 99)
#' identical(sim2$trap_df, sim3$trap_df)      # TRUE
#' identical(sim2$ind_df$dest, sim3$ind_df$dest)   # FALSE
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
                         surface = FALSE, aRange = 10, smoothness = 0.5,
                         trap_layout = c("random", "grid"), jitter = 0,
                         trap_cov = NULL, births = c("traps", "uniform"),
                         cost = NULL, beta_cost = NULL, cost_res = 1,
                         op_mat = NULL, exclude_natal = FALSE,
                         landscape_seed = NULL, design = NULL, seed = NULL) {

  if (!is.null(seed)) set.seed(seed)
  births <- match.arg(births)
  trap_layout <- match.arg(trap_layout)

  ## ------------------------------------------------------------ design
  if (is.null(design)) {
    design <- make_design(ntraps, nyears, nborn, prage, L, ncell, buffer,
                          surface, aRange, smoothness, births, trap_layout,
                          jitter, trap_cov, landscape_seed)
  }
  landscape <- design$landscape
  traps     <- design$traps
  animals   <- design$animals
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
  if (!is.null(cost)) {
    if (!surface && is.null(design))
      stop("cost requires a landscape surface: set surface = TRUE", call. = FALSE)
    cf <- make_cost_fun(landscape, cost, cost_res,
                        animals[, c("src_x", "src_y")], traps[, c("x", "y")])
    if (length(beta_cost) != length(cf$names))
      stop("beta_cost needs ", length(cf$names), " value(s): ",
           paste(cf$names, collapse = ", "), call. = FALSE)
    d2 <- cf$fun(beta_cost)                    # squared least-cost distance
  }
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

## Build landscape, traps and cohorts. With a seed, the RNG is seeded inside a
## local scope and restored afterwards, so the design is reproducible while the
## caller's random stream (and hence the dispersal draw) is untouched.
make_design <- function(ntraps, nyears, nborn, prage, L, ncell, buffer,
                        surface, aRange, smoothness, births, trap_layout,
                        jitter, trap_cov, landscape_seed) {
  if (!is.null(landscape_seed)) {
    had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
    if (had_seed) old_seed <- get(".Random.seed", envir = globalenv())
    on.exit({
      if (had_seed) assign(".Random.seed", old_seed, envir = globalenv())
      else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
        rm(".Random.seed", envir = globalenv())
    }, add = TRUE)
    set.seed(landscape_seed)
  }
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

  if (trap_layout == "random") {
    cells  <- as.data.frame(terra::crds(landscape))
    inside <- cells$x > buffer & cells$x < L - buffer &
              cells$y > buffer & cells$y < L - buffer
    pick   <- sample(which(inside), ntraps)
    tx <- cells$x[pick]; ty <- cells$y[pick]
  } else {
    ## rectangular lattice nx * ny = ntraps, as close to square as the
    ## factorisation allows, optionally jittered by a fraction of the spacing
    nx <- max(which(ntraps %% seq_len(floor(sqrt(ntraps))) == 0))
    ny <- ntraps / nx
    if (ny < nx) { tmp <- nx; nx <- ny; ny <- tmp }
    spx <- (L - 2 * buffer) / nx; spy <- (L - 2 * buffer) / ny
    gg  <- expand.grid(x = buffer + spx * (seq_len(nx) - 0.5),
                       y = buffer + spy * (seq_len(ny) - 0.5))
    tx <- gg$x; ty <- gg$y
    if (jitter > 0) {
      tx <- tx + stats::runif(ntraps, -jitter * spx, jitter * spx)
      ty <- ty + stats::runif(ntraps, -jitter * spy, jitter * spy)
    }
  }
  traps <- data.frame(trap_id = seq_len(ntraps), x = tx, y = ty)
  traps$cov <- terra::extract(landscape, traps[, c("x", "y")])$cov
  if (!is.null(trap_cov)) {
    if (is.function(trap_cov)) trap_cov <- trap_cov(ntraps)
    if (is.numeric(trap_cov) && is.null(dim(trap_cov))) trap_cov <- data.frame(cov = trap_cov)
    trap_cov <- as.data.frame(trap_cov)
    if (nrow(trap_cov) != ntraps)
      stop("trap_cov must have one row (or value) per trap", call. = FALSE)
    if (any(c("trap_id", "x", "y") %in% names(trap_cov)))
      stop("trap_cov columns cannot be named trap_id, x or y", call. = FALSE)
    traps[names(trap_cov)] <- trap_cov
  }

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

  list(landscape = landscape, traps = traps, animals = animals)
}
