#' Assemble redistribution data for a destination selection function
#'
#' Combines the individual table, the destination table, the operational
#' matrix and any time-varying destination covariates into a single object
#' that [dsf()] fits and [simulate_dsf()] draws from. All structural checks
#' (destinations indexed correctly, observed destination operational in the
#' capture year, no empty choice sets) are made here.
#'
#' @details
#' The destination hazard may vary by year, \eqn{\log \theta_{jt} = Z_{jt}
#' \beta_\theta}, so that detection effort, bait status, habitat and so on
#' can change between years. Every individual belongs to one year (its
#' `obs_year`), so a term that is constant across destinations within a year
#' -- an intercept, or a year main effect -- cancels in the normalisation and
#' is not identifiable; [dsf()] refuses such terms. Time-invariant columns of
#' `traps` are recycled across years; `trap_tv` supplies columns that differ
#' between years.
#'
#' A column named `effort` (in `traps` or in `trap_tv`) is used by [dsf()] as
#' an offset \eqn{\log(\mathrm{effort}_{jt})}; it must be positive wherever
#' the destination is operational. Covariate values at destinations that are
#' not operational in a year are never used; `NA`s there are replaced by 0
#' so that the design matrix can be built.
#'
#' @param ind Data frame, one row per individual. Required: `src_x`, `src_y`
#'   (origin); `obs_year` (year of observation) unless there is a single
#'   year; `dest` (row index of the destination in `traps`) for fitting --
#'   it may be absent when the frame is used only for simulation. Any
#'   individual covariates for the `sigma` formula.
#' @param traps Data frame, one row per candidate destination. Required: `x`,
#'   `y`. Optional `trap_id`; any time-invariant covariates for the `theta`
#'   formula.
#' @param op Optional matrix of destinations by year (rows in `traps` order,
#'   columns named by year), positive where the destination was operational.
#'   `NULL` makes every destination available in every observed year.
#' @param trap_tv Optional time-varying destination covariates: a 3-d array
#'   `J x T x p` with `dimnames[[2]]` the years and `dimnames[[3]]` the
#'   covariate names (year dimnames may be omitted if the columns are in the
#'   same order as `op`), or a data frame in long format with columns
#'   `trap`, `year` and one column per covariate.
#' @param exclude_natal Logical; remove each individual's nearest destination
#'   from its choice set (use when origins are themselves destinations).
#'
#' @return An object of class `"dsf_frame"`: a list with `ind`, `traps`,
#'   `op`, `years`, `t` (integer year index per individual), `d2` (n-by-J
#'   squared Euclidean distances), `avail` (n-by-J logical), `obs` (n-by-2
#'   index of observed destinations, or `NULL`), `Z` (long data frame of
#'   J*T rows, destination within year, holding all destination covariates),
#'   `offset` (J-by-T matrix of `log(effort)`, or zero), `n`, `J`, `T`.
#'
#' @seealso [dsf()], [simulate_dsf()]
#' @export
dsf_frame <- function(ind, traps, op = NULL, trap_tv = NULL, exclude_natal = FALSE) {

  ## ------------------------------------------------------------ basic checks
  if (!all(c("src_x", "src_y") %in% names(ind)))
    stop("ind must contain src_x and src_y", call. = FALSE)
  if (!all(c("x", "y") %in% names(traps)))
    stop("traps must contain x and y", call. = FALSE)
  ind <- as.data.frame(ind); traps <- as.data.frame(traps)
  n <- nrow(ind); J <- nrow(traps)
  if (is.null(traps$trap_id)) traps$trap_id <- seq_len(J)
  has_dest <- !is.null(ind$dest)
  if (has_dest && !all(ind$dest %in% seq_len(J)))
    stop("ind$dest must index rows of traps", call. = FALSE)

  ## ------------------------------------------------------------ years
  if (is.null(op)) {
    if (is.null(ind$obs_year)) ind$obs_year <- 1L
    years <- sort(unique(ind$obs_year))
    op <- matrix(1L, J, length(years), dimnames = list(NULL, as.character(years)))
  } else {
    if (is.null(ind$obs_year))
      stop("ind$obs_year is required when op is supplied", call. = FALSE)
    if (nrow(op) != J) stop("op must have one row per destination", call. = FALSE)
    if (is.null(colnames(op))) stop("op must have column names giving the years", call. = FALSE)
    years <- colnames(op)
    miss <- setdiff(unique(as.character(ind$obs_year)), years)
    if (length(miss))
      stop("obs_year values missing from colnames(op): ", paste(miss, collapse = ", "),
           call. = FALSE)
  }
  years <- as.character(years); T <- length(years)
  t_i   <- match(as.character(ind$obs_year), years)

  ## ------------------------------------------------------------ distance & availability
  d2    <- outer(ind$src_x, traps$x, "-")^2 + outer(ind$src_y, traps$y, "-")^2
  avail <- t(op[, t_i, drop = FALSE]) > 0
  if (exclude_natal) avail[cbind(seq_len(n), apply(d2, 1, which.min))] <- FALSE
  if (any(rowSums(avail) == 0))
    stop("individuals with no available destination: ",
         paste(which(rowSums(avail) == 0), collapse = ", "), call. = FALSE)
  obs <- NULL
  if (has_dest) {
    obs <- cbind(seq_len(n), ind$dest)
    if (!all(avail[obs]))
      stop("individuals observed at a destination not operational in their obs_year: ",
           paste(which(!avail[obs]), collapse = ", "), call. = FALSE)
  }

  ## ------------------------------------------------------------ long destination table
  ## rows ordered destination within year, so column t of a J x T matrix is
  ## rows (t - 1) * J + 1:J
  Z <- traps[rep(seq_len(J), times = T), , drop = FALSE]
  Z$year <- rep(years, each = J)
  rownames(Z) <- NULL

  if (!is.null(trap_tv)) {
    if (is.array(trap_tv) && length(dim(trap_tv)) == 3) {
      dn <- dimnames(trap_tv)
      if (dim(trap_tv)[1] != J || dim(trap_tv)[2] != T)
        stop("trap_tv must be J x T x p with J = ", J, " and T = ", T, call. = FALSE)
      if (is.null(dn[[3]])) stop("trap_tv needs covariate names in dimnames[[3]]", call. = FALSE)
      if (!is.null(dn[[2]]) && !identical(as.character(dn[[2]]), years))
        trap_tv <- trap_tv[, years, , drop = FALSE]
      for (v in dn[[3]]) Z[[v]] <- as.vector(trap_tv[, , v])
    } else if (is.data.frame(trap_tv)) {
      if (!all(c("trap", "year") %in% names(trap_tv)))
        stop("long trap_tv needs columns trap and year", call. = FALSE)
      key_Z  <- paste(Z$trap_id, Z$year)
      key_tv <- paste(trap_tv$trap, trap_tv$year)
      if (anyDuplicated(key_tv)) stop("trap_tv has duplicated trap-year rows", call. = FALSE)
      m <- match(key_Z, key_tv)
      for (v in setdiff(names(trap_tv), c("trap", "year"))) Z[[v]] <- trap_tv[[v]][m]
    } else stop("trap_tv must be a J x T x p array or a long data frame", call. = FALSE)
  }

  ## covariates at non-operational destination-years are unused: fill NAs
  open <- as.vector(op[, years] > 0)
  for (v in setdiff(names(Z), c("trap_id", "x", "y", "year"))) {
    bad <- is.na(Z[[v]])
    if (any(bad & open))
      stop("NA in destination covariate '", v, "' at operational destination-years",
           call. = FALSE)
    if (any(bad)) Z[[v]][bad] <- 0
  }

  ## ------------------------------------------------------------ effort offset
  offset <- matrix(0, J, T, dimnames = list(NULL, years))
  if (!is.null(Z$effort)) {
    e <- Z$effort
    if (any(e[open] <= 0)) stop("effort must be positive where operational", call. = FALSE)
    e[!open] <- 1
    offset <- matrix(log(e), J, T, dimnames = list(NULL, years))
  }

  structure(list(ind = ind, traps = traps, op = op, years = years, t = t_i,
                 d2 = d2, avail = avail, obs = obs, Z = Z, offset = offset,
                 exclude_natal = exclude_natal, n = n, J = J, T = T),
            class = "dsf_frame")
}

#' @describeIn dsf_frame Print a summary of the frame.
#' @param x A `"dsf_frame"` object.
#' @param ... Ignored.
#' @export
print.dsf_frame <- function(x, ...) {
  cat("dsf_frame:", x$n, "individuals,", x$J, "destinations,", x$T,
      if (x$T == 1) "year\n" else "years\n")
  cat("  years:", paste(x$years, collapse = ", "), "\n")
  cat("  destinations per individual (mean):", round(mean(rowSums(x$avail)), 1), "\n")
  tv <- setdiff(names(x$Z), c(names(x$traps), "year"))
  cat("  destination covariates:",
      paste(setdiff(names(x$traps), c("trap_id", "x", "y")), collapse = ", "),
      if (length(tv)) paste0("  time-varying: ", paste(tv, collapse = ", ")) else "", "\n")
  invisible(x)
}

## Build the theta design matrix on the long table and check identifiability:
## no intercept, and no column constant across destinations within every year.
theta_design <- function(theta, frame) {
  if (is_null_formula(theta)) return(matrix(0, frame$J * frame$T, 0))
  Xt <- stats::model.matrix(theta, frame$Z)
  if ("(Intercept)" %in% colnames(Xt))
    stop("theta must have no intercept (not identifiable): use ~ 0 + ...", call. = FALSE)
  yr <- rep(seq_len(frame$T), each = frame$J)
  const <- vapply(seq_len(ncol(Xt)), function(k)
    all(tapply(Xt[, k], yr, function(z) diff(range(z)) == 0)), logical(1))
  if (any(const))
    stop("theta term(s) constant across destinations within year, not identifiable: ",
         paste(colnames(Xt)[const], collapse = ", "), call. = FALSE)
  Xt
}

## n x J matrix of log theta_{j, t_i} from the long linear predictor
theta_rows <- function(lp_long, frame) {
  lt <- matrix(lp_long, frame$J, frame$T)
  t(lt[, frame$t, drop = FALSE])
}
