#' Build a least-cost distance calculator from a landscape raster
#'
#' Builds the 8-neighbour graph of a resistance raster once and returns a
#' function that computes squared least-cost distances from a set of origins
#' to a set of destinations for any value of the resistance parameters. Used
#' internally by [dsf()] and [simulate_dsf()] for the non-Euclidean model.
#'
#' @details
#' The cost of cell \eqn{x} is \eqn{\exp(\alpha^\top c(x))} where \eqn{c(x)}
#' is the vector of raster values named in `cost` at that cell. The weight of
#' the edge between adjacent cells \eqn{x, x'} is the mean of their costs
#' times the Euclidean step length, and the ecological distance is the
#' least-cost path length over the graph ([igraph::distances()]). With
#' \eqn{\alpha = 0} every cell has cost 1 and the result approximates
#' Euclidean distance on the grid.
#'
#' The result is identical (to floating-point precision) to
#' `gdistance::costDistance()` on a transition layer built with conductance
#' `1 / mean(cost)` over 8 directions and `geoCorrection(type = "c")`; the
#' difference is that the adjacency is built once here and only the edge
#' weights are recomputed when \eqn{\alpha} changes, which matters because
#' the surface changes at every likelihood evaluation.
#'
#' The raster is used at the resolution supplied. It is the user's
#' responsibility to prepare it: it must contain every origin and
#' destination, ideally with a buffer of about two dispersal scales beyond
#' them so that least-cost paths are not forced along the raster edge; it
#' must have no missing cells; and it should be coarse enough that Dijkstra
#' from every origin at every likelihood evaluation is affordable (use
#' [terra::aggregate()] beforehand). These conditions are checked and an
#' informative error is raised if they fail.
#'
#' @param landscape A `SpatRaster` with one layer per cost covariate.
#' @param cost One-sided formula naming the layers, e.g. `~ cov`. No intercept
#'   is used: a constant multiplier on all costs rescales every distance by
#'   the same factor and is absorbed by \eqn{\sigma}. `~ 1` gives a flat
#'   surface (every cell cost 1, no resistance parameters): grid-path
#'   distance on the same lattice, the like-for-like null for a resistance
#'   model.
#' @param src Two-column matrix or data frame of origin coordinates.
#' @param dst Two-column matrix or data frame of destination coordinates.
#' @param standardise Logical; centre and scale each cost covariate (mean 0,
#'   unit SD over the raster cells) before building the cost surface?
#'   Default `TRUE`. This keeps \eqn{\exp(\alpha^\top c)} well behaved during
#'   optimisation and makes \eqn{\alpha} an effect per SD of the covariate;
#'   the means and SDs used are returned in `scale` so that \eqn{\alpha} can
#'   be converted to raw units (`alpha / sd`). A message reports the scaling.
#'
#' @return A list with `fun(alpha)`, returning the n-by-J matrix of squared
#'   least-cost distances for resistance vector `alpha`; `names`, the cost
#'   covariate names; `raster`, the `SpatRaster` used; and `scale`, a data
#'   frame of the mean and SD applied to each covariate (all 0 and 1 if
#'   `standardise = FALSE`).
#'
#' @importFrom terra res cellFromXY values ext ncell
#' @importFrom igraph make_graph distances
#' @importFrom stats sd
#' @keywords internal
make_cost_fun <- function(landscape, cost, src, dst, standardise = TRUE) {
  if (!requireNamespace("igraph", quietly = TRUE))
    stop("package 'igraph' is required for the non-Euclidean DSF", call. = FALSE)
  r <- landscape

  ## ------------------------------------------------------------ checks
  vars <- all.vars(cost)
  miss <- setdiff(vars, names(r))
  if (length(miss)) stop("cost covariates not in landscape: ", paste(miss, collapse = ", "),
                         call. = FALSE)

  src <- as.matrix(src); dst <- as.matrix(dst)
  src_cell <- terra::cellFromXY(r, src)
  dst_cell <- terra::cellFromXY(r, dst)
  n_out <- sum(is.na(src_cell)) + sum(is.na(dst_cell))
  if (n_out > 0) {
    e  <- as.vector(terra::ext(r))
    bb <- apply(rbind(src, dst), 2, range)
    stop(n_out, " origin(s)/destination(s) fall outside the landscape raster. ",
         "The resistance surface must cover every origin and destination, ideally ",
         "with a buffer of about two dispersal scales beyond them so that ",
         "least-cost paths are not forced along the raster edge. Raster extent: x ",
         signif(e[1], 5), " to ", signif(e[2], 5), ", y ", signif(e[3], 5), " to ",
         signif(e[4], 5), "; points span x ", signif(bb[1, 1], 5), " to ",
         signif(bb[2, 1], 5), ", y ", signif(bb[1, 2], 5), " to ", signif(bb[2, 2], 5),
         ".", call. = FALSE)
  }

  ncells  <- terra::ncell(r)
  cell_df <- if (length(vars)) as.data.frame(terra::values(r[[vars]])) else
    data.frame(row.names = seq_len(ncells))
  for (v in vars) {
    nbad <- sum(!is.finite(cell_df[[v]]))
    if (nbad > 0)
      stop("cost covariate '", v, "' has ", nbad, " of ", nrow(cell_df),
           " cells with missing values. The resistance surface must be complete; ",
           "fill or mask missing cells before fitting.", call. = FALSE)
  }

  ## ------------------------------------------------------------ standardise
  scale_df <- data.frame(var = vars, mean = rep(0, length(vars)), sd = rep(1, length(vars)))
  if (standardise && length(vars)) {
    for (i in seq_along(vars)) {
      v <- vars[i]; m <- mean(cell_df[[v]]); s <- stats::sd(cell_df[[v]])
      if (!is.finite(s) || s == 0) stop("cost covariate '", v, "' is constant", call. = FALSE)
      cell_df[[v]] <- (cell_df[[v]] - m) / s
      scale_df$mean[i] <- m; scale_df$sd[i] <- s
    }
    message("cost covariate(s) standardised over raster cells: ",
            paste0(vars, " (mean ", signif(scale_df$mean, 3), ", sd ",
                   signif(scale_df$sd, 3), ")", collapse = "; "),
            "; cost coefficients are per SD")
  }
  Xc <- if (length(vars)) stats::model.matrix(stats::update(cost, ~ . + 0), cell_df) else
    matrix(0, ncells, 0)                       # flat surface: cost = ~ 1

  ## ------------------------------------------------------------ graph
  ## 8-neighbour edge list on the grid (cell ids are row-major in terra)
  nr <- nrow(r); nc <- ncol(r); cs <- terra::res(r)[1]
  id <- matrix(seq_len(nr * nc), nr, nc, byrow = TRUE)     # id[row, col]
  from <- to <- integer(0); len <- numeric(0)
  for (d in list(c(0, 1), c(1, 0), c(1, 1), c(1, -1))) {
    dr <- d[1]; dc <- d[2]
    rows <- seq_len(nr - abs(dr)); cols <- seq_len(nc - abs(dc))
    a <- id[rows + max(0, -dr), cols + max(0, -dc)]
    b <- id[rows + max(0,  dr), cols + max(0,  dc)]
    from <- c(from, a); to <- c(to, b)
    len  <- c(len, rep(cs * sqrt(dr^2 + dc^2), length(a)))
  }
  g <- igraph::make_graph(edges = rbind(from, to), n = nr * nc, directed = FALSE)

  usrc <- sort(unique(src_cell)); spos <- match(src_cell, usrc)
  udst <- sort(unique(dst_cell)); dpos <- match(dst_cell, udst)

  fun <- function(alpha = numeric(0)) {
    if (length(alpha) != ncol(Xc))
      stop("alpha must have length ", ncol(Xc), call. = FALSE)
    if (any(!is.finite(alpha))) stop("non-finite resistance parameter", call. = FALSE)
    ## cap the linear predictor so exp() cannot overflow to Inf (which would
    ## give NaN edge weights); costs above exp(50) are effectively barriers
    lp <- if (ncol(Xc)) drop(Xc %*% alpha) else rep(0, ncells)
    cell_cost <- exp(pmin(lp, 50))
    w <- 0.5 * (cell_cost[from] + cell_cost[to]) * len
    D <- igraph::distances(g, v = usrc, to = udst, weights = w)
    D[spos, dpos, drop = FALSE]^2
  }
  list(fun = fun, names = colnames(Xc), raster = r, scale = scale_df)
}
