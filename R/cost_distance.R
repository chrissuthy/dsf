#' Build a least-cost distance calculator from a landscape raster
#'
#' Aggregates a landscape raster to the resolution used for path
#' calculations, builds the 8-neighbour graph once, and returns a function
#' that computes squared least-cost distances from a set of origins to a set
#' of destinations for any value of the resistance parameters. Used
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
#' @param landscape A `SpatRaster` with one layer per cost covariate.
#' @param cost One-sided formula naming the layers, e.g. `~ cov`. No intercept
#'   is used: a constant multiplier on all costs rescales every distance by
#'   the same factor and is absorbed by \eqn{\sigma}.
#' @param cost_res Cell size (map units) for the path graph. The landscape is
#'   aggregated by block mean to this resolution before the graph is built.
#' @param src Two-column matrix or data frame of origin coordinates.
#' @param dst Two-column matrix or data frame of destination coordinates.
#' @param standardise Logical; centre and scale each cost covariate (mean 0,
#'   unit SD over the aggregated raster cells) before building the cost
#'   surface? Default `TRUE`. This keeps \eqn{\exp(\alpha^\top c)} well behaved
#'   during optimisation and makes \eqn{\alpha} an effect per SD of the
#'   covariate; the means and SDs used are returned in `scale` so that
#'   \eqn{\alpha} can be converted to raw units (`alpha / sd`). A message
#'   reports the scaling.
#'
#' @return A list with `fun(alpha)`, returning the n-by-J matrix of squared
#'   least-cost distances for resistance vector `alpha`; `names`, the cost
#'   covariate names; `raster`, the aggregated `SpatRaster`; and `scale`, a
#'   data frame of the mean and SD applied to each covariate (all 0 and 1 if
#'   `standardise = FALSE`).
#'
#' @importFrom terra aggregate res cellFromXY values
#' @importFrom igraph make_graph distances
#' @importFrom stats sd
#' @keywords internal
make_cost_fun <- function(landscape, cost, cost_res, src, dst, standardise = TRUE) {
  if (!requireNamespace("igraph", quietly = TRUE))
    stop("package 'igraph' is required for the non-Euclidean DSF", call. = FALSE)

  ## aggregate to the path resolution
  fac <- round(cost_res / terra::res(landscape)[1])
  r   <- if (fac > 1) terra::aggregate(landscape, fact = fac, fun = "mean") else landscape

  ## cost covariates: design matrix over cells, no intercept
  vars <- all.vars(cost)
  miss <- setdiff(vars, names(r))
  if (length(miss)) stop("cost covariates not in landscape: ", paste(miss, collapse = ", "),
                         call. = FALSE)
  cell_df <- as.data.frame(terra::values(r[[vars]]))
  if (any(!is.finite(as.matrix(cell_df))))
    stop("NA in cost covariates after aggregation", call. = FALSE)
  scale_df <- data.frame(var = vars, mean = 0, sd = 1)
  if (standardise) {
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
  Xc <- stats::model.matrix(stats::update(cost, ~ . + 0), cell_df)

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

  ## map origins and destinations to cells
  src_cell <- terra::cellFromXY(r, as.matrix(src))
  dst_cell <- terra::cellFromXY(r, as.matrix(dst))
  if (anyNA(src_cell) || anyNA(dst_cell))
    stop("origins or destinations fall outside the landscape", call. = FALSE)
  usrc <- sort(unique(src_cell)); spos <- match(src_cell, usrc)
  udst <- sort(unique(dst_cell)); dpos <- match(dst_cell, udst)

  fun <- function(alpha) {
    if (any(!is.finite(alpha))) stop("non-finite resistance parameter", call. = FALSE)
    ## cap the linear predictor so exp() cannot overflow to Inf (which would
    ## give NaN edge weights); costs above exp(50) are effectively barriers
    cell_cost <- exp(pmin(drop(Xc %*% alpha), 50))
    w <- 0.5 * (cell_cost[from] + cell_cost[to]) * len
    D <- igraph::distances(g, v = usrc, to = udst, weights = w)
    D[spos, dpos, drop = FALSE]^2
  }
  list(fun = fun, names = colnames(Xc), raster = r, scale = scale_df)
}
