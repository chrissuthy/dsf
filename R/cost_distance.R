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
#' should be coarse enough that Dijkstra from every origin at every
#' likelihood evaluation is affordable (use [terra::aggregate()]
#' beforehand); and missing cells are either an error (default) or, with
#' `na_barrier = TRUE`, impassable: they are removed from the graph, so
#' paths go around them under every model (including the flat surface) and
#' an origin and destination that are not connected through non-missing
#' cells raise an error. Use this for sea, lakes or any region movement
#' cannot cross; it is not a substitute for extending a covariate beyond the
#' study boundary.
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
#' @param na_barrier Logical; treat missing cells as impassable rather than
#'   as an error? Default `FALSE`. A cell is missing if any cost covariate
#'   is `NA` there (for `cost = ~ 1`, the first layer of `landscape` is
#'   used).
#' @param window Optional focal-window radius in map units. For each origin,
#'   only destinations within `window` (Euclidean) receive a finite distance;
#'   the rest are `Inf` (zero hazard), and the least-cost calculation runs on
#'   the subgraph of a box of half-width `1.25 * window` around the origin,
#'   built once. This cuts the per-evaluation cost from (origins × whole
#'   raster) to (origins × box). Choose `window` as several dispersal scales
#'   (4–5 σ) so that the dropped tail mass is negligible; [dsf()] refuses a
#'   window that excludes any observed destination. `NULL` (default)
#'   computes every origin to every destination on the whole raster.
#' @param standardise Logical; centre and scale each cost covariate (mean 0,
#'   unit SD over the raster cells) before building the cost surface?
#'   Default `TRUE`. This keeps \eqn{\exp(\alpha^\top c)} well behaved during
#'   optimisation and makes \eqn{\alpha} an effect per SD of the covariate;
#'   the means and SDs used are returned in `scale` so that \eqn{\alpha} can
#'   be converted to raw units (`alpha / sd`). A message reports the scaling.
#'
#' @return A list with `fun(alpha)`, returning the n-by-J matrix of squared
#'   least-cost distances for resistance vector `alpha` (`Inf` beyond the
#'   window); `names`, the cost covariate names; `raster`, the `SpatRaster`
#'   used; `scale`, a data frame of the mean and SD applied to each covariate
#'   (all 0 and 1 if `standardise = FALSE`); `window`; and `in_window`, the
#'   n-by-J logical matrix of pairs with finite distance.
#'
#' @importFrom terra res cellFromXY values ext ncell xyFromCell rowColFromCell
#' @importFrom igraph make_graph distances
#' @importFrom stats sd
#' @keywords internal
make_cost_fun <- function(landscape, cost, src, dst, standardise = TRUE,
                          na_barrier = FALSE, window = NULL) {
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
  ## missing cells: the layers named in cost, or the first layer for ~ 1
  miss_cell <- if (length(vars)) {
    Reduce(`|`, lapply(vars, function(v) !is.finite(cell_df[[v]])))
  } else !is.finite(terra::values(r[[1]])[, 1])
  if (any(miss_cell)) {
    if (!na_barrier)
      stop(sum(miss_cell), " of ", ncells, " raster cells have missing values. ",
           "The resistance surface must be complete; fill or mask missing cells ",
           "before fitting, or set na_barrier = TRUE to treat them as impassable.",
           call. = FALSE)
    if (any(miss_cell[c(src_cell, dst_cell)]))
      stop(sum(miss_cell[src_cell]), " origin(s) and ", sum(miss_cell[dst_cell]),
           " destination(s) fall on missing (impassable) cells.", call. = FALSE)
    for (v in vars) cell_df[[v]][miss_cell] <- 0     # never used: edges removed below
  }

  ## ------------------------------------------------------------ standardise
  scale_df <- data.frame(var = vars, mean = rep(0, length(vars)), sd = rep(1, length(vars)))
  if (standardise && length(vars)) {
    for (i in seq_along(vars)) {
      v <- vars[i]
      m <- mean(cell_df[[v]][!miss_cell]); s <- stats::sd(cell_df[[v]][!miss_cell])
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
  ## impassable cells: drop every edge that touches one
  keep <- !(miss_cell[from] | miss_cell[to])
  from <- from[keep]; to <- to[keep]; len <- len[keep]
  usrc <- sort(unique(src_cell)); spos <- match(src_cell, usrc)
  udst <- sort(unique(dst_cell)); dpos <- match(dst_cell, udst)
  n_src <- nrow(src); n_dst <- nrow(dst)

  ## ------------------------------------------------------------ Dijkstra driver
  if (is.null(window)) {
    ## one graph, every origin to every destination
    g <- igraph::make_graph(edges = rbind(from, to), n = nr * nc, directed = FALSE)
    dist_fun <- function(w) igraph::distances(g, v = usrc, to = udst, weights = w)[spos, dpos, drop = FALSE]
    in_window <- matrix(TRUE, n_src, n_dst)
  } else {
    ## focal windows: for each unique origin cell, destinations within
    ## `window` (Euclidean) and the induced subgraph of a box of half-width
    ## 1.25 * window around the origin, built once; distances to
    ## destinations beyond the window are Inf (zero hazard)
    if (!is.numeric(window) || length(window) != 1 || window <= 0)
      stop("window must be a single positive distance in map units", call. = FALSE)
    half <- 1.25 * window
    e_row <- (from - 1) %/% nc + 1; e_col <- (from - 1) %% nc + 1
    t_row <- (to   - 1) %/% nc + 1; t_col <- (to   - 1) %% nc + 1
    xy_src <- terra::xyFromCell(r, usrc)
    in_window <- matrix(FALSE, n_src, n_dst)
    focal <- vector("list", length(usrc))
    for (u in seq_along(usrc)) {
      rows_u <- which(spos == u)
      d_euc  <- sqrt((dst[, 1] - xy_src[u, 1])^2 + (dst[, 2] - xy_src[u, 2])^2)
      dcols  <- which(d_euc <= window)
      in_window[rows_u, dcols] <- TRUE
      rc <- terra::rowColFromCell(r, usrc[u])
      hr <- ceiling(half / cs)
      r0 <- max(1, rc[1] - hr); r1 <- min(nr, rc[1] + hr)
      c0 <- max(1, rc[2] - hr); c1 <- min(nc, rc[2] + hr)
      eidx <- which(e_row >= r0 & e_row <= r1 & e_col >= c0 & e_col <= c1 &
                    t_row >= r0 & t_row <= r1 & t_col >= c0 & t_col <= c1)
      cells_u <- as.vector(id[r0:r1, c0:c1])
      loc <- match(c(from[eidx], to[eidx]), cells_u)
      g_u <- igraph::make_graph(edges = rbind(loc[seq_along(eidx)], loc[-seq_along(eidx)]),
                                n = length(cells_u), directed = FALSE)
      ucell <- sort(unique(dst_cell[dcols]))
      focal[[u]] <- list(rows = rows_u, dcols = dcols, eidx = eidx, g = g_u,
                         v = match(usrc[u], cells_u),
                         to = match(ucell, cells_u),
                         dpos = match(dst_cell[dcols], ucell))
    }
    dist_fun <- function(w) {
      D <- matrix(Inf, n_src, n_dst)
      for (f in focal) {
        if (!length(f$dcols)) next
        d <- igraph::distances(f$g, v = f$v, to = f$to, weights = w[f$eidx])[1, f$dpos]
        D[f$rows, f$dcols] <- matrix(d, length(f$rows), length(f$dcols), byrow = TRUE)
      }
      D
    }
  }

  fun <- function(alpha = numeric(0)) {
    if (length(alpha) != ncol(Xc))
      stop("alpha must have length ", ncol(Xc), call. = FALSE)
    if (any(!is.finite(alpha))) stop("non-finite resistance parameter", call. = FALSE)
    ## cap the linear predictor so exp() cannot overflow to Inf (which would
    ## give NaN edge weights); costs above exp(50) are effectively barriers
    lp <- if (ncol(Xc)) as.vector(Xc %*% alpha) else rep(0, ncells)   # as.vector drops row names
    cell_cost <- exp(pmin(lp, 50))
    w <- 0.5 * (cell_cost[from] + cell_cost[to]) * len
    dist_fun(w)^2
  }

  if (any(miss_cell)) {
    ## every origin must reach every in-window destination through passable cells
    D0 <- fun(rep(0, ncol(Xc)))
    bad <- in_window & !is.finite(D0)
    if (any(bad))
      stop(sum(bad), " origin-destination pair(s) are not connected through ",
           "non-missing cells; with na_barrier = TRUE the passable region must ",
           "connect every origin to every destination", 
           if (!is.null(window)) " within the window", ".", call. = FALSE)
  }
  list(fun = fun, names = colnames(Xc), raster = r, scale = scale_df,
       window = window, in_window = in_window)
}
