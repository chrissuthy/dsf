test_that("dsf recovers simulation parameters", {
  skip_on_cran()
  sim <- simulate_dsf(~ age_c, ~ 0 + cov, beta_sigma = c(log(5), 0.3), beta_theta = 0.8,
                      ncell = 200, surface = TRUE, exclude_natal = TRUE, seed = 1)
  fit <- dsf(~ age_c, ~ 0 + cov, sim$ind_df, sim$trap_df, sim$op_mat, exclude_natal = TRUE)
  expect_equal(fit$convergence, 0)
  truth <- c(log(5), 0.3, 0.8)
  expect_true(all(abs(fit$coef - truth) < 3 * fit$se))
})

test_that("theta intercept is refused", {
  sim <- simulate_dsf(~ 1, ~ 1, beta_sigma = log(5), ncell = 100, surface = TRUE, seed = 2)
  expect_error(dsf(~ 1, ~ cov, sim$ind_df, sim$trap_df, sim$op_mat), "intercept")
})

test_that("observed destination must be available", {
  sim <- simulate_dsf(~ 1, ~ 1, beta_sigma = log(5), ncell = 100, seed = 3)
  op <- sim$op_mat; op[sim$ind_df$dest[1], as.character(sim$ind_df$obs_year[1])] <- 0L
  expect_error(dsf(~ 1, ~ 1, sim$ind_df, sim$trap_df, op), "not operational")
})

test_that("dsf_frame reproduces the direct fit and refuses year effects", {
  sim <- simulate_dsf(~ 1, ~ 0 + cov, beta_sigma = log(5), beta_theta = 0.5,
                      ncell = 100, trap_cov = function(n) runif(n, -2, 2), seed = 11)
  fr <- dsf_frame(sim$ind_df, sim$trap_df, sim$op_mat)
  expect_s3_class(fr, "dsf_frame")
  expect_equal(nrow(fr$Z), fr$J * fr$T)
  f1 <- dsf(~ 1, ~ 0 + cov, sim$ind_df, sim$trap_df, sim$op_mat)
  f2 <- dsf(~ 1, ~ 0 + cov, fr)
  f3 <- dsf(~ 1, ~ 0 + cov, sim$frame)          # frame returned by the simulator
  expect_equal(f1$coef, f2$coef)
  expect_equal(f1$logLik, f2$logLik)
  expect_equal(f1$coef, f3$coef)
  expect_false(is.null(sim$frame$obs))
  ## a year main effect is constant across destinations within year
  expect_error(dsf(~ 1, ~ 0 + year, fr), "not identifiable")
  expect_error(dsf(~ 1, ~ 0 + cov + year, fr), "not identifiable")
})

test_that("time-invariant trap_tv equals the time-invariant fit", {
  sim <- simulate_dsf(~ 1, ~ 0 + cov, beta_sigma = log(5), beta_theta = 0.5,
                      ncell = 100, trap_cov = function(n) runif(n, -2, 2), seed = 12)
  yrs <- colnames(sim$op_mat)
  tv  <- array(sim$trap_df$cov, dim = c(nrow(sim$trap_df), length(yrs), 1),
               dimnames = list(NULL, yrs, "cov_tv"))
  f1 <- dsf(~ 1, ~ 0 + cov,    sim$ind_df, sim$trap_df, sim$op_mat)
  f2 <- dsf(~ 1, ~ 0 + cov_tv, sim$ind_df, sim$trap_df, sim$op_mat, trap_tv = tv)
  expect_equal(unname(f1$coef), unname(f2$coef), tolerance = 1e-6)
})

test_that("a time-varying destination covariate is simulated and recovered", {
  skip_on_cran()
  tvf <- function(J, years) {
    a <- array(runif(J * length(years), -2, 2), dim = c(J, length(years), 1),
               dimnames = list(NULL, as.character(years), "bait"))
    a
  }
  sim <- simulate_dsf(~ 1, ~ 0 + bait, beta_sigma = log(5), beta_theta = 0.7,
                      ncell = 100, trap_tv = tvf, exclude_natal = TRUE, seed = 13)
  expect_equal(dim(sim$trap_tv)[3], 1L)
  fit <- dsf(~ 1, ~ 0 + bait, sim$ind_df, sim$trap_df, sim$op_mat,
             trap_tv = sim$trap_tv, exclude_natal = TRUE)
  expect_equal(fit$convergence, 0)
  expect_true(abs(fit$coef["theta:bait"] - 0.7) < 3 * fit$se["theta:bait"])
  expect_true(abs(fit$coef["sigma:(Intercept)"] - log(5)) < 3 * fit$se["sigma:(Intercept)"])
  ## long-format trap_tv gives the same frame
  long <- data.frame(trap = rep(seq_len(dim(sim$trap_tv)[1]), times = dim(sim$trap_tv)[2]),
                     year = rep(dimnames(sim$trap_tv)[[2]], each = dim(sim$trap_tv)[1]),
                     bait = as.vector(sim$trap_tv[, , "bait"]))
  fr_long <- dsf_frame(sim$ind_df, sim$trap_df, sim$op_mat, trap_tv = long, exclude_natal = TRUE)
  expect_equal(fr_long$Z$bait, fit$frame$Z$bait)
})

test_that("time-varying effort enters as an offset", {
  sim <- simulate_dsf(~ 1, ~ 1, beta_sigma = log(5), ncell = 100, seed = 14)
  yrs <- colnames(sim$op_mat); J <- nrow(sim$trap_df)
  eff <- array(exp(rnorm(J * length(yrs))), dim = c(J, length(yrs), 1),
               dimnames = list(NULL, yrs, "effort"))
  fr <- dsf_frame(sim$ind_df, sim$trap_df, sim$op_mat, trap_tv = eff)
  expect_equal(fr$offset[, 1], log(eff[, 1, "effort"]))
  eff[1, 1, "effort"] <- -1
  expect_error(dsf_frame(sim$ind_df, sim$trap_df, sim$op_mat, trap_tv = eff), "positive")
})

test_that("landscape_seed fixes the design but not the dispersal draw", {
  a <- simulate_dsf(~ 1, ~ 1, beta_sigma = log(5), ncell = 100, landscape_seed = 7)
  b <- simulate_dsf(~ 1, ~ 1, beta_sigma = log(5), ncell = 100, landscape_seed = 7)
  expect_identical(a$trap_df, b$trap_df)
  expect_identical(a$ind_df[, c("src_x", "src_y", "obs_year")],
                   b$ind_df[, c("src_x", "src_y", "obs_year")])
  expect_false(identical(a$ind_df$dest, b$ind_df$dest))
})

test_that("non-Euclidean DSF recovers alpha (coarse grid)", {
  skip_on_cran()
  skip_if_not_installed("igraph")
  sim <- simulate_dsf(~ 1, ~ 1, beta_sigma = log(5),
                      ncell = 50, surface = TRUE, aRange = 10,
                      cost = ~ cov, beta_cost = 1,
                      exclude_natal = TRUE, seed = 5)
  fit <- dsf(~ 1, ~ 1, sim$ind_df, sim$trap_df, sim$op_mat, exclude_natal = TRUE,
             cost = ~ cov, landscape = sim$landscape)
  expect_equal(fit$convergence, 0)
  expect_true(abs(fit$coef["cost:cov"] - 1) < 3 * fit$se["cost:cov"])
  expect_true(abs(fit$coef["sigma:(Intercept)"] - log(5)) < 3 * fit$se["sigma:(Intercept)"])
})

test_that("alpha = 0 cost distance approximates Euclidean", {
  skip_if_not_installed("igraph")
  sim <- simulate_dsf(~ 1, ~ 1, beta_sigma = log(5), ncell = 100, surface = TRUE, seed = 6)
  cf <- make_cost_fun(sim$landscape, ~ cov,
                      sim$ind_df[, c("src_x", "src_y")], sim$trap_df[, c("x", "y")])
  d_cost <- sqrt(cf$fun(0))
  d_euc  <- sqrt(outer(sim$ind_df$src_x, sim$trap_df$x, "-")^2 +
                 outer(sim$ind_df$src_y, sim$trap_df$y, "-")^2)
  ## 8-neighbour grid paths are at most ~8% longer than straight lines; cell
  ## snapping adds up to one cell of slack, so only test well-separated pairs
  ratio <- (d_cost / d_euc)[d_euc > 10]
  expect_true(all(ratio > 0.9 & ratio < 1.15))
})

test_that("cost raster must contain all points and have no NA", {
  skip_if_not_installed("igraph")
  sim <- simulate_dsf(~ 1, ~ 1, beta_sigma = log(5), ncell = 50, surface = TRUE, seed = 15)
  small <- terra::crop(sim$landscape, terra::ext(20, 80, 20, 80))
  expect_error(dsf(~ 1, ~ 1, sim$frame, cost = ~ cov, landscape = small), "outside")
  holey <- sim$landscape; holey[1] <- NA
  expect_error(dsf(~ 1, ~ 1, sim$frame, cost = ~ cov, landscape = holey), "missing values")
})

test_that("trap_cov supplies destination covariates and is part of the fixed design", {
  f <- function(n) runif(n, -2, 2)
  a <- simulate_dsf(~ 1, ~ 0 + cov, beta_sigma = log(5), beta_theta = 0.5,
                    ncell = 100, trap_cov = f, landscape_seed = 9)
  b <- simulate_dsf(~ 1, ~ 0 + cov, beta_sigma = log(5), beta_theta = 0.5,
                    ncell = 100, trap_cov = f, landscape_seed = 9)
  expect_true(all(a$trap_df$cov >= -2 & a$trap_df$cov <= 2))
  expect_identical(a$trap_df$cov, b$trap_df$cov)
  ## data frame with a new column name
  d <- simulate_dsf(~ 1, ~ 0 + bait, beta_sigma = log(5), beta_theta = 1, ncell = 100,
                    trap_cov = data.frame(bait = rbinom(250, 1, 0.5)), seed = 10)
  expect_true("bait" %in% names(d$trap_df))
  expect_error(simulate_dsf(~ 1, ~ 1, beta_sigma = log(5), ncell = 100, trap_cov = 1:3),
               "one row")
})

test_that("grid trap layout gives the requested density and respects buffer = 0", {
  sim <- simulate_dsf(~ 1, ~ 1, beta_sigma = log(5), ntraps = 16, L = 40, ncell = 200,
                      buffer = 0, trap_layout = "grid", births = "uniform", seed = 8)
  expect_equal(nrow(sim$trap_df), 16)
  expect_true(all(sim$trap_df$x > 0 & sim$trap_df$x < 40))
  ## regular lattice: nearest-neighbour spacing exactly L / 4 = 10 km
  dd <- as.matrix(dist(sim$trap_df[, c("x", "y")])); diag(dd) <- Inf
  expect_equal(unname(apply(dd, 1, min)), rep(10, 16))
  ## jittered lattice: close to it
  simj <- simulate_dsf(~ 1, ~ 1, beta_sigma = log(5), ntraps = 16, L = 40, ncell = 200,
                       buffer = 0, trap_layout = "grid", jitter = 0.25,
                       births = "uniform", seed = 8)
  ddj <- as.matrix(dist(simj$trap_df[, c("x", "y")])); diag(ddj) <- Inf
  expect_true(abs(mean(apply(ddj, 1, min)) - 10) < 2.5)
})
