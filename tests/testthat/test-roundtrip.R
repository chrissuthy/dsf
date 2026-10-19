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
  expect_error(dsf(~ 1, ~ 1, sim$ind_df, sim$trap_df, op), "unavailable")
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
                      ncell = 100, surface = TRUE, aRange = 10,
                      cost = ~ cov, beta_cost = 1, cost_res = 2,
                      exclude_natal = TRUE, seed = 5)
  fit <- dsf(~ 1, ~ 1, sim$ind_df, sim$trap_df, sim$op_mat, exclude_natal = TRUE,
             cost = ~ cov, landscape = sim$landscape, cost_res = 2)
  expect_equal(fit$convergence, 0)
  expect_true(abs(fit$coef["cost:cov"] - 1) < 3 * fit$se["cost:cov"])
  expect_true(abs(fit$coef["sigma:(Intercept)"] - log(5)) < 3 * fit$se["sigma:(Intercept)"])
})

test_that("alpha = 0 cost distance approximates Euclidean", {
  skip_if_not_installed("igraph")
  sim <- simulate_dsf(~ 1, ~ 1, beta_sigma = log(5), ncell = 100, surface = TRUE, seed = 6)
  cf <- make_cost_fun(sim$landscape, ~ cov, cost_res = 1,
                      sim$ind_df[, c("src_x", "src_y")], sim$trap_df[, c("x", "y")])
  d_cost <- sqrt(cf$fun(0))
  d_euc  <- sqrt(outer(sim$ind_df$src_x, sim$trap_df$x, "-")^2 +
                 outer(sim$ind_df$src_y, sim$trap_df$y, "-")^2)
  ## 8-neighbour grid paths are at most ~8% longer than straight lines; cell
  ## snapping adds up to one cell of slack, so only test well-separated pairs
  ratio <- (d_cost / d_euc)[d_euc > 10]
  expect_true(all(ratio > 0.9 & ratio < 1.15))
})

test_that("grid trap layout gives the requested density and respects buffer = 0", {
  sim <- simulate_dsf(~ 1, ~ 1, beta_sigma = log(5), ntraps = 16, L = 40, ncell = 200,
                      buffer = 0, trap_layout = "grid", births = "uniform", seed = 8)
  expect_equal(nrow(sim$trap_df), 16)
  expect_true(all(sim$trap_df$x > 0 & sim$trap_df$x < 40))
  ## jittered lattice: nearest-neighbour spacing close to L / 4 = 10 km
  dd <- as.matrix(dist(sim$trap_df[, c("x", "y")])); diag(dd) <- Inf
  expect_true(abs(mean(apply(dd, 1, min)) - 10) < 2.5)
})
