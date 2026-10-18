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
