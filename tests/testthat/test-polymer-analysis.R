# ==============================================================================
# Tests for polymer analysis functions
# ==============================================================================

test_that("measure_mh_parameters fits Mark-Houwink equation", {
  # Test with synthetic data: [eta] = 0.0001 * M^0.7
  mw <- c(10000, 25000, 50000, 100000, 250000)
  iv <- 0.0001 * mw^0.7

  mh <- measure_mh_parameters(mw, iv)

  expect_s3_class(mh, "mh_parameters")
  expect_equal(mh$K, 0.0001, tolerance = 0.0001)
  expect_equal(mh$a, 0.7, tolerance = 0.01)
  expect_gt(mh$r_squared, 0.99)
})

test_that("measure_mh_parameters respects MW range", {
  mw <- c(1000, 5000, 10000, 50000, 100000, 500000)
  iv <- 0.0001 * mw^0.7

  mh_full <- measure_mh_parameters(mw, iv)
  mh_subset <- measure_mh_parameters(mw, iv, mw_range = c(10000, 100000))

  expect_equal(mh_subset$n_points, 3)
  expect_equal(mh_subset$mw_range, c(10000, 100000))
})

test_that("measure_mh_parameters handles weighted regression", {
  mw <- c(10000, 25000, 50000, 100000, 250000)
  iv <- 0.0001 * mw^0.7
  weights <- c(0.1, 0.5, 1.0, 0.5, 0.1)

  mh <- measure_mh_parameters(mw, iv, weights = weights)

  expect_s3_class(mh, "mh_parameters")
  expect_true(!is.null(mh$K))
})

test_that("measure_mh_parameters requires minimum data points", {
  expect_error(
    measure_mh_parameters(c(10000, 20000), c(0.1, 0.2)),
    "3 valid"
  )
})

test_that("measure_mh_parameters print method works", {
  mw <- c(10000, 50000, 100000)
  iv <- 0.0001 * mw^0.7

  mh <- measure_mh_parameters(mw, iv)

  expect_output(print(mh), "Mark-Houwink")
  expect_output(print(mh), "K =")
  expect_output(print(mh), "a =")
})

test_that("measure_branching_index calculates g from Rg data", {
  mw <- c(100000, 200000, 500000)
  rg_branched <- c(10, 15, 25)
  rg_linear <- c(12, 18, 32)

  reference <- data.frame(mw = mw, rg = rg_linear)

  g <- measure_branching_index(
    mw = mw,
    rg = rg_branched,
    reference = reference,
    method = "g"
  )

  expect_s3_class(g, "branching_index")
  expect_true("g" %in% names(g))

  # g should be < 1 for branched (smaller Rg)
  expect_true(all(g$g < 1))

  # g = (10/12)^2, (15/18)^2, (25/32)^2
  expect_equal(g$g[1], (10/12)^2, tolerance = 0.01)
})

test_that("measure_branching_index calculates g_prime from IV data", {
  mw <- c(100000, 200000, 500000)
  iv_branched <- c(0.3, 0.4, 0.6)

  mh_linear <- list(K = 0.0001, a = 0.7)

  g <- measure_branching_index(
    mw = mw,
    intrinsic_visc = iv_branched,
    mh_linear = mh_linear,
    method = "g_prime"
  )

  expect_s3_class(g, "branching_index")
  expect_true("g_prime" %in% names(g))

  # Linear IV = K * M^a
  iv_linear_expected <- mh_linear$K * mw^mh_linear$a
  g_prime_expected <- iv_branched / iv_linear_expected

  expect_equal(g$g_prime, g_prime_expected, tolerance = 0.01)
})

test_that("measure_branching_index estimates branch points", {
  mw <- c(100000, 200000, 500000)
  rg_branched <- c(10, 15, 25)
  rg_linear <- c(12, 18, 32)

  g <- measure_branching_index(
    mw = mw,
    rg = rg_branched,
    reference = data.frame(mw = mw, rg = rg_linear),
    method = "g"
  )

  expect_true("branches_per_molecule" %in% names(g))
  expect_true(all(g$branches_per_molecule >= 0))
})

test_that("measure_branching_index requires appropriate data", {
  expect_error(
    measure_branching_index(mw = c(100000), method = "g"),
    "rg"
  )

  expect_error(
    measure_branching_index(
      mw = c(100000),
      intrinsic_visc = c(0.5),
      method = "g_prime"
    ),
    "mh_linear"
  )
})

test_that("measure_conformation_data prepares plot data", {
  mw <- c(10000, 50000, 100000, 500000)
  iv <- c(0.15, 0.35, 0.50, 0.95)

  plot_data <- measure_conformation_data(mw, iv, y_type = "iv")

  expect_s3_class(plot_data, "conformation_data")
  expect_true("log_mw" %in% names(plot_data))
  expect_true("log_y" %in% names(plot_data))
  expect_true("fitted" %in% names(plot_data))

  # Check log values
  expect_equal(plot_data$log_mw, log10(mw))
  expect_equal(plot_data$log_y, log10(iv))
})

test_that("measure_conformation_data includes fit attributes", {
  mw <- c(10000, 50000, 100000, 500000)
  iv <- c(0.15, 0.35, 0.50, 0.95)

  plot_data <- measure_conformation_data(mw, iv, fit_line = TRUE)

  expect_true(!is.null(attr(plot_data, "slope")))
  expect_true(!is.null(attr(plot_data, "intercept")))
  expect_true(!is.null(attr(plot_data, "r_squared")))
})

test_that("measure_conformation_data handles missing values", {
  mw <- c(10000, NA, 100000, 500000)
  iv <- c(0.15, 0.35, NA, 0.95)

  plot_data <- measure_conformation_data(mw, iv)

  # Should only have 2 valid points
  expect_equal(nrow(plot_data), 2)
})
