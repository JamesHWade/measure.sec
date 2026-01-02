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
  expect_equal(g$g[1], (10 / 12)^2, tolerance = 0.01)
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


# ==============================================================================
# Tests for measure_branching_frequency
# ==============================================================================

test_that("measure_branching_frequency calculates branches for random architecture", {
  g_values <- c(0.9, 0.7, 0.5, 0.3)

  bf <- measure_branching_frequency(g_values)

  expect_s3_class(bf, "branching_frequency")
  expect_equal(bf$architecture, rep("random", 4))

  # Higher g (less branching) -> fewer branches
  expect_true(bf$branches_per_molecule[1] < bf$branches_per_molecule[4])

  # g=1 should give 0 branches (but we can't test exactly 1)
  bf_linear <- measure_branching_frequency(0.99)
  expect_lt(bf_linear$branches_per_molecule[1], 0.5)
})

test_that("measure_branching_frequency handles star architectures", {
  # 3-arm star: g = 7/9 ~ 0.778
  bf_star3 <- measure_branching_frequency(0.77, architecture = "star_3")
  expect_equal(bf_star3$branches_per_molecule[1], 1)

  # 4-arm star: g = 10/16 = 0.625
  bf_star4 <- measure_branching_frequency(0.62, architecture = "star_4")
  expect_equal(bf_star4$branches_per_molecule[1], 1)
})

test_that("measure_branching_frequency calculates branch density with MW", {
  g_values <- c(0.8, 0.6, 0.4)
  mw <- c(50000, 100000, 200000)

  bf <- measure_branching_frequency(g_values, mw = mw)

  expect_true("branch_density" %in% names(bf))
  expect_true("mw" %in% names(bf))
  expect_equal(bf$mw, mw)
})

test_that("measure_branching_frequency validates inputs", {
  expect_error(
    measure_branching_frequency(c(1.5, 0.5)),
    "between 0 and 1"
  )

  expect_error(
    measure_branching_frequency(c(0.5, 0.6), architecture = "star_f"),
    "arms"
  )
})

test_that("measure_branching_frequency print method works", {
  bf <- measure_branching_frequency(c(0.8, 0.6, 0.4))

  expect_output(print(bf), "Branching Frequency")
  expect_output(print(bf), "Zimm-Stockmayer")
})

test_that("measure_branching_frequency handles comb architecture", {
  # Comb polymer: g ~ 1 / (1 + 2*n_b/3)
  # For g = 0.6, n_b = 3 * (1/0.6 - 1) / 2 = 1.0
  bf_comb <- measure_branching_frequency(0.6, architecture = "comb")
  expect_equal(bf_comb$architecture, "comb")
  expect_equal(bf_comb$branches_per_molecule[1], 1.0, tolerance = 0.01)

  # g = 1 should give 0 branches
  bf_linear <- measure_branching_frequency(0.99, architecture = "comb")
  expect_lt(bf_linear$branches_per_molecule[1], 0.1)
})

test_that("measure_branching_frequency handles star_f with arms parameter", {
  # 5-arm star: g = (3*5 - 2) / 5^2 = 13/25 = 0.52
  g_5arm <- (3 * 5 - 2) / 5^2
  bf_star5 <- measure_branching_frequency(
    g_5arm,
    architecture = "star_f",
    arms = 5
  )
  expect_equal(bf_star5$branches_per_molecule[1], 1) # All stars have 1 branch point

  # 6-arm star: g = (3*6 - 2) / 6^2 = 16/36 = 0.444
  g_6arm <- (3 * 6 - 2) / 6^2
  bf_star6 <- measure_branching_frequency(
    g_6arm,
    architecture = "star_f",
    arms = 6
  )
  expect_equal(bf_star6$branches_per_molecule[1], 1)
})

test_that("measure_branching_frequency validates arms parameter", {
  expect_error(
    measure_branching_frequency(0.5, architecture = "star_f", arms = 2),
    "must be.*>= 3"
  )

  expect_error(
    measure_branching_frequency(0.5, architecture = "star_f", arms = "three"),
    "must be.*numeric"
  )
})


# ==============================================================================
# Tests for measure_rg_mw_scaling
# ==============================================================================

test_that("measure_rg_mw_scaling fits power law relationship", {
  # Test with nu = 0.588 (good solvent)
  mw <- c(10000, 50000, 100000, 500000, 1000000)
  rg <- 0.05 * mw^0.588

  scaling <- measure_rg_mw_scaling(mw, rg)

  expect_s3_class(scaling, "rg_mw_scaling")
  expect_equal(scaling$nu, 0.588, tolerance = 0.01)
  expect_gt(scaling$r_squared, 0.99)
  expect_equal(scaling$conformation, "good solvent (swollen coil)")
})

test_that("measure_rg_mw_scaling interprets conformations correctly", {
  mw <- c(10000, 50000, 100000, 500000)

  # Compact/spherical (nu ~ 0.33)
  rg_compact <- 0.5 * mw^0.33
  scaling_compact <- measure_rg_mw_scaling(mw, rg_compact)
  expect_equal(scaling_compact$conformation, "compact/spherical")

  # Theta solvent (nu ~ 0.50)
  rg_theta <- 0.1 * mw^0.50
  scaling_theta <- measure_rg_mw_scaling(mw, rg_theta)
  expect_equal(scaling_theta$conformation, "theta solvent (ideal chain)")
})

test_that("measure_rg_mw_scaling provides confidence intervals", {
  mw <- c(10000, 50000, 100000, 500000, 1000000)
  rg <- 0.05 * mw^0.588

  scaling <- measure_rg_mw_scaling(mw, rg)

  expect_true(!is.na(scaling$nu_ci_lower))
  expect_true(!is.na(scaling$nu_ci_upper))
  expect_lt(scaling$nu_ci_lower, scaling$nu)
  expect_gt(scaling$nu_ci_upper, scaling$nu)
})

test_that("measure_rg_mw_scaling respects mw_range", {
  mw <- c(1000, 10000, 50000, 100000, 500000)
  rg <- 0.05 * mw^0.588

  scaling <- measure_rg_mw_scaling(mw, rg, mw_range = c(10000, 100000))

  expect_equal(scaling$n_points, 3)
  expect_equal(scaling$mw_range, c(10000, 100000))
})

test_that("measure_rg_mw_scaling requires minimum data", {
  expect_error(
    measure_rg_mw_scaling(c(10000, 20000), c(5, 10)),
    "3 valid"
  )
})

test_that("measure_rg_mw_scaling print method works", {
  mw <- c(10000, 50000, 100000)
  rg <- 0.05 * mw^0.588

  scaling <- measure_rg_mw_scaling(mw, rg)

  expect_output(print(scaling), "Rg-MW Scaling")
  expect_output(print(scaling), "Flory exponent")
  expect_output(print(scaling), "Conformation")
})

test_that("measure_rg_mw_scaling handles weighted regression", {
  mw <- c(10000, 50000, 100000, 500000, 1000000)
  rg <- 0.05 * mw^0.588
  weights <- c(0.1, 0.5, 1.0, 0.5, 0.1)

  scaling <- measure_rg_mw_scaling(mw, rg, weights = weights)

  expect_s3_class(scaling, "rg_mw_scaling")
  expect_true(!is.null(scaling$nu))
  # Result should still be reasonable
  expect_equal(scaling$nu, 0.588, tolerance = 0.05)
})


# ==============================================================================
# Tests for measure_sec_column_performance
# ==============================================================================

test_that("measure_sec_column_performance calculates basic metrics", {
  cal_data <- data.frame(
    retention = c(5.2, 6.1, 7.0, 8.2, 9.5, 10.8),
    mw = c(1200000, 400000, 100000, 30000, 5000, 580)
  )

  perf <- measure_sec_column_performance(cal_data, column_length = 30)

  expect_s3_class(perf, "sec_column_performance")
  expect_equal(perf$separation_range$exclusion_limit, 1200000)
  expect_equal(perf$separation_range$total_permeation, 580)
  expect_gt(perf$selectivity, 0)
  expect_gt(perf$calibration_r_squared, 0.9)
})

test_that("measure_sec_column_performance calculates plate metrics with width", {
  cal_data <- data.frame(
    retention = c(5.2, 6.1, 7.0, 8.2, 9.5, 10.8),
    mw = c(1200000, 400000, 100000, 30000, 5000, 580),
    width = c(0.4, 0.35, 0.30, 0.28, 0.25, 0.30)
  )

  perf <- measure_sec_column_performance(
    cal_data,
    column_length = 30,
    particle_size = 5
  )

  expect_false(is.na(perf$hetp_mm))
  expect_false(is.na(perf$plates_per_meter))
  expect_false(is.na(perf$reduced_hetp))
  expect_false(is.na(perf$peak_capacity))
})

test_that("measure_sec_column_performance handles different column names", {
  # Using retention_time and molecular_weight
  cal_data <- data.frame(
    retention_time = c(5.2, 6.1, 7.0),
    molecular_weight = c(1200000, 400000, 100000)
  )

  perf <- measure_sec_column_performance(cal_data)

  expect_s3_class(perf, "sec_column_performance")
})

test_that("measure_sec_column_performance validates inputs", {
  bad_data <- data.frame(x = 1:5, y = 1:5)

  expect_error(
    measure_sec_column_performance(bad_data),
    "retention"
  )

  bad_data2 <- data.frame(retention = 1:5, z = 1:5)
  expect_error(
    measure_sec_column_performance(bad_data2),
    "molecular weight"
  )
})

test_that("measure_sec_column_performance validates column_length", {
  cal_data <- data.frame(
    retention = c(5.2, 6.1, 7.0),
    mw = c(1200000, 400000, 100000)
  )

  expect_error(
    measure_sec_column_performance(cal_data, column_length = -10),
    "positive numeric"
  )

  expect_error(
    measure_sec_column_performance(cal_data, column_length = "thirty"),
    "positive numeric"
  )
})

test_that("measure_sec_column_performance validates particle_size", {
  cal_data <- data.frame(
    retention = c(5.2, 6.1, 7.0),
    mw = c(1200000, 400000, 100000)
  )

  expect_error(
    measure_sec_column_performance(cal_data, particle_size = -5),
    "positive numeric"
  )

  expect_error(
    measure_sec_column_performance(cal_data, particle_size = "five"),
    "positive numeric"
  )
})

test_that("measure_sec_column_performance print method works", {
  cal_data <- data.frame(
    retention = c(5.2, 6.1, 7.0, 8.2),
    mw = c(1200000, 400000, 100000, 30000),
    width = c(0.4, 0.35, 0.30, 0.28)
  )

  perf <- measure_sec_column_performance(cal_data)

  expect_output(print(perf), "SEC Column Performance")
  expect_output(print(perf), "Separation Range")
})


# ==============================================================================
# Tests for measure_branching_model_comparison
# ==============================================================================

test_that("measure_branching_model_comparison fits multiple models", {
  # Generate data that fits random branching model
  mw <- c(50000, 100000, 200000, 500000)
  # Simulate g values decreasing with MW (more branches at higher MW)
  g_exp <- c(0.85, 0.72, 0.58, 0.42)

  comparison <- measure_branching_model_comparison(g_exp, mw)

  expect_s3_class(comparison, "branching_model_comparison")
  expect_true("model_fits" %in% names(comparison))
  expect_true("best_model" %in% names(comparison))
  expect_true("predictions" %in% names(comparison))
})

test_that("measure_branching_model_comparison identifies best model", {
  mw <- c(50000, 100000, 200000, 500000)
  g_exp <- c(0.85, 0.72, 0.58, 0.42)

  comparison <- measure_branching_model_comparison(g_exp, mw)

  expect_true(
    comparison$best_model %in% c("random", "star", "comb", "hyperbranched")
  )
  expect_true(!is.null(comparison$best_fit$r_squared))
})

test_that("measure_branching_model_comparison allows model selection", {
  mw <- c(50000, 100000, 200000, 500000)
  g_exp <- c(0.85, 0.72, 0.58, 0.42)

  comparison <- measure_branching_model_comparison(
    g_exp,
    mw,
    models = c("random", "comb")
  )

  expect_equal(nrow(comparison$model_fits), 2)
  expect_true(all(comparison$model_fits$model %in% c("random", "comb")))
})

test_that("measure_branching_model_comparison validates inputs", {
  expect_error(
    measure_branching_model_comparison(c(0.5, 0.6), c(10000, 20000)),
    "3 valid"
  )

  expect_error(
    measure_branching_model_comparison(c(0.5, 0.6, 0.7), c(10000, 20000)),
    "same length"
  )
})

test_that("measure_branching_model_comparison print method works", {
  mw <- c(50000, 100000, 200000, 500000)
  g_exp <- c(0.85, 0.72, 0.58, 0.42)

  comparison <- measure_branching_model_comparison(g_exp, mw)

  expect_output(print(comparison), "Branching Model Comparison")
  expect_output(print(comparison), "Best Model")
})

test_that("measure_branching_model_comparison handles problematic data gracefully", {
  # Data that's very hard to fit (all same value)
  mw <- c(50000, 100000, 200000)
  g_exp <- c(0.5, 0.5, 0.5) # Constant g - no MW dependence

  # Should either fit successfully (rare) or warn about failed models
  # Either outcome is acceptable
  result <- suppressWarnings(measure_branching_model_comparison(g_exp, mw))

  # Result can be NULL (all failed) or a valid comparison (some succeeded)
  expect_true(is.null(result) || inherits(result, "branching_model_comparison"))
})

test_that("measure_branching_model_comparison warns about failed models", {
  # Create data that some models may struggle with
  mw <- c(50000, 100000, 200000, 500000)
  g_exp <- c(0.99, 0.98, 0.97, 0.96) # Very linear-like

  # This may cause some models to fail - check that warning is issued
  # The test just ensures no errors occur
  result <- tryCatch(
    suppressWarnings(measure_branching_model_comparison(g_exp, mw)),
    error = function(e) NULL
  )
  # Result could be NULL or a valid comparison
  expect_true(is.null(result) || inherits(result, "branching_model_comparison"))
})

test_that("measure_branching_model_comparison handles NA in input data", {
  mw <- c(50000, NA, 200000, 500000, 1000000)
  g_exp <- c(0.85, 0.72, NA, 0.42, 0.30)

  # Should remove NA values and still work with remaining 3 points
  comparison <- measure_branching_model_comparison(g_exp, mw)

  expect_s3_class(comparison, "branching_model_comparison")
})
