# ==============================================================================
# Tests for molecular weight calculation steps
# ==============================================================================

test_that("step_sec_mw_averages calculates Mn, Mw, Mz correctly", {
  skip_if_not_installed("measure")

  # Create test data with known MW distribution
  # Use a unimodal distribution for testing
  # Assume x-axis is log10(MW)
  log_mw <- seq(3, 6, by = 0.05) # log10(MW) from 1000 to 1M
  n <- length(log_mw)

  # Concentration profile (Gaussian)
  conc <- dnorm(log_mw, mean = 4.5, sd = 0.5)
  conc <- conc / sum(conc) # Normalize

  test_data <- tibble::tibble(sample_id = "test")

  # Create measure object with log10(MW) as location and concentration as value
  test_data$ri <- measure::new_measure_list(
    list(measure::new_measure_tbl(location = log_mw, value = conc))
  )

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_mw_averages(measures = "ri")

  prepped <- recipes::prep(rec)
  result <- recipes::bake(prepped, new_data = NULL)

  expect_s3_class(result, "tbl_df")
  expect_true("mw_mn" %in% names(result))
  expect_true("mw_mw" %in% names(result))
  expect_true("mw_mz" %in% names(result))
  expect_true("mw_dispersity" %in% names(result))

  # Check that Mw > Mn (always true for polydisperse samples)
  expect_gt(result$mw_mw[1], result$mw_mn[1])

  # Check that Mz > Mw
  expect_gt(result$mw_mz[1], result$mw_mw[1])

  # Check dispersity = Mw/Mn
  expect_equal(
    result$mw_dispersity[1],
    result$mw_mw[1] / result$mw_mn[1],
    tolerance = 0.01
  )
})

test_that("step_sec_mw_averages respects integration_range", {
  skip_if_not_installed("measure")

  log_mw <- seq(3, 6, by = 0.05)
  conc <- dnorm(log_mw, mean = 4.5, sd = 0.5)

  test_data <- tibble::tibble(sample_id = "test")

  test_data$ri <- measure::new_measure_list(
    list(measure::new_measure_tbl(location = log_mw, value = conc))
  )

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_mw_averages(measures = "ri", integration_range = c(4.0, 5.0))

  prepped <- recipes::prep(rec)
  result <- recipes::bake(prepped, new_data = NULL)

  expect_s3_class(result, "tbl_df")
  expect_true("mw_mn" %in% names(result))
})

test_that("step_sec_mw_fractions calculates fractions correctly", {
  skip_if_not_installed("measure")

  log_mw <- seq(3, 6, by = 0.05)
  conc <- dnorm(log_mw, mean = 4.5, sd = 0.5)

  test_data <- tibble::tibble(sample_id = "test")

  test_data$ri <- measure::new_measure_list(
    list(measure::new_measure_tbl(location = log_mw, value = conc))
  )

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_mw_fractions(measures = "ri", cutoffs = c(10000, 100000))

  prepped <- recipes::prep(rec)
  result <- recipes::bake(prepped, new_data = NULL)

  expect_s3_class(result, "tbl_df")

  # Should have fraction columns
  frac_cols <- names(result)[grepl("frac", names(result), fixed = TRUE)]
  expect_true(length(frac_cols) > 0)
})

test_that("step_sec_mw_distribution creates differential distribution", {
  skip_if_not_installed("measure")

  log_mw <- seq(3, 6, by = 0.05)
  conc <- dnorm(log_mw, mean = 4.5, sd = 0.5)

  test_data <- tibble::tibble(sample_id = "test")

  test_data$ri <- measure::new_measure_list(
    list(measure::new_measure_tbl(location = log_mw, value = conc))
  )

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_mw_distribution(measures = "ri", type = "differential")

  prepped <- recipes::prep(rec)
  result <- recipes::bake(prepped, new_data = NULL)

  expect_s3_class(result, "tbl_df")
  # The distribution replaces the measure column with distribution data
  expect_true("ri" %in% names(result))
  # Check that the result contains measure data
  expect_s3_class(result$ri, "measure_list")
})

test_that("step_sec_baseline corrects baseline", {
  skip_if_not_installed("measure")

  time <- seq(5, 25, by = 0.1)

  # Signal with baseline drift
  peak <- dnorm(time, mean = 15, sd = 1)
  baseline <- 0.1 + 0.005 * time # Linear drift
  signal <- peak + baseline

  test_data <- tibble::tibble(sample_id = "test")

  test_data$ri <- measure::new_measure_list(
    list(measure::new_measure_tbl(location = time, value = signal))
  )

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_baseline(
      measures = "ri",
      left_frac = 0.1,
      right_frac = 0.1
    )

  prepped <- recipes::prep(rec)
  result <- recipes::bake(prepped, new_data = NULL)

  expect_s3_class(result, "tbl_df")

  # Baseline should be approximately removed at edges
  corrected <- result$ri[[1]]$value
  expect_lt(abs(corrected[1]), abs(signal[1]))
})

test_that("step_sec_mw_averages validates output_cols", {
  skip_if_not_installed("measure")

  log_mw <- seq(3, 6, by = 0.05)
  conc <- dnorm(log_mw, mean = 4.5, sd = 0.5)

  test_data <- tibble::tibble(sample_id = "test")

  test_data$ri <- measure::new_measure_list(
    list(measure::new_measure_tbl(location = log_mw, value = conc))
  )

  expect_error(
    recipes::recipe(~., data = test_data) |>
      step_sec_mw_averages(measures = "ri", output_cols = c("invalid")),
    "Invalid output columns"
  )
})

# -- Uncertainty propagation tests ----------------------------------------------

test_that("step_sec_mw_averages requires calibration_error when include_uncertainty is TRUE", {
  skip_if_not_installed("measure")

  log_mw <- seq(3, 6, by = 0.05)
  conc <- dnorm(log_mw, mean = 4.5, sd = 0.5)

  test_data <- tibble::tibble(sample_id = "test")
  test_data$ri <- measure::new_measure_list(
    list(measure::new_measure_tbl(location = log_mw, value = conc))
  )

  expect_error(
    recipes::recipe(~., data = test_data) |>
      step_sec_mw_averages(measures = "ri", include_uncertainty = TRUE),
    "calibration_error.*required"
  )
})

test_that("step_sec_mw_averages validates calibration_error", {
  skip_if_not_installed("measure")

  log_mw <- seq(3, 6, by = 0.05)
  conc <- dnorm(log_mw, mean = 4.5, sd = 0.5)

  test_data <- tibble::tibble(sample_id = "test")
  test_data$ri <- measure::new_measure_list(
    list(measure::new_measure_tbl(location = log_mw, value = conc))
  )

  # Must be numeric
  expect_error(
    recipes::recipe(~., data = test_data) |>
      step_sec_mw_averages(
        measures = "ri",
        include_uncertainty = TRUE,
        calibration_error = "not_numeric"
      ),
    "single numeric value"
  )

  # Must be positive
  expect_error(
    recipes::recipe(~., data = test_data) |>
      step_sec_mw_averages(
        measures = "ri",
        include_uncertainty = TRUE,
        calibration_error = -0.01
      ),
    "must be positive"
  )
})

test_that("step_sec_mw_averages calculates uncertainties with include_uncertainty = TRUE", {
  skip_if_not_installed("measure")

  log_mw <- seq(3, 6, by = 0.05)
  conc <- dnorm(log_mw, mean = 4.5, sd = 0.5)

  test_data <- tibble::tibble(sample_id = "test")
  test_data$ri <- measure::new_measure_list(
    list(measure::new_measure_tbl(location = log_mw, value = conc))
  )

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_mw_averages(
      measures = "ri",
      include_uncertainty = TRUE,
      calibration_error = 0.02 # Typical RMSE in log10(MW)
    )

  prepped <- recipes::prep(rec)
  result <- recipes::bake(prepped, new_data = NULL)

  # Should have uncertainty columns
  expect_true("mw_mn_uncertainty" %in% names(result))
  expect_true("mw_mw_uncertainty" %in% names(result))
  expect_true("mw_mz_uncertainty" %in% names(result))
  expect_true("mw_dispersity_uncertainty" %in% names(result))

  # Uncertainties should be positive
  expect_gt(result$mw_mn_uncertainty[1], 0)
  expect_gt(result$mw_mw_uncertainty[1], 0)
  expect_gt(result$mw_mz_uncertainty[1], 0)
  expect_gt(result$mw_dispersity_uncertainty[1], 0)

  # Uncertainties should be smaller than the values themselves (reasonable)
  expect_lt(result$mw_mn_uncertainty[1], result$mw_mn[1])
  expect_lt(result$mw_mw_uncertainty[1], result$mw_mw[1])
  expect_lt(result$mw_mz_uncertainty[1], result$mw_mz[1])
})

test_that("step_sec_mw_averages uncertainty scales with calibration error", {
  skip_if_not_installed("measure")

  log_mw <- seq(3, 6, by = 0.05)
  conc <- dnorm(log_mw, mean = 4.5, sd = 0.5)

  test_data <- tibble::tibble(sample_id = "test")
  test_data$ri <- measure::new_measure_list(
    list(measure::new_measure_tbl(location = log_mw, value = conc))
  )

  # Small calibration error
  rec_small <- recipes::recipe(~., data = test_data) |>
    step_sec_mw_averages(
      measures = "ri",
      include_uncertainty = TRUE,
      calibration_error = 0.01
    )

  # Larger calibration error
  rec_large <- recipes::recipe(~., data = test_data) |>
    step_sec_mw_averages(
      measures = "ri",
      include_uncertainty = TRUE,
      calibration_error = 0.05
    )

  result_small <- recipes::bake(recipes::prep(rec_small), new_data = NULL)
  result_large <- recipes::bake(recipes::prep(rec_large), new_data = NULL)

  # Larger calibration error should give larger uncertainties
  expect_gt(
    result_large$mw_mn_uncertainty[1],
    result_small$mw_mn_uncertainty[1]
  )
  expect_gt(
    result_large$mw_mw_uncertainty[1],
    result_small$mw_mw_uncertainty[1]
  )
})

test_that("step_sec_mw_averages tidy includes uncertainty info", {
  skip_if_not_installed("measure")

  log_mw <- seq(3, 6, by = 0.05)
  conc <- dnorm(log_mw, mean = 4.5, sd = 0.5)

  test_data <- tibble::tibble(sample_id = "test")
  test_data$ri <- measure::new_measure_list(
    list(measure::new_measure_tbl(location = log_mw, value = conc))
  )

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_mw_averages(
      measures = "ri",
      include_uncertainty = TRUE,
      calibration_error = 0.02
    )

  prepped <- recipes::prep(rec)
  tidy_result <- recipes::tidy(prepped, number = 1)

  expect_true("include_uncertainty" %in% names(tidy_result))
  expect_true("calibration_error" %in% names(tidy_result))
  expect_equal(tidy_result$include_uncertainty, TRUE)
  expect_equal(tidy_result$calibration_error, 0.02)
})

test_that("step_sec_mw_averages print method shows uncertainty status", {
  skip_if_not_installed("measure")

  log_mw <- seq(3, 6, by = 0.05)
  conc <- dnorm(log_mw, mean = 4.5, sd = 0.5)

  test_data <- tibble::tibble(sample_id = "test")
  test_data$ri <- measure::new_measure_list(
    list(measure::new_measure_tbl(location = log_mw, value = conc))
  )

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_mw_averages(
      measures = "ri",
      include_uncertainty = TRUE,
      calibration_error = 0.02
    )

  expect_output(print(rec), "uncertainties")
})

# -- RF baseline tests ---------------------------------------------------------

test_that("step_sec_baseline with method='rf' corrects curved baseline", {
  skip_if_not_installed("measure")
  skip_if_not_installed("IDPmisc")

  time <- seq(5, 25, by = 0.1)

  # Signal with curved baseline
  peak <- dnorm(time, mean = 15, sd = 1)
  baseline <- 0.1 + 0.002 * (time - 15)^2  # Curved drift
  signal <- peak + baseline

  test_data <- tibble::tibble(sample_id = "test")

  test_data$ri <- measure::new_measure_list(
    list(measure::new_measure_tbl(location = time, value = signal))
  )

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_baseline(measures = "ri", method = "rf")

  prepped <- recipes::prep(rec)
  result <- recipes::bake(prepped, new_data = NULL)

  expect_s3_class(result, "tbl_df")

  # RF should remove baseline at edges
  corrected <- result$ri[[1]]$value
  # The edges should be closer to zero than the original signal
  expect_lt(abs(corrected[1]), abs(signal[1] - min(signal)))
})

test_that("step_sec_baseline RF method validates rf_span parameter", {
  skip_if_not_installed("measure")

  time <- seq(5, 25, by = 0.1)
  signal <- dnorm(time, mean = 15, sd = 1)

  test_data <- tibble::tibble(sample_id = "test")
  test_data$ri <- measure::new_measure_list(
    list(measure::new_measure_tbl(location = time, value = signal))
  )

  # Invalid rf_span (> 1)
  expect_error(
    recipes::recipe(~., data = test_data) |>
      step_sec_baseline(measures = "ri", method = "rf", rf_span = 1.5) |>
      recipes::prep(),
    "rf_span.*between 0 and 1"
  )

  # Invalid rf_span (negative)
  expect_error(
    recipes::recipe(~., data = test_data) |>
      step_sec_baseline(measures = "ri", method = "rf", rf_span = -0.5) |>
      recipes::prep(),
    "rf_span.*between 0 and 1"
  )
})

test_that("step_sec_baseline RF method validates rf_maxit parameter", {
  skip_if_not_installed("measure")

  time <- seq(5, 25, by = 0.1)
  signal <- dnorm(time, mean = 15, sd = 1)

  test_data <- tibble::tibble(sample_id = "test")
  test_data$ri <- measure::new_measure_list(
    list(measure::new_measure_tbl(location = time, value = signal))
  )

  # Invalid rf_maxit (wrong length)
  expect_error(
    recipes::recipe(~., data = test_data) |>
      step_sec_baseline(measures = "ri", method = "rf", rf_maxit = 20) |>
      recipes::prep(),
    "rf_maxit.*length 2"
  )
})

test_that("step_sec_baseline tidy includes RF parameters", {
  skip_if_not_installed("measure")

  time <- seq(5, 25, by = 0.1)
  signal <- dnorm(time, mean = 15, sd = 1)

  test_data <- tibble::tibble(sample_id = "test")
  test_data$ri <- measure::new_measure_list(
    list(measure::new_measure_tbl(location = time, value = signal))
  )

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_baseline(
      measures = "ri",
      method = "rf",
      rf_span = 0.5,
      rf_maxit = c(10, 10)
    )

  prepped <- recipes::prep(rec)
  tidy_result <- recipes::tidy(prepped, number = 1)

  expect_true("rf_span" %in% names(tidy_result))
  expect_true("rf_maxit" %in% names(tidy_result))
  expect_equal(tidy_result$rf_span, 0.5)
  expect_equal(tidy_result$rf_maxit[[1]], c(10, 10))
})

test_that("step_sec_baseline RF method skips fraction validation", {
  skip_if_not_installed("measure")
  skip_if_not_installed("IDPmisc")

  time <- seq(5, 25, by = 0.1)
  signal <- dnorm(time, mean = 15, sd = 1)

  test_data <- tibble::tibble(sample_id = "test")
  test_data$ri <- measure::new_measure_list(
    list(measure::new_measure_tbl(location = time, value = signal))
  )

  # This should not error even though left_frac + right_frac > 0.5
  # because RF method doesn't use these parameters
  expect_no_error(
    recipes::recipe(~., data = test_data) |>
      step_sec_baseline(
        measures = "ri",
        method = "rf",
        left_frac = 0.3,
        right_frac = 0.3
      ) |>
      recipes::prep()
  )
})

test_that("step_sec_baseline required_pkgs includes IDPmisc for RF method", {
  skip_if_not_installed("measure")

  time <- seq(5, 25, by = 0.1)
  signal <- dnorm(time, mean = 15, sd = 1)

  test_data <- tibble::tibble(sample_id = "test")
  test_data$ri <- measure::new_measure_list(
    list(measure::new_measure_tbl(location = time, value = signal))
  )

  # RF method requires IDPmisc
  rec_rf <- recipes::recipe(~., data = test_data) |>
    step_sec_baseline(measures = "ri", method = "rf")
  prepped_rf <- recipes::prep(rec_rf)

  pkgs_rf <- recipes::required_pkgs(prepped_rf)
  expect_true("IDPmisc" %in% pkgs_rf)

  # Linear method does not require IDPmisc
  rec_linear <- recipes::recipe(~., data = test_data) |>
    step_sec_baseline(measures = "ri", method = "linear")
  prepped_linear <- recipes::prep(rec_linear)

  pkgs_linear <- recipes::required_pkgs(prepped_linear)
  expect_false("IDPmisc" %in% pkgs_linear)
})
