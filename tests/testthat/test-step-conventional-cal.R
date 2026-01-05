# ==============================================================================
# Tests for step_sec_conventional_cal
# ==============================================================================

# Helper to create test data with measure columns
create_test_sec_data <- function() {
  time <- seq(10, 20, by = 0.1)
  signal <- dnorm(time, mean = 15, sd = 1.5)

  test_data <- tibble::tibble(sample_id = "test")
  test_data$ri <- measure::new_measure_list(
    list(measure::new_measure_tbl(location = time, value = signal))
  )
  test_data
}

# Helper to create calibration standards
create_ps_standards <- function() {
  data.frame(
    retention = c(11.0, 12.5, 14.0, 15.5, 17.0, 18.5),
    log_mw = c(6.0, 5.5, 5.0, 4.5, 4.0, 3.5)
  )
}

# -- Constructor validation tests ----------------------------------------------

test_that("step_sec_conventional_cal requires standards", {
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()

  expect_error(
    recipes::recipe(~., data = test_data) |>
      step_sec_conventional_cal(),
    "standards.*required"
  )
})

test_that("step_sec_conventional_cal validates standards is a data frame", {
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()

  expect_error(
    recipes::recipe(~., data = test_data) |>
      step_sec_conventional_cal(standards = c(1, 2, 3)),
    "must be a data frame"
  )
})

test_that("step_sec_conventional_cal validates standards has location column", {
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()

  bad_standards <- data.frame(
    x = c(11, 12, 13, 14),
    log_mw = c(6, 5.5, 5, 4.5)
  )

  expect_error(
    recipes::recipe(~., data = test_data) |>
      step_sec_conventional_cal(standards = bad_standards),
    "location column"
  )
})

test_that("step_sec_conventional_cal validates standards has MW column", {
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()

  bad_standards <- data.frame(
    retention = c(11, 12, 13, 14),
    y = c(6, 5.5, 5, 4.5)
  )

  expect_error(
    recipes::recipe(~., data = test_data) |>
      step_sec_conventional_cal(standards = bad_standards),
    "molecular weight column"
  )
})

test_that("step_sec_conventional_cal validates sufficient standards for fit type", {
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()

  # Only 3 standards - not enough for cubic
  few_standards <- data.frame(
    retention = c(11, 14, 17),
    log_mw = c(6, 5, 4)
  )

  expect_error(
    recipes::recipe(~., data = test_data) |>
      step_sec_conventional_cal(
        standards = few_standards,
        fit_type = "cubic"
      ) |>
      recipes::prep(),
    "Insufficient standards"
  )

  # Should work with linear fit
  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_conventional_cal(standards = few_standards, fit_type = "linear")

  expect_no_error(suppressWarnings(recipes::prep(rec)))
})

# -- Prep tests ----------------------------------------------------------------

test_that("step_sec_conventional_cal fits calibration curve during prep", {
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()
  standards <- create_ps_standards()

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_conventional_cal(
      standards = standards,
      fit_type = "cubic",
      extrapolation = "none"
    )

  suppressWarnings(prepped <- recipes::prep(rec))

  # Check that calibration fit was created
  step <- prepped$steps[[1]]
  expect_true(step$trained)
  expect_true(inherits(step$calibration_fit, "lm"))
  expect_true(!is.null(step$calibration_diagnostics))
  expect_true(step$calibration_diagnostics$r_squared > 0.99)
  expect_equal(length(step$calibration_range), 2)

  # Check new diagnostics fields
  diag <- step$calibration_diagnostics
  expect_true(!is.null(diag$rmse_log_mw))
  expect_true(!is.null(diag$adj_r_squared))
  expect_true(!is.null(diag$residual_std_error))
  expect_true(!is.null(diag$standard_results))
  expect_s3_class(diag$standard_results, "tbl_df")
  expect_equal(nrow(diag$standard_results), 6) # 6 standards
})

test_that("step_sec_conventional_cal warns on poor fit quality", {
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()

  # Standards with poor linear relationship (non-monotonic)
  bad_standards <- data.frame(
    retention = c(11, 12, 13, 14, 15, 16),
    log_mw = c(5, 6, 4, 5.5, 3.5, 5.2) # Non-monotonic
  )

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_conventional_cal(standards = bad_standards, fit_type = "cubic")

  expect_warning(
    recipes::prep(rec),
    "fit quality is low"
  )
})

# -- Bake tests ----------------------------------------------------------------

test_that("step_sec_conventional_cal creates MW column on bake", {
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()
  standards <- create_ps_standards()

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_conventional_cal(
      standards = standards,
      output_col = "calibrated_mw",
      extrapolation = "none" # Avoid warnings and extrapolation issues
    )

  suppressWarnings(prepped <- recipes::prep(rec))
  result <- recipes::bake(prepped, new_data = NULL)

  expect_true("calibrated_mw" %in% names(result))
  expect_s3_class(result$calibrated_mw, "measure_list")

  # Check values are reasonable (log MW between 3 and 7) for non-NA values
  mw_values <- result$calibrated_mw[[1]]$value
  valid_values <- mw_values[!is.na(mw_values)]
  expect_true(length(valid_values) > 0)
  expect_true(all(valid_values >= 3 & valid_values <= 7))
})

test_that("step_sec_conventional_cal returns MW in Daltons when log_output = FALSE", {
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()
  standards <- create_ps_standards()

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_conventional_cal(
      standards = standards,
      log_output = FALSE,
      extrapolation = "none" # Avoid extrapolation issues
    )

  suppressWarnings(prepped <- recipes::prep(rec))
  result <- recipes::bake(prepped, new_data = NULL)

  # Values should be in Daltons (> 1000) for in-range data
  mw_values <- result$mw[[1]]$value
  valid_values <- mw_values[!is.na(mw_values)]
  expect_true(length(valid_values) > 0)
  expect_true(all(valid_values > 1000))
})

test_that("step_sec_conventional_cal warns for out-of-range data", {
  skip_if_not_installed("measure")

  # Create data that extends beyond calibration range
  time <- seq(5, 25, by = 0.1) # Extends beyond 11-18.5 calibration range
  signal <- dnorm(time, mean = 15, sd = 3)

  test_data <- tibble::tibble(sample_id = "test")
  test_data$ri <- measure::new_measure_list(
    list(measure::new_measure_tbl(location = time, value = signal))
  )

  standards <- create_ps_standards()

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_conventional_cal(standards = standards, extrapolation = "warn")

  # Warning is emitted during prep (when training data is baked)
  expect_warning(
    recipes::prep(rec),
    "outside calibration range"
  )
})

test_that("step_sec_conventional_cal returns NA for out-of-range when extrapolation='none'", {
  skip_if_not_installed("measure")

  # Create data that extends beyond calibration range
  time <- seq(5, 25, by = 0.1)
  signal <- dnorm(time, mean = 15, sd = 3)

  test_data <- tibble::tibble(sample_id = "test")
  test_data$ri <- measure::new_measure_list(
    list(measure::new_measure_tbl(location = time, value = signal))
  )

  standards <- create_ps_standards()

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_conventional_cal(standards = standards, extrapolation = "none")

  suppressWarnings(prepped <- recipes::prep(rec))
  result <- recipes::bake(prepped, new_data = NULL)

  mw_values <- result$mw[[1]]$value

  # Values at beginning (time < 11) should be NA
  expect_true(all(is.na(mw_values[time < 11])))

  # Values at end (time > 18.5) should be NA
  expect_true(all(is.na(mw_values[time > 18.5])))
})

# -- Different fit types -------------------------------------------------------

test_that("step_sec_conventional_cal works with different fit types", {
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()
  standards <- create_ps_standards()

  fit_types <- c("linear", "quadratic", "cubic")

  for (fit in fit_types) {
    rec <- recipes::recipe(~., data = test_data) |>
      step_sec_conventional_cal(
        standards = standards,
        fit_type = fit,
        extrapolation = "none"
      )

    suppressWarnings(prepped <- recipes::prep(rec))
    result <- recipes::bake(prepped, new_data = NULL)

    expect_true("mw" %in% names(result), info = paste("fit_type:", fit))
  }
})

# -- GAM fit type tests ---------------------------------------------------------

test_that("step_sec_conventional_cal works with GAM fit type", {
  skip_if_not_installed("measure")
  skip_if_not_installed("mgcv")

  test_data <- create_test_sec_data()
  standards <- create_ps_standards()

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_conventional_cal(
      standards = standards,
      fit_type = "gam",
      extrapolation = "none"
    )

  suppressWarnings(prepped <- recipes::prep(rec))

  # Check that GAM calibration fit was created
  step <- prepped$steps[[1]]
  expect_true(step$trained)
  expect_true(inherits(step$calibration_fit, "gam"))
  expect_true(!is.null(step$calibration_diagnostics))
  expect_true(step$calibration_diagnostics$r_squared > 0.99)

  # Bake should produce valid results

  result <- recipes::bake(prepped, new_data = NULL)
  expect_true("mw" %in% names(result))

  # Check values are reasonable (log MW between 3 and 7) for non-NA values
  mw_values <- result$mw[[1]]$value
  valid_values <- mw_values[!is.na(mw_values)]
  expect_true(length(valid_values) > 0)
  expect_true(all(valid_values >= 3 & valid_values <= 7))
})

test_that("step_sec_conventional_cal GAM rejects insufficient standards", {
  skip_if_not_installed("measure")
  skip_if_not_installed("mgcv")

  test_data <- create_test_sec_data()

  # Only 3 unique standards - GAM requires 4
  few_standards <- data.frame(
    retention = c(11, 14, 17),
    log_mw = c(6, 5, 4)
  )

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_conventional_cal(
      standards = few_standards,
      fit_type = "gam",
      extrapolation = "none"
    )

  # Should error due to insufficient standards (min 4 required for GAM)
  expect_error(
    recipes::prep(rec),
    "Insufficient standards"
  )
})

test_that("step_sec_conventional_cal GAM tidy method returns coefficients", {
  skip_if_not_installed("measure")
  skip_if_not_installed("mgcv")

  test_data <- create_test_sec_data()
  standards <- create_ps_standards()

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_conventional_cal(
      standards = standards,
      fit_type = "gam",
      extrapolation = "none"
    )

  suppressWarnings(prepped <- recipes::prep(rec))
  tidy_result <- recipes::tidy(prepped, number = 1)

  expect_s3_class(tidy_result, "tbl_df")
  expect_equal(tidy_result$fit_type, "gam")
  expect_true(tidy_result$r_squared > 0.99)
  expect_true(!is.null(tidy_result$coefficients[[1]]))
  expect_true(length(tidy_result$coefficients[[1]]) > 0)
})

test_that("step_sec_conventional_cal GAM requires mgcv package", {
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()
  standards <- create_ps_standards()

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_conventional_cal(
      standards = standards,
      fit_type = "gam",
      extrapolation = "none"
    )

  step <- rec$steps[[1]]
  required <- required_pkgs.step_sec_conventional_cal(step)

  expect_true("mgcv" %in% required)
})

# -- Column name flexibility ---------------------------------------------------

test_that("step_sec_conventional_cal accepts various column names for standards", {
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()

  # Test with 'time' column
  standards1 <- data.frame(
    time = c(11, 12.5, 14, 15.5, 17, 18.5),
    mw = c(1e6, 3e5, 1e5, 3e4, 1e4, 3e3)
  )

  rec1 <- recipes::recipe(~., data = test_data) |>
    step_sec_conventional_cal(standards = standards1, extrapolation = "none")

  expect_no_error(suppressWarnings(recipes::prep(rec1)))

  # Test with 'volume' column
  standards2 <- data.frame(
    volume = c(11, 12.5, 14, 15.5, 17, 18.5),
    molecular_weight = c(1e6, 3e5, 1e5, 3e4, 1e4, 3e3)
  )

  rec2 <- recipes::recipe(~., data = test_data) |>
    step_sec_conventional_cal(standards = standards2, extrapolation = "none")

  expect_no_error(suppressWarnings(recipes::prep(rec2)))
})

# -- Print and tidy methods ----------------------------------------------------

test_that("step_sec_conventional_cal print method works", {
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()
  standards <- create_ps_standards()

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_conventional_cal(standards = standards, extrapolation = "none")

  expect_output(print(rec), "conventional")
  expect_output(print(rec), "cubic")

  suppressWarnings(prepped <- recipes::prep(rec))
  expect_output(print(prepped), "R\u00b2")
  expect_output(print(prepped), "RMSE")
})

test_that("step_sec_conventional_cal tidy method returns expected structure", {
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()
  standards <- create_ps_standards()

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_conventional_cal(standards = standards, extrapolation = "none")

  suppressWarnings(prepped <- recipes::prep(rec))
  tidy_result <- recipes::tidy(prepped, number = 1)

  expect_s3_class(tidy_result, "tbl_df")

  # Core fit quality columns
  expect_true("fit_type" %in% names(tidy_result))
  expect_true("r_squared" %in% names(tidy_result))
  expect_true("adj_r_squared" %in% names(tidy_result))
  expect_true("rmse_log_mw" %in% names(tidy_result))
  expect_true("residual_std_error" %in% names(tidy_result))

  # Standard deviation metrics
  expect_true("max_abs_pct_deviation" %in% names(tidy_result))
  expect_true("mean_abs_pct_deviation" %in% names(tidy_result))

  # Calibration range and metadata
  expect_true("calibration_min" %in% names(tidy_result))
  expect_true("calibration_max" %in% names(tidy_result))
  expect_true("n_standards" %in% names(tidy_result))
  expect_true("degrees_of_freedom" %in% names(tidy_result))
  expect_true("coefficients" %in% names(tidy_result))
  expect_true("standard_results" %in% names(tidy_result))

  # Validate values
  expect_equal(tidy_result$fit_type, "cubic")
  expect_equal(tidy_result$n_standards, 6)
  expect_true(tidy_result$r_squared > 0.99)
  expect_true(tidy_result$rmse_log_mw >= 0)

  # Check standard_results is a nested tibble with expected columns
  std_results <- tidy_result$standard_results[[1]]
  expect_s3_class(std_results, "tbl_df")
  expect_equal(nrow(std_results), 6)
  expect_true(all(
    c(
      "location",
      "actual_log_mw",
      "predicted_log_mw",
      "residual_log_mw",
      "actual_mw",
      "predicted_mw",
      "pct_deviation",
      "prediction_se",
      "ci_lower_log_mw",
      "ci_upper_log_mw"
    ) %in%
      names(std_results)
  ))
})

# -- Integration with other steps ----------------------------------------------

test_that("step_sec_conventional_cal integrates with step_sec_mw_averages", {
  skip_if_not_installed("measure")

  # Create data that fits entirely within calibration range
  time <- seq(12, 17, by = 0.1) # Within 11-18.5 calibration range

  signal <- dnorm(time, mean = 14.5, sd = 1)

  test_data <- tibble::tibble(sample_id = "test")
  test_data$ri <- measure::new_measure_list(
    list(measure::new_measure_tbl(location = time, value = signal))
  )

  standards <- create_ps_standards()

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_conventional_cal(
      standards = standards,
      output_col = "log_mw",
      extrapolation = "none"
    ) |>
    step_sec_mw_averages(measures = "log_mw")

  suppressWarnings(prepped <- recipes::prep(rec))
  result <- recipes::bake(prepped, new_data = NULL)

  # Should have MW average columns
  expect_true("mw_mn" %in% names(result))
  expect_true("mw_mw" %in% names(result))
  expect_true("mw_dispersity" %in% names(result))
})
