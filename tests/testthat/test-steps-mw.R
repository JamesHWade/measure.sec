# ==============================================================================
# Tests for molecular weight calculation steps
# ==============================================================================

test_that("step_sec_mw_averages calculates Mn, Mw, Mz correctly", {
  skip_if_not_installed("measure")

  # Create test data with known MW distribution
  # Use a unimodal distribution for testing
  # Assume x-axis is log10(MW)
  log_mw <- seq(3, 6, by = 0.05)  # log10(MW) from 1000 to 1M
  n <- length(log_mw)

  # Concentration profile (Gaussian)
  conc <- dnorm(log_mw, mean = 4.5, sd = 0.5)
  conc <- conc / sum(conc)  # Normalize

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
  expect_equal(result$mw_dispersity[1], result$mw_mw[1] / result$mw_mn[1], tolerance = 0.01)
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
  frac_cols <- names(result)[grepl("frac", names(result))]
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
  baseline <- 0.1 + 0.005 * time  # Linear drift
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
