# ==============================================================================
# Tests for step_sec_integration_window
# ==============================================================================

# Helper to create test data with measure columns
create_test_sec_data <- function() {
  # Create SEC-like chromatogram
  time <- seq(8, 20, by = 0.05)
  # Single Gaussian peak centered at 14
  signal <- 100 *
    dnorm(time, mean = 14, sd = 1.5) +
    rnorm(length(time), sd = 0.5)

  test_data <- tibble::tibble(sample_id = "test")
  test_data$ri <- measure::new_measure_list(
    list(measure::new_measure_tbl(location = time, value = signal))
  )
  test_data
}

# Helper to create multi-sample data
create_multi_sample_data <- function() {
  time <- seq(8, 20, by = 0.05)

  # Two samples with different peak positions
  signal1 <- 100 *
    dnorm(time, mean = 13, sd = 1.2) +
    rnorm(length(time), sd = 0.3)
  signal2 <- 100 *
    dnorm(time, mean = 15, sd = 1.5) +
    rnorm(length(time), sd = 0.3)

  test_data <- tibble::tibble(sample_id = c("sample1", "sample2"))
  test_data$ri <- measure::new_measure_list(
    list(
      measure::new_measure_tbl(location = time, value = signal1),
      measure::new_measure_tbl(location = time, value = signal2)
    )
  )
  test_data
}

# -- Constructor validation tests ----------------------------------------------

test_that("step_sec_integration_window validates start parameter", {
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()

  expect_error(
    recipes::recipe(~., data = test_data) |>
      step_sec_integration_window(start = "abc"),
    "start.*single numeric"
  )

  expect_error(
    recipes::recipe(~., data = test_data) |>
      step_sec_integration_window(start = c(1, 2)),
    "start.*single numeric"
  )
})

test_that("step_sec_integration_window validates end parameter", {
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()

  expect_error(
    recipes::recipe(~., data = test_data) |>
      step_sec_integration_window(end = "abc"),
    "end.*single numeric"
  )
})

test_that("step_sec_integration_window validates start < end", {
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()

  expect_error(
    recipes::recipe(~., data = test_data) |>
      step_sec_integration_window(start = 15, end = 10),
    "start.*less than.*end"
  )

  expect_error(
    recipes::recipe(~., data = test_data) |>
      step_sec_integration_window(start = 10, end = 10),
    "start.*less than.*end"
  )
})

test_that("step_sec_integration_window validates signal_threshold", {
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()

  expect_error(
    recipes::recipe(~., data = test_data) |>
      step_sec_integration_window(signal_threshold = 0),
    "signal_threshold.*between 0 and 1"
  )

  expect_error(
    recipes::recipe(~., data = test_data) |>
      step_sec_integration_window(signal_threshold = 1),
    "signal_threshold.*between 0 and 1"
  )

  expect_error(
    recipes::recipe(~., data = test_data) |>
      step_sec_integration_window(signal_threshold = -0.5),
    "signal_threshold.*between 0 and 1"
  )
})

test_that("step_sec_integration_window validates calibration_range", {
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()

  expect_error(
    recipes::recipe(~., data = test_data) |>
      step_sec_integration_window(calibration_range = c(15, 10)),
    "calibration_range.*first element.*less than"
  )

  expect_error(
    recipes::recipe(~., data = test_data) |>
      step_sec_integration_window(calibration_range = c(10)),
    "calibration_range.*length 2"
  )
})

test_that("step_sec_integration_window validates min_window_width", {
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()

  expect_error(
    recipes::recipe(~., data = test_data) |>
      step_sec_integration_window(min_window_width = 0),
    "min_window_width.*positive"
  )

  expect_error(
    recipes::recipe(~., data = test_data) |>
      step_sec_integration_window(min_window_width = -1),
    "min_window_width.*positive"
  )
})

# -- Prep and bake tests -------------------------------------------------------

test_that("step_sec_integration_window can prep with explicit bounds", {
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_integration_window(start = 10, end = 18)

  prepped <- recipes::prep(rec)
  expect_s3_class(prepped, "recipe")
})

test_that("step_sec_integration_window can prep with auto-detection", {
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_integration_window(auto_detect = TRUE)

  prepped <- recipes::prep(rec)
  expect_s3_class(prepped, "recipe")
})

test_that("step_sec_integration_window bake creates .integration_window column", {
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_integration_window(start = 10, end = 18)

  result <- recipes::prep(rec) |>
    recipes::bake(new_data = NULL)

  expect_true(".integration_window" %in% names(result))
  expect_equal(length(result$.integration_window), 1)
})

test_that("step_sec_integration_window output has correct structure", {
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_integration_window(start = 10, end = 18)

  result <- recipes::prep(rec) |>
    recipes::bake(new_data = NULL)

  window <- result$.integration_window[[1]]
  expect_s3_class(window, "tbl_df")
  expect_true("start" %in% names(window))
  expect_true("end" %in% names(window))
  expect_equal(window$start, 10)
  expect_equal(window$end, 18)
})

test_that("step_sec_integration_window works with multiple samples", {
  skip_if_not_installed("measure")

  test_data <- create_multi_sample_data()

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_integration_window(start = 10, end = 18)

  result <- recipes::prep(rec) |>
    recipes::bake(new_data = NULL)

  expect_equal(nrow(result), 2)
  expect_equal(length(result$.integration_window), 2)

  # Both samples should have the same window
  expect_equal(result$.integration_window[[1]]$start, 10)
  expect_equal(result$.integration_window[[2]]$start, 10)
  expect_equal(result$.integration_window[[1]]$end, 18)
  expect_equal(result$.integration_window[[2]]$end, 18)
})

test_that("step_sec_integration_window auto-detect finds signal region", {
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_integration_window(auto_detect = TRUE, signal_threshold = 0.01)

  result <- recipes::prep(rec) |>
    recipes::bake(new_data = NULL)

  window <- result$.integration_window[[1]]

  # Peak is centered at 14 with sd=1.5, should detect roughly 10-18
  expect_lt(window$start, 14)
  expect_gt(window$end, 14)
  expect_gt(window$end - window$start, 2) # Should be wider than 2 units
})

test_that("step_sec_integration_window respects calibration_range", {
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_integration_window(
      auto_detect = TRUE,
      calibration_range = c(10, 16),
      extend_beyond_cal = 0.5
    )

  result <- recipes::prep(rec) |>
    recipes::bake(new_data = NULL)

  window <- result$.integration_window[[1]]

  # Start should not go below calibration minimum
  expect_gte(window$start, 10)

  # End can extend up to 50% beyond calibration max (16 + 0.5*6 = 19)
  expect_lte(window$end, 19)
})

test_that("step_sec_integration_window enforces min_window_width", {
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()

  # Use narrow explicit bounds that should be expanded
  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_integration_window(
      start = 14,
      end = 14.5,
      min_window_width = 2.0
    )

  result <- suppressMessages(
    recipes::prep(rec) |>
      recipes::bake(new_data = NULL)
  )

  window <- result$.integration_window[[1]]

  # Window should be at least 2.0 wide
  expect_gte(window$end - window$start, 2.0)
})

test_that("step_sec_integration_window errors without auto_detect and missing bounds", {
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_integration_window(auto_detect = FALSE, start = 10)

  expect_error(
    recipes::prep(rec),
    "Cannot determine integration window"
  )
})

# -- Tidy method tests ---------------------------------------------------------

test_that("tidy.step_sec_integration_window returns expected tibble", {
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_integration_window(start = 10, end = 18)

  # Before prep
  tidy_result <- recipes::tidy(rec, number = 1)
  expect_s3_class(tidy_result, "tbl_df")
  expect_true("start" %in% names(tidy_result))
  expect_true("end" %in% names(tidy_result))
  expect_true("auto_detect" %in% names(tidy_result))
  expect_equal(tidy_result$start, 10)
  expect_equal(tidy_result$end, 18)

  # After prep
  prepped <- recipes::prep(rec)
  tidy_prepped <- recipes::tidy(prepped, number = 1)
  expect_s3_class(tidy_prepped, "tbl_df")
  expect_equal(tidy_prepped$start, 10)
  expect_equal(tidy_prepped$end, 18)
})

test_that("tidy shows computed values after auto-detect prep", {
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_integration_window(auto_detect = TRUE)

  # Before prep - start and end should be NA
  tidy_result <- recipes::tidy(rec, number = 1)
  expect_true(is.na(tidy_result$start))
  expect_true(is.na(tidy_result$end))

  # After prep - should have computed values
  prepped <- recipes::prep(rec)
  tidy_prepped <- recipes::tidy(prepped, number = 1)
  expect_false(is.na(tidy_prepped$start))
  expect_false(is.na(tidy_prepped$end))
})

# -- Print method tests --------------------------------------------------------

test_that("print.step_sec_integration_window works", {
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()

  # Explicit bounds
  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_integration_window(start = 10, end = 18)

  expect_output(print(rec), "10.*18")

  # Auto-detect
  rec_auto <- recipes::recipe(~., data = test_data) |>
    step_sec_integration_window(auto_detect = TRUE)

  expect_output(print(rec_auto), "auto-detect")
})

# -- required_pkgs tests -------------------------------------------------------

test_that("required_pkgs.step_sec_integration_window lists measure packages", {
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_integration_window(start = 10, end = 18)

  pkgs <- recipes::required_pkgs(rec)
  expect_true("measure.sec" %in% pkgs)
  expect_true("measure" %in% pkgs)
})
