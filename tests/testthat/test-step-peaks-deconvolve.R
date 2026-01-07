# ==============================================================================
# Tests for step_sec_peaks_deconvolve
# ==============================================================================

# Helper to create test data with measure and peaks columns
create_test_peaks_data <- function() {
  time <- seq(8, 20, by = 0.05)

  # Two overlapping peaks
  signal <- 100 *
    dnorm(time, mean = 12, sd = 1.0) +
    60 * dnorm(time, mean = 14, sd = 1.2) +
    rnorm(length(time), sd = 0.3)

  test_data <- tibble::tibble(sample_id = "test")
  test_data$.measures <- measure::new_measure_list(
    list(measure::new_measure_tbl(location = time, value = signal))
  )

  # Add detected peaks (simulated)
  peaks <- measure:::new_peaks_tbl(
    peak_id = c(1L, 2L),
    location = c(12.0, 14.0),
    height = c(40.0, 24.0),
    left_base = c(10.0, 12.5),
    right_base = c(13.5, 16.0),
    area = c(100.0, 72.0)
  )
  test_data$.peaks <- measure:::new_peaks_list(list(peaks))

  test_data
}

# Helper for multi-sample data
create_multi_sample_peaks_data <- function() {
  time <- seq(8, 20, by = 0.05)

  signal1 <- 100 *
    dnorm(time, mean = 12, sd = 1.0) +
    60 * dnorm(time, mean = 14, sd = 1.2) +
    rnorm(length(time), sd = 0.3)
  signal2 <- 80 *
    dnorm(time, mean = 13, sd = 1.1) +
    rnorm(length(time), sd = 0.3)

  test_data <- tibble::tibble(sample_id = c("sample1", "sample2"))
  test_data$.measures <- measure::new_measure_list(
    list(
      measure::new_measure_tbl(location = time, value = signal1),
      measure::new_measure_tbl(location = time, value = signal2)
    )
  )

  # Peaks for each sample
  peaks1 <- measure:::new_peaks_tbl(
    peak_id = c(1L, 2L),
    location = c(12.0, 14.0),
    height = c(40.0, 24.0),
    left_base = c(10.0, 12.5),
    right_base = c(13.5, 16.0),
    area = c(100.0, 72.0)
  )
  peaks2 <- measure:::new_peaks_tbl(
    peak_id = 1L,
    location = 13.0,
    height = 32.0,
    left_base = 10.5,
    right_base = 15.5,
    area = 88.0
  )
  test_data$.peaks <- measure:::new_peaks_list(list(peaks1, peaks2))

  test_data
}

# -- Constructor validation tests ----------------------------------------------

test_that("step_sec_peaks_deconvolve validates model parameter", {
  skip_if_not_installed("measure")

  test_data <- create_test_peaks_data()

  expect_error(
    recipes::recipe(~., data = test_data) |>
      step_sec_peaks_deconvolve(model = "invalid"),
    "model"
  )
})

test_that("step_sec_peaks_deconvolve validates optimizer parameter", {
  skip_if_not_installed("measure")

  test_data <- create_test_peaks_data()

  expect_error(
    recipes::recipe(~., data = test_data) |>
      step_sec_peaks_deconvolve(optimizer = "invalid"),
    "optimizer"
  )
})

test_that("step_sec_peaks_deconvolve validates max_iter parameter", {
  skip_if_not_installed("measure")

  test_data <- create_test_peaks_data()

  expect_error(
    recipes::recipe(~., data = test_data) |>
      step_sec_peaks_deconvolve(max_iter = -10),
    "max_iter"
  )

  expect_error(
    recipes::recipe(~., data = test_data) |>
      step_sec_peaks_deconvolve(max_iter = "100"),
    "max_iter"
  )
})

test_that("step_sec_peaks_deconvolve validates quality_threshold parameter", {
  skip_if_not_installed("measure")

  test_data <- create_test_peaks_data()

  expect_error(
    recipes::recipe(~., data = test_data) |>
      step_sec_peaks_deconvolve(quality_threshold = 1.5),
    "quality_threshold"
  )

  expect_error(
    recipes::recipe(~., data = test_data) |>
      step_sec_peaks_deconvolve(quality_threshold = -0.1),
    "quality_threshold"
  )
})

test_that("step_sec_peaks_deconvolve validates smart_init parameter", {
  skip_if_not_installed("measure")

  test_data <- create_test_peaks_data()

  expect_error(
    recipes::recipe(~., data = test_data) |>
      step_sec_peaks_deconvolve(smart_init = "yes"),
    "smart_init"
  )
})

test_that("step_sec_peaks_deconvolve validates constrain_positions parameter", {
  skip_if_not_installed("measure")

  test_data <- create_test_peaks_data()

  expect_error(
    recipes::recipe(~., data = test_data) |>
      step_sec_peaks_deconvolve(constrain_positions = "yes"),
    "constrain_positions"
  )
})

# -- Default value tests -------------------------------------------------------

test_that("step_sec_peaks_deconvolve has SEC-appropriate defaults", {
  skip_if_not_installed("measure")

  test_data <- create_test_peaks_data()

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_peaks_deconvolve()

  # Check defaults via tidy
  tidy_result <- recipes::tidy(rec, number = 1)

  # EMG is the default model for SEC (handles chromatographic tailing)
  expect_equal(tidy_result$model, "emg")

  # Higher quality threshold for SEC (0.9 vs measure's 0.8)
  expect_equal(tidy_result$quality_threshold, 0.9)

  # Smart init should be enabled by default
  expect_equal(tidy_result$smart_init, TRUE)
})

# -- Prep tests ----------------------------------------------------------------

test_that("step_sec_peaks_deconvolve requires peaks column", {
  skip_if_not_installed("measure")

  # Data without peaks
  time <- seq(8, 20, by = 0.05)
  signal <- 100 * dnorm(time, mean = 14, sd = 1.5)
  test_data <- tibble::tibble(sample_id = "test")
  test_data$.measures <- measure::new_measure_list(
    list(measure::new_measure_tbl(location = time, value = signal))
  )

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_peaks_deconvolve()

  expect_error(
    recipes::prep(rec),
    "peaks"
  )
})

test_that("step_sec_peaks_deconvolve can prep with valid data", {
  skip_if_not_installed("measure")

  test_data <- create_test_peaks_data()

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_peaks_deconvolve()

  prepped <- recipes::prep(rec)
  expect_s3_class(prepped, "recipe")
})

test_that("step_sec_peaks_deconvolve auto-detects measure column", {
  skip_if_not_installed("measure")

  # Create data with a named measure column instead of .measures
  time <- seq(8, 20, by = 0.05)
  signal <- 100 *
    dnorm(time, mean = 12, sd = 1.0) +
    60 * dnorm(time, mean = 14, sd = 1.2)

  test_data <- tibble::tibble(sample_id = "test")
  test_data$ri <- measure::new_measure_list(
    list(measure::new_measure_tbl(location = time, value = signal))
  )

  peaks <- measure:::new_peaks_tbl(
    peak_id = c(1L, 2L),
    location = c(12.0, 14.0),
    height = c(40.0, 24.0),
    left_base = c(10.0, 12.5),
    right_base = c(13.5, 16.0),
    area = c(100.0, 72.0)
  )
  test_data$.peaks <- measure:::new_peaks_list(list(peaks))

  # Should auto-detect 'ri' column
  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_peaks_deconvolve()

  # Should prep without error (may give informative message)
  expect_message(
    prepped <- recipes::prep(rec),
    regexp = "ri|measure",
    ignore.case = TRUE
  )
})

# -- Model options tests -------------------------------------------------------

test_that("step_sec_peaks_deconvolve accepts gaussian model", {
  skip_if_not_installed("measure")

  test_data <- create_test_peaks_data()

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_peaks_deconvolve(model = "gaussian")

  prepped <- recipes::prep(rec)
  expect_s3_class(prepped, "recipe")

  tidy_result <- recipes::tidy(prepped, number = 1)
  expect_equal(tidy_result$model, "gaussian")
})

test_that("step_sec_peaks_deconvolve accepts bigaussian model", {
  skip_if_not_installed("measure")

  test_data <- create_test_peaks_data()

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_peaks_deconvolve(model = "bigaussian")

  prepped <- recipes::prep(rec)
  expect_s3_class(prepped, "recipe")

  tidy_result <- recipes::tidy(prepped, number = 1)
  expect_equal(tidy_result$model, "bigaussian")
})

# -- Optimizer options tests ---------------------------------------------------

test_that("step_sec_peaks_deconvolve accepts all optimizer options", {
  skip_if_not_installed("measure")

  test_data <- create_test_peaks_data()

  for (opt in c("auto", "lbfgsb", "multistart", "nelder_mead")) {
    rec <- recipes::recipe(~., data = test_data) |>
      step_sec_peaks_deconvolve(optimizer = opt)

    prepped <- recipes::prep(rec)
    expect_s3_class(prepped, "recipe")

    tidy_result <- recipes::tidy(prepped, number = 1)
    expect_equal(tidy_result$optimizer, opt)
  }
})

# -- Tidy method tests ---------------------------------------------------------

test_that("tidy.step_sec_peaks_deconvolve returns expected tibble", {
  skip_if_not_installed("measure")

  test_data <- create_test_peaks_data()

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_peaks_deconvolve(
      model = "emg",
      optimizer = "multistart",
      max_iter = 1000L,
      quality_threshold = 0.95,
      smart_init = FALSE
    )

  # Before prep
  tidy_result <- recipes::tidy(rec, number = 1)
  expect_s3_class(tidy_result, "tbl_df")
  expect_true(all(
    c(
      "model",
      "optimizer",
      "max_iter",
      "quality_threshold",
      "smart_init",
      "id"
    ) %in%
      names(tidy_result)
  ))
  expect_equal(tidy_result$model, "emg")
  expect_equal(tidy_result$optimizer, "multistart")
  expect_equal(tidy_result$max_iter, 1000L)
  expect_equal(tidy_result$quality_threshold, 0.95)
  expect_equal(tidy_result$smart_init, FALSE)

  # After prep
  prepped <- recipes::prep(rec)
  tidy_prepped <- recipes::tidy(prepped, number = 1)
  expect_s3_class(tidy_prepped, "tbl_df")
  expect_equal(tidy_prepped$model, "emg")
  expect_equal(tidy_prepped$optimizer, "multistart")
})

# -- Print method tests --------------------------------------------------------

test_that("print.step_sec_peaks_deconvolve works", {
  skip_if_not_installed("measure")

  test_data <- create_test_peaks_data()

  # With EMG model
  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_peaks_deconvolve(model = "emg", optimizer = "multistart")

  expect_output(print(rec), "emg")
  expect_output(print(rec), "multistart")

  # With Gaussian model
  rec2 <- recipes::recipe(~., data = test_data) |>
    step_sec_peaks_deconvolve(model = "gaussian", optimizer = "auto")

  expect_output(print(rec2), "gaussian")
  expect_output(print(rec2), "auto")
})

# -- required_pkgs tests -------------------------------------------------------

test_that("required_pkgs.step_sec_peaks_deconvolve lists required packages", {
  skip_if_not_installed("measure")

  test_data <- create_test_peaks_data()

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_peaks_deconvolve()

  pkgs <- recipes::required_pkgs(rec)
  expect_true("measure.sec" %in% pkgs)
  expect_true("measure" %in% pkgs)
})

# -- Multi-sample tests --------------------------------------------------------

test_that("step_sec_peaks_deconvolve works with multiple samples", {
  skip_if_not_installed("measure")

  test_data <- create_multi_sample_peaks_data()

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_peaks_deconvolve()

  prepped <- recipes::prep(rec)
  expect_s3_class(prepped, "recipe")
})
