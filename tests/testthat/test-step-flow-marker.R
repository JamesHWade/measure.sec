# ==============================================================================
# Tests for step_sec_flow_marker
# ==============================================================================

# Helper to create test data with a flow marker peak
# Uses deterministic data (no random noise) for reproducible tests
create_flow_marker_data <- function(flow_marker_position = 19.0) {
  time <- seq(8, 22, by = 0.05)

  # Main polymer peak at ~14 mL
  polymer_peak <- 100 * dnorm(time, mean = 14, sd = 1.5)

  # Flow marker peak (sharp, small molecule) at specified position
  # Use larger amplitude and sharper peak for clear detection
  flow_marker <- 80 * dnorm(time, mean = flow_marker_position, sd = 0.1)

  # Deterministic signal - no noise for reproducible tests

  signal <- polymer_peak + flow_marker

  test_data <- tibble::tibble(sample_id = "test")
  test_data$ri <- measure::new_measure_list(
    list(measure::new_measure_tbl(location = time, value = signal))
  )
  test_data
}

# Helper to create multi-sample data with varying flow marker positions
# Uses deterministic data for reproducible tests
create_multi_sample_flow_data <- function() {
  time <- seq(8, 22, by = 0.05)

  # Sample 1: flow marker at 19.0 mL
  polymer1 <- 100 * dnorm(time, mean = 14, sd = 1.5)
  flow1 <- 80 * dnorm(time, mean = 19.0, sd = 0.1)
  signal1 <- polymer1 + flow1

  # Sample 2: flow marker shifted to 19.1 mL (simulating flow variation)
  polymer2 <- 100 * dnorm(time, mean = 14.1, sd = 1.5)
  flow2 <- 80 * dnorm(time, mean = 19.1, sd = 0.1)
  signal2 <- polymer2 + flow2

  # Sample 3: flow marker at 18.9 mL
  polymer3 <- 100 * dnorm(time, mean = 13.9, sd = 1.5)
  flow3 <- 80 * dnorm(time, mean = 18.9, sd = 0.1)
  signal3 <- polymer3 + flow3

  test_data <- tibble::tibble(sample_id = c("sample1", "sample2", "sample3"))
  test_data$ri <- measure::new_measure_list(
    list(
      measure::new_measure_tbl(location = time, value = signal1),
      measure::new_measure_tbl(location = time, value = signal2),
      measure::new_measure_tbl(location = time, value = signal3)
    )
  )
  test_data
}

# -- Constructor validation tests ----------------------------------------------

test_that("step_sec_flow_marker requires marker_range", {
  skip_if_not_installed("measure")

  test_data <- create_flow_marker_data()

  expect_error(
    recipes::recipe(~., data = test_data) |>
      step_sec_flow_marker(),
    "marker_range.*required"
  )
})

test_that("step_sec_flow_marker validates marker_range is numeric length 2", {
  skip_if_not_installed("measure")

  test_data <- create_flow_marker_data()

  expect_error(
    recipes::recipe(~., data = test_data) |>
      step_sec_flow_marker(marker_range = c(18)),
    "marker_range.*length 2"
  )

  expect_error(
    recipes::recipe(~., data = test_data) |>
      step_sec_flow_marker(marker_range = c("a", "b")),
    "marker_range.*numeric"
  )
})

test_that("step_sec_flow_marker validates marker_range order", {
  skip_if_not_installed("measure")

  test_data <- create_flow_marker_data()

  expect_error(
    recipes::recipe(~., data = test_data) |>
      step_sec_flow_marker(marker_range = c(20, 18)),
    "marker_range.*order"
  )
})

test_that("step_sec_flow_marker validates target_volume is single numeric", {
  skip_if_not_installed("measure")

  test_data <- create_flow_marker_data()

  expect_error(
    recipes::recipe(~., data = test_data) |>
      step_sec_flow_marker(marker_range = c(18, 20), target_volume = c(18, 19)),
    "target_volume.*single numeric"
  )
})

# -- Basic functionality tests -------------------------------------------------

test_that("step_sec_flow_marker can be prepped and baked", {
  skip_if_not_installed("measure")

  test_data <- create_flow_marker_data(flow_marker_position = 19.0)

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_flow_marker(
      marker_range = c(18, 20),
      target_volume = 19.0
    ) |>
    recipes::prep()

  result <- recipes::bake(rec, new_data = NULL)

  expect_s3_class(result, "tbl_df")
  expect_true("ri" %in% names(result))
  expect_true("flow_marker_correction" %in% names(result))
})

test_that("step_sec_flow_marker applies correction when flow marker is shifted", {
  skip_if_not_installed("measure")

  # Flow marker at 19.1 instead of target 19.0
  test_data <- create_flow_marker_data(flow_marker_position = 19.1)

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_flow_marker(
      marker_range = c(18, 20),
      target_volume = 19.0
    ) |>
    recipes::prep()

  result <- recipes::bake(rec, new_data = NULL)

  # Correction should be approximately 0.1 (allow some tolerance for peak detection)
  expect_equal(result$flow_marker_correction, 0.1, tolerance = 0.1)

  # Check that locations were shifted
  original_locs <- test_data$ri[[1]]$location
  corrected_locs <- result$ri[[1]]$location
  shift <- mean(original_locs - corrected_locs)
  expect_equal(shift, 0.1, tolerance = 0.1)
})

test_that("step_sec_flow_marker applies no correction when marker is at target", {
  skip_if_not_installed("measure")

  test_data <- create_flow_marker_data(flow_marker_position = 19.0)

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_flow_marker(
      marker_range = c(18, 20),
      target_volume = 19.0
    ) |>
    recipes::prep()

  result <- recipes::bake(rec, new_data = NULL)

  # Correction should be approximately 0
  expect_equal(result$flow_marker_correction, 0, tolerance = 0.05)
})

test_that("step_sec_flow_marker handles multiple samples", {
  skip_if_not_installed("measure")

  test_data <- create_multi_sample_flow_data()

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_flow_marker(
      marker_range = c(18, 20),
      target_volume = 19.0
    ) |>
    recipes::prep()

  result <- recipes::bake(rec, new_data = NULL)

  # Should have 3 correction values
  expect_length(result$flow_marker_correction, 3)

  # Corrections should reflect the shift from target
  # Sample 1: 19.0 -> 19.0 = ~0
  # Sample 2: 19.1 -> 19.0 = ~0.1
  # Sample 3: 18.9 -> 19.0 = ~-0.1
  expect_equal(result$flow_marker_correction[1], 0, tolerance = 0.05)
  expect_equal(result$flow_marker_correction[2], 0.1, tolerance = 0.05)
  expect_equal(result$flow_marker_correction[3], -0.1, tolerance = 0.05)
})

# -- Auto-detection tests ------------------------------------------------------

test_that("step_sec_flow_marker auto-detects flow marker with second derivative", {
  skip_if_not_installed("measure")

  test_data <- create_flow_marker_data(flow_marker_position = 19.2)

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_flow_marker(
      marker_range = c(18, 20),
      target_volume = 19.0,
      auto_detect = TRUE
    ) |>
    recipes::prep()

  result <- recipes::bake(rec, new_data = NULL)

  # Should detect the sharp flow marker peak (allow for peak detection variance)
  expect_equal(result$flow_marker_correction, 0.2, tolerance = 0.15)
})

test_that("step_sec_flow_marker works with auto_detect = FALSE", {
  skip_if_not_installed("measure")

  test_data <- create_flow_marker_data(flow_marker_position = 19.15)

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_flow_marker(
      marker_range = c(18, 20),
      target_volume = 19.0,
      auto_detect = FALSE
    ) |>
    recipes::prep()

  result <- recipes::bake(rec, new_data = NULL)

  # Should still find the peak maximum
  expect_equal(result$flow_marker_correction, 0.15, tolerance = 0.1)
})

# -- store_correction option ---------------------------------------------------

test_that("step_sec_flow_marker respects store_correction = FALSE", {
  skip_if_not_installed("measure")

  test_data <- create_flow_marker_data()

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_flow_marker(
      marker_range = c(18, 20),
      target_volume = 19.0,
      store_correction = FALSE
    ) |>
    recipes::prep()

  result <- recipes::bake(rec, new_data = NULL)

  expect_false("flow_marker_correction" %in% names(result))
})

# -- tidy method ---------------------------------------------------------------

test_that("tidy method returns expected columns", {
  skip_if_not_installed("measure")

  test_data <- create_flow_marker_data()

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_flow_marker(
      marker_range = c(18, 20),
      target_volume = 19.0
    ) |>
    recipes::prep()

  tidy_result <- recipes::tidy(rec, number = 1)

  expect_true("marker_range_min" %in% names(tidy_result))
  expect_true("marker_range_max" %in% names(tidy_result))
  expect_true("target_volume" %in% names(tidy_result))
  expect_true("auto_detect" %in% names(tidy_result))
  expect_equal(tidy_result$marker_range_min, 18)
  expect_equal(tidy_result$marker_range_max, 20)
  expect_equal(tidy_result$target_volume, 19.0)
})

# -- print method --------------------------------------------------------------

test_that("print method works for trained and untrained steps", {
  skip_if_not_installed("measure")

  test_data <- create_flow_marker_data()

  # Untrained
  rec_untrained <- recipes::recipe(~., data = test_data) |>
    step_sec_flow_marker(marker_range = c(18, 20))

  expect_output(print(rec_untrained$steps[[1]]), "flow marker")

  # Trained
  rec_trained <- recipes::prep(rec_untrained)
  expect_output(print(rec_trained$steps[[1]]), "flow marker")
})
