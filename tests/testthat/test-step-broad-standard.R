# ==============================================================================
# Tests for step_sec_broad_standard
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

# Helper to create a broad standard chromatogram
# This simulates a polydisperse polymer with known Mn and Mw
create_broad_standard <- function(mn = 50000, mw = 150000) {
  # Generate a broad distribution
  # Use a log-normal distribution centered appropriately
  dispersity <- mw / mn

  # Simulate chromatogram spanning typical SEC range
  time <- seq(10, 20, by = 0.05)

  # Peak position should be related to Mw
  # For typical SEC: higher MW elutes earlier
  # Center around 14-16 min range
  peak_time <- 15

  # Width relates to dispersity - broader dispersity = wider peak
  sd_time <- 0.5 + 0.5 * (dispersity - 1)

  # Create asymmetric peak (typical for broad polymers)
  signal <- dnorm(time, mean = peak_time, sd = sd_time)

  # Add slight tail
  signal <- signal + 0.3 * dnorm(time, mean = peak_time + 1, sd = sd_time * 1.2)

  # Normalize

  signal <- signal / max(signal)

  data.frame(
    time = time,
    signal = signal
  )
}

# -- Constructor validation tests ----------------------------------------------

test_that("step_sec_broad_standard requires broad_standard", {
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()

  expect_error(
    recipes::recipe(~., data = test_data) |>
      step_sec_broad_standard(known_mn = 50000, known_mw = 150000),
    "broad_standard.*required"
  )
})

test_that("step_sec_broad_standard requires known_mn and known_mw", {
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()
  broad_std <- create_broad_standard()

  expect_error(
    recipes::recipe(~., data = test_data) |>
      step_sec_broad_standard(broad_standard = broad_std, known_mw = 150000),
    "known_mn.*known_mw.*required"
  )

  expect_error(
    recipes::recipe(~., data = test_data) |>
      step_sec_broad_standard(broad_standard = broad_std, known_mn = 50000),
    "known_mn.*known_mw.*required"
  )
})

test_that("step_sec_broad_standard validates broad_standard is a data frame", {
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()

  expect_error(
    recipes::recipe(~., data = test_data) |>
      step_sec_broad_standard(
        broad_standard = c(1, 2, 3),
        known_mn = 50000,
        known_mw = 150000
      ),
    "must be a data frame"
  )
})

test_that("step_sec_broad_standard validates broad_standard has location column", {
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()

  bad_standard <- data.frame(
    x = seq(10, 20, by = 0.1),
    signal = dnorm(seq(10, 20, by = 0.1), mean = 15, sd = 1.5)
  )

  expect_error(
    recipes::recipe(~., data = test_data) |>
      step_sec_broad_standard(
        broad_standard = bad_standard,
        known_mn = 50000,
        known_mw = 150000
      ),
    "location column"
  )
})

test_that("step_sec_broad_standard validates broad_standard has signal column", {
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()

  bad_standard <- data.frame(
    time = seq(10, 20, by = 0.1),
    y = dnorm(seq(10, 20, by = 0.1), mean = 15, sd = 1.5)
  )

  expect_error(
    recipes::recipe(~., data = test_data) |>
      step_sec_broad_standard(
        broad_standard = bad_standard,
        known_mn = 50000,
        known_mw = 150000
      ),
    "signal column"
  )
})

test_that("step_sec_broad_standard validates known_mn is positive", {
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()
  broad_std <- create_broad_standard()

  expect_error(
    recipes::recipe(~., data = test_data) |>
      step_sec_broad_standard(
        broad_standard = broad_std,
        known_mn = -50000,
        known_mw = 150000
      ),
    "positive number"
  )
})

test_that("step_sec_broad_standard validates known_mw >= known_mn", {
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()
  broad_std <- create_broad_standard()

  expect_error(
    recipes::recipe(~., data = test_data) |>
      step_sec_broad_standard(
        broad_standard = broad_std,
        known_mn = 150000,
        known_mw = 50000
      ),
    "greater than or equal"
  )
})

# -- Prep tests ----------------------------------------------------------------

test_that("step_sec_broad_standard fits calibration during prep", {
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()
  broad_std <- create_broad_standard(mn = 50000, mw = 150000)

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_broad_standard(
      broad_standard = broad_std,
      known_mn = 50000,
      known_mw = 150000,
      fit_type = "linear",
      extrapolation = "none"
    )

  suppressWarnings(prepped <- recipes::prep(rec))

  step <- prepped$steps[[1]]
  expect_true(step$trained)
  expect_true(!is.null(step$calibration_coefficients))
  expect_true(!is.null(step$calibration_range))
  expect_true(!is.null(step$calibration_diagnostics))

  # Check coefficient names
  coefs <- step$calibration_coefficients
  expect_true("intercept" %in% names(coefs))
  expect_true("slope" %in% names(coefs))

  # Slope should be negative (higher MW elutes first)
  expect_true(coefs["slope"] < 0)
})

test_that("step_sec_broad_standard produces reasonable diagnostics", {
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()
  broad_std <- create_broad_standard(mn = 50000, mw = 150000)

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_broad_standard(
      broad_standard = broad_std,
      known_mn = 50000,
      known_mw = 150000,
      extrapolation = "none"
    )

  suppressWarnings(prepped <- recipes::prep(rec))

  diag <- prepped$steps[[1]]$calibration_diagnostics

  # Check that calculated values are close to known values
  # Allow some tolerance since this is an optimization
  expect_true(abs(diag$mn_error_pct) < 10) # Within 10%
  expect_true(abs(diag$mw_error_pct) < 10)
  expect_true(abs(diag$dispersity_error_pct) < 10)

  # Check diagnostics structure
  expect_true(!is.null(diag$known_mn))
  expect_true(!is.null(diag$known_mw))
  expect_true(!is.null(diag$known_dispersity))
  expect_true(!is.null(diag$calculated_mn))
  expect_true(!is.null(diag$calculated_mw))
  expect_true(!is.null(diag$calculated_dispersity))
  expect_equal(diag$known_dispersity, 150000 / 50000) # 3.0
})

test_that("step_sec_broad_standard warns on poor fit", {
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()

  # Create a very narrow standard that won't fit well with given Mn/Mw
  narrow_std <- data.frame(
    time = seq(14.5, 15.5, by = 0.01),
    signal = dnorm(seq(14.5, 15.5, by = 0.01), mean = 15, sd = 0.1)
  )

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_broad_standard(
      broad_standard = narrow_std,
      known_mn = 10000,
      known_mw = 500000, # Very high dispersity
      extrapolation = "none"
    )

  # Should warn about poor fit
  expect_warning(
    recipes::prep(rec),
    ">5% error"
  )
})

# -- Bake tests ----------------------------------------------------------------

test_that("step_sec_broad_standard creates MW column on bake", {
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()
  broad_std <- create_broad_standard(mn = 50000, mw = 150000)

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_broad_standard(
      broad_standard = broad_std,
      known_mn = 50000,
      known_mw = 150000,
      output_col = "calibrated_mw",
      extrapolation = "none"
    )

  suppressWarnings(prepped <- recipes::prep(rec))
  result <- recipes::bake(prepped, new_data = NULL)

  expect_true("calibrated_mw" %in% names(result))
  expect_s3_class(result$calibrated_mw, "measure_list")

  # Check values are reasonable (log MW between 3 and 8) for non-NA values
  mw_values <- result$calibrated_mw[[1]]$value
  valid_values <- mw_values[!is.na(mw_values)]
  expect_true(length(valid_values) > 0)
  expect_true(all(valid_values >= 2 & valid_values <= 8))
})

test_that("step_sec_broad_standard returns MW in Daltons when log_output = FALSE", {
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()
  broad_std <- create_broad_standard(mn = 50000, mw = 150000)

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_broad_standard(
      broad_standard = broad_std,
      known_mn = 50000,
      known_mw = 150000,
      log_output = FALSE,
      extrapolation = "none"
    )

  suppressWarnings(prepped <- recipes::prep(rec))
  result <- recipes::bake(prepped, new_data = NULL)

  # Values should be in Daltons (typically > 100)
  mw_values <- result$mw[[1]]$value
  valid_values <- mw_values[!is.na(mw_values)]
  expect_true(length(valid_values) > 0)
  expect_true(all(valid_values > 100))
})

test_that("step_sec_broad_standard warns for out-of-range data", {
  skip_if_not_installed("measure")

  # Create data that extends beyond calibration range
  time <- seq(5, 25, by = 0.1)
  signal <- dnorm(time, mean = 15, sd = 3)

  test_data <- tibble::tibble(sample_id = "test")
  test_data$ri <- measure::new_measure_list(
    list(measure::new_measure_tbl(location = time, value = signal))
  )

  broad_std <- create_broad_standard(mn = 50000, mw = 150000)

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_broad_standard(
      broad_standard = broad_std,
      known_mn = 50000,
      known_mw = 150000,
      extrapolation = "warn"
    )

  # Warning is emitted during prep (when it internally bakes training data)
  expect_warning(
    recipes::prep(rec),
    "outside calibration range"
  )
})

test_that("step_sec_broad_standard returns NA for out-of-range with extrapolation='none'", {
  skip_if_not_installed("measure")

  # Create data that extends beyond typical range
  time <- seq(5, 25, by = 0.1)
  signal <- dnorm(time, mean = 15, sd = 3)

  test_data <- tibble::tibble(sample_id = "test")
  test_data$ri <- measure::new_measure_list(
    list(measure::new_measure_tbl(location = time, value = signal))
  )

  broad_std <- create_broad_standard(mn = 50000, mw = 150000)

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_broad_standard(
      broad_standard = broad_std,
      known_mn = 50000,
      known_mw = 150000,
      extrapolation = "none"
    )

  suppressWarnings(prepped <- recipes::prep(rec))
  result <- recipes::bake(prepped, new_data = NULL)

  mw_values <- result$mw[[1]]$value
  cal_range <- prepped$steps[[1]]$calibration_range

  # Values outside calibration range should be NA
  expect_true(all(is.na(mw_values[time < cal_range[1]])))
  expect_true(all(is.na(mw_values[time > cal_range[2]])))
})

# -- Different fit types -------------------------------------------------------

test_that("step_sec_broad_standard works with linear fit", {
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()
  broad_std <- create_broad_standard(mn = 50000, mw = 150000)

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_broad_standard(
      broad_standard = broad_std,
      known_mn = 50000,
      known_mw = 150000,
      fit_type = "linear",
      extrapolation = "none"
    )

  suppressWarnings(prepped <- recipes::prep(rec))
  result <- recipes::bake(prepped, new_data = NULL)

  expect_true("mw" %in% names(result))

  # Check coefficients
  coefs <- prepped$steps[[1]]$calibration_coefficients
  expect_equal(length(coefs), 2)
  expect_true(all(c("intercept", "slope") %in% names(coefs)))
})

test_that("step_sec_broad_standard works with quadratic fit", {
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()
  broad_std <- create_broad_standard(mn = 50000, mw = 150000)

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_broad_standard(
      broad_standard = broad_std,
      known_mn = 50000,
      known_mw = 150000,
      fit_type = "quadratic",
      extrapolation = "none"
    )

  suppressWarnings(prepped <- recipes::prep(rec))
  result <- recipes::bake(prepped, new_data = NULL)

  expect_true("mw" %in% names(result))

  # Check coefficients
  coefs <- prepped$steps[[1]]$calibration_coefficients
  expect_equal(length(coefs), 3)
  expect_true(all(c("intercept", "slope", "quadratic") %in% names(coefs)))
})

# -- Column name flexibility ---------------------------------------------------

test_that("step_sec_broad_standard accepts various column names", {
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()

  # Test with 'volume' and 'response' columns
  broad_std1 <- data.frame(
    volume = seq(10, 20, by = 0.1),
    response = dnorm(seq(10, 20, by = 0.1), mean = 15, sd = 1.5)
  )

  rec1 <- recipes::recipe(~., data = test_data) |>
    step_sec_broad_standard(
      broad_standard = broad_std1,
      known_mn = 50000,
      known_mw = 150000,
      extrapolation = "none"
    )

  expect_no_error(suppressWarnings(recipes::prep(rec1)))

  # Test with 'retention' and 'ri' columns
  broad_std2 <- data.frame(
    retention = seq(10, 20, by = 0.1),
    ri = dnorm(seq(10, 20, by = 0.1), mean = 15, sd = 1.5)
  )

  rec2 <- recipes::recipe(~., data = test_data) |>
    step_sec_broad_standard(
      broad_standard = broad_std2,
      known_mn = 50000,
      known_mw = 150000,
      extrapolation = "none"
    )

  expect_no_error(suppressWarnings(recipes::prep(rec2)))
})

# -- Print and tidy methods ----------------------------------------------------

test_that("step_sec_broad_standard print method works", {
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()
  broad_std <- create_broad_standard(mn = 50000, mw = 150000)

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_broad_standard(
      broad_standard = broad_std,
      known_mn = 50000,
      known_mw = 150000,
      extrapolation = "none"
    )

  expect_output(print(rec), "broad standard")
  expect_output(print(rec), "hamielec")
  expect_output(print(rec), "linear")

  suppressWarnings(prepped <- recipes::prep(rec))
  expect_output(print(prepped), "Mn err")
  expect_output(print(prepped), "Mw err")
})

test_that("step_sec_broad_standard tidy method returns expected structure", {
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()
  broad_std <- create_broad_standard(mn = 50000, mw = 150000)

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_broad_standard(
      broad_standard = broad_std,
      known_mn = 50000,
      known_mw = 150000,
      extrapolation = "none"
    )

  suppressWarnings(prepped <- recipes::prep(rec))
  tidy_result <- recipes::tidy(prepped, number = 1)

  expect_s3_class(tidy_result, "tbl_df")

  # Expected columns
  expect_true("method" %in% names(tidy_result))
  expect_true("fit_type" %in% names(tidy_result))
  expect_true("known_mn" %in% names(tidy_result))
  expect_true("known_mw" %in% names(tidy_result))
  expect_true("known_dispersity" %in% names(tidy_result))
  expect_true("calculated_mn" %in% names(tidy_result))
  expect_true("calculated_mw" %in% names(tidy_result))
  expect_true("calculated_dispersity" %in% names(tidy_result))
  expect_true("mn_error_pct" %in% names(tidy_result))
  expect_true("mw_error_pct" %in% names(tidy_result))
  expect_true("dispersity_error_pct" %in% names(tidy_result))
  expect_true("calibration_min" %in% names(tidy_result))
  expect_true("calibration_max" %in% names(tidy_result))
  expect_true("coefficients" %in% names(tidy_result))

  # Validate values
  expect_equal(tidy_result$method, "hamielec")
  expect_equal(tidy_result$fit_type, "linear")
  expect_equal(tidy_result$known_mn, 50000)
  expect_equal(tidy_result$known_mw, 150000)
  expect_equal(tidy_result$known_dispersity, 3.0)
})

# -- Integration with other steps ----------------------------------------------

test_that("step_sec_broad_standard integrates with step_sec_mw_averages", {
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()
  broad_std <- create_broad_standard(mn = 50000, mw = 150000)

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_broad_standard(
      broad_standard = broad_std,
      known_mn = 50000,
      known_mw = 150000,
      output_col = "log_mw",
      extrapolation = "warn"
    ) |>
    step_sec_mw_averages(measures = "log_mw")

  suppressWarnings(prepped <- recipes::prep(rec))
  suppressWarnings(result <- recipes::bake(prepped, new_data = NULL))

  # Should have MW average columns
  expect_true("mw_mn" %in% names(result))
  expect_true("mw_mw" %in% names(result))
  expect_true("mw_dispersity" %in% names(result))
})

# -- Integral method tests -----------------------------------------------------

test_that("step_sec_broad_standard integral method requires reference_mwd", {
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()
  broad_std <- create_broad_standard(mn = 50000, mw = 150000)

  expect_error(
    recipes::recipe(~., data = test_data) |>
      step_sec_broad_standard(
        broad_standard = broad_std,
        known_mn = 50000,
        known_mw = 150000,
        method = "integral"
      ),
    "requires.*reference_mwd"
  )
})

test_that("step_sec_broad_standard validates reference_mwd columns", {
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()
  broad_std <- create_broad_standard(mn = 50000, mw = 150000)

  # Missing cumulative column
  bad_ref <- data.frame(mw = c(1000, 10000, 100000))
  expect_error(
    recipes::recipe(~., data = test_data) |>
      step_sec_broad_standard(
        broad_standard = broad_std,
        known_mn = 50000,
        known_mw = 150000,
        method = "integral",
        reference_mwd = bad_ref
      ),
    "mw.*cumulative"
  )
})

test_that("step_sec_broad_standard validates reference_mwd cumulative range", {
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()
  broad_std <- create_broad_standard(mn = 50000, mw = 150000)

  # Cumulative values > 1
  bad_ref <- data.frame(
    mw = c(1000, 10000, 100000),
    cumulative = c(0.1, 0.5, 1.5)
  )
  expect_error(
    recipes::recipe(~., data = test_data) |>
      step_sec_broad_standard(
        broad_standard = broad_std,
        known_mn = 50000,
        known_mw = 150000,
        method = "integral",
        reference_mwd = bad_ref
      ),
    "cumulative.*0 and 1"
  )
})

test_that("step_sec_broad_standard integral method works with valid reference", {
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()
  broad_std <- create_broad_standard(mn = 50000, mw = 150000)

  # Create a reference cumulative MWD matching the broad standard
  # Use log-normal distribution parameterized by Mn and Mw
  ref_mw <- 10^seq(3, 6, length.out = 50)
  # Approximation: for log-normal, sigma^2 = log(Mw/Mn)
  sigma_sq <- log(150000 / 50000)
  mu <- log(50000) + sigma_sq / 2 # Mn = exp(mu - sigma^2/2)
  ref_cum <- stats::plnorm(ref_mw, meanlog = mu, sdlog = sqrt(sigma_sq))

  reference_mwd <- data.frame(mw = ref_mw, cumulative = ref_cum)

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_broad_standard(
      broad_standard = broad_std,
      known_mn = 50000,
      known_mw = 150000,
      method = "integral",
      reference_mwd = reference_mwd
    )

  prepped <- recipes::prep(rec)
  expect_s3_class(prepped, "recipe")

  # Check that tidy returns correct method
  tidy_result <- recipes::tidy(prepped, number = 1)
  expect_equal(tidy_result$method, "integral")
})

# -- Additional coverage tests -------------------------------------------------

test_that("step_sec_broad_standard errors on insufficient data points", {
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()

  # Create standard with only 5 points (need at least 10)
  sparse_std <- data.frame(
    time = seq(14, 16, by = 0.5),
    signal = dnorm(seq(14, 16, by = 0.5), mean = 15, sd = 0.5)
  )

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_broad_standard(
      broad_standard = sparse_std,
      known_mn = 50000,
      known_mw = 150000
    )

  expect_error(
    recipes::prep(rec),
    "Insufficient data points"
  )
})

test_that("step_sec_broad_standard respects integration_range", {
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()
  broad_std <- create_broad_standard(mn = 50000, mw = 150000)

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_broad_standard(
      broad_standard = broad_std,
      known_mn = 50000,
      known_mw = 150000,
      integration_range = c(13, 17),
      extrapolation = "none"
    )

  suppressWarnings(prepped <- recipes::prep(rec))

  # Calibration range should match integration_range
  cal_range <- prepped$steps[[1]]$calibration_range
  expect_equal(cal_range[1], 13, tolerance = 0.1)
  expect_equal(cal_range[2], 17, tolerance = 0.1)
})

test_that("step_sec_broad_standard reports convergence diagnostics", {
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()
  broad_std <- create_broad_standard(mn = 50000, mw = 150000)

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_broad_standard(
      broad_standard = broad_std,
      known_mn = 50000,
      known_mw = 150000,
      extrapolation = "none"
    )

  suppressWarnings(prepped <- recipes::prep(rec))

  # Convergence should be 0 for successful optimization
  diag <- prepped$steps[[1]]$calibration_diagnostics
  expect_equal(diag$convergence, 0)
  expect_true(!is.null(diag$iterations))
  expect_true(!is.null(diag$final_objective))
})
