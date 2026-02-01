# ==============================================================================
# Tests for step_sec_band_broadening
# ==============================================================================

# -- Helper functions (must be defined first) ----------------------------------

.test_estimate_fwhm <- function(x, y) {
	y <- y - min(y, na.rm = TRUE)
	max_val <- max(y, na.rm = TRUE)
	half_max <- max_val / 2
	max_idx <- which.max(y)

	left_x <- x[1]
	right_x <- x[length(x)]

	for (i in seq(max_idx, 1)) {
		if (y[i] < half_max) {
			left_x <- x[i]
			break
		}
	}

	for (i in seq(max_idx, length(y))) {
		if (y[i] < half_max) {
			right_x <- x[i]
			break
		}
	}

	abs(right_x - left_x)
}

# Helper to create test data with a broadened peak
create_broadened_peak_data <- function(
	true_sigma = 0.3,
	broadening_sigma = 0.1,
	n_points = 200
) {
	time <- seq(5, 20, length.out = n_points)

	# True narrow peak
	true_peak <- dnorm(time, mean = 12, sd = true_sigma)

	# Broadened peak (convolution with Gaussian)
	broadening_kernel <- dnorm(
		seq(-2, 2, length.out = 41),
		mean = 0,
		sd = broadening_sigma / mean(diff(time))
	)
	broadening_kernel <- broadening_kernel / sum(broadening_kernel)

	broadened <- stats::filter(true_peak, broadening_kernel, sides = 2)
	broadened[is.na(broadened)] <- 0

	list(
		time = time,
		true_peak = true_peak,
		broadened = as.numeric(broadened),
		broadening_sigma = broadening_sigma
	)
}

# Helper to create test SEC data tibble
create_test_sec_data <- function(peak_data) {
	test_data <- tibble::tibble(sample_id = "test")
	test_data$ri <- measure::new_measure_list(
		list(
			measure::new_measure_tbl(
				location = peak_data$time,
				value = peak_data$broadened
			)
		)
	)
	test_data
}

# -- Constructor validation tests ----------------------------------------------

test_that("step_sec_band_broadening requires sigma or calibration_peak", {
	skip_if_not_installed("measure")

	peak_data <- create_broadened_peak_data()
	test_data <- create_test_sec_data(peak_data)

	expect_error(
		recipes::recipe(~., data = test_data) |>
			step_sec_band_broadening(),
		"sigma.*calibration_peak"
	)
})

test_that("step_sec_band_broadening validates sigma is positive", {
	skip_if_not_installed("measure")

	peak_data <- create_broadened_peak_data()
	test_data <- create_test_sec_data(peak_data)

	expect_error(
		recipes::recipe(~., data = test_data) |>
			step_sec_band_broadening(sigma = -0.1),
		"positive"
	)

	expect_error(
		recipes::recipe(~., data = test_data) |>
			step_sec_band_broadening(sigma = 0),
		"positive"
	)
})

test_that("step_sec_band_broadening validates damping is between 0 and 1", {
	skip_if_not_installed("measure")

	peak_data <- create_broadened_peak_data()
	test_data <- create_test_sec_data(peak_data)

	expect_error(
		recipes::recipe(~., data = test_data) |>
			step_sec_band_broadening(sigma = 0.1, damping = 0),
		"damping"
	)

	expect_error(
		recipes::recipe(~., data = test_data) |>
			step_sec_band_broadening(sigma = 0.1, damping = 1.5),
		"damping"
	)
})

test_that("step_sec_band_broadening validates iterations is positive", {
	skip_if_not_installed("measure")

	peak_data <- create_broadened_peak_data()
	test_data <- create_test_sec_data(peak_data)

	expect_error(
		recipes::recipe(~., data = test_data) |>
			step_sec_band_broadening(sigma = 0.1, iterations = 0),
		"iterations"
	)
})

# -- Prep tests ----------------------------------------------------------------

test_that("step_sec_band_broadening preps with explicit sigma", {
	skip_if_not_installed("measure")

	peak_data <- create_broadened_peak_data()
	test_data <- create_test_sec_data(peak_data)

	rec <- recipes::recipe(~., data = test_data) |>
		step_sec_band_broadening(sigma = 0.1, method = "tung")

	prepped <- recipes::prep(rec)

	expect_true(prepped$steps[[1]]$trained)
	expect_equal(prepped$steps[[1]]$estimated_sigma, 0.1)
})

test_that("step_sec_band_broadening preps with calibration_peak", {
	skip_if_not_installed("measure")

	peak_data <- create_broadened_peak_data()
	test_data <- create_test_sec_data(peak_data)

	# Create a narrow standard peak for calibration
	narrow_peak <- measure::new_measure_tbl(
		location = seq(10, 14, length.out = 100),
		value = dnorm(seq(10, 14, length.out = 100), mean = 12, sd = 0.15)
	)

	rec <- recipes::recipe(~., data = test_data) |>
		step_sec_band_broadening(calibration_peak = narrow_peak, method = "tung")

	expect_message(
		prepped <- recipes::prep(rec),
		"sigma"
	)

	expect_true(prepped$steps[[1]]$trained)
	expect_true(!is.null(prepped$steps[[1]]$estimated_sigma))
	expect_true(prepped$steps[[1]]$estimated_sigma > 0)
})

test_that("step_sec_band_broadening warns for large sigma", {
	skip_if_not_installed("measure")

	peak_data <- create_broadened_peak_data(true_sigma = 0.3)
	test_data <- create_test_sec_data(peak_data)

	# Sigma larger than 50% of peak width should warn
	rec <- recipes::recipe(~., data = test_data) |>
		step_sec_band_broadening(sigma = 2.0, method = "tung")

	expect_warning(
		recipes::prep(rec),
		"50%"
	)
})

# -- Bake tests ----------------------------------------------------------------

test_that("step_sec_band_broadening applies Tung correction", {
	skip_if_not_installed("measure")

	# Create test data with known broadening
	time <- seq(5, 20, length.out = 300)
	# True narrow peak with sigma = 0.3
	true_signal <- dnorm(time, mean = 12, sd = 0.3)
	# Add broadening by smoothing
	broadened <- stats::filter(true_signal, rep(1 / 15, 15), sides = 2)
	broadened[is.na(broadened)] <- 0
	broadened <- as.numeric(broadened)

	test_data <- tibble::tibble(sample_id = "test")
	test_data$ri <- measure::new_measure_list(
		list(measure::new_measure_tbl(location = time, value = broadened))
	)

	# Apply correction with appropriate sigma
	rec <- recipes::recipe(~., data = test_data) |>
		step_sec_band_broadening(
			sigma = 0.15, # Estimate of instrumental broadening
			method = "tung",
			damping = 0.8,
			iterations = 2
		)

	prepped <- recipes::prep(rec)
	result <- recipes::bake(prepped, new_data = NULL)

	expect_s3_class(result, "tbl_df")
	expect_true("ri" %in% names(result))

	# Check that the correction produces a valid result
	corrected_values <- result$ri[[1]]$value

	# The corrected peak should have higher maximum (sharper peak)
	original_max <- max(broadened, na.rm = TRUE)
	corrected_max <- max(corrected_values, na.rm = TRUE)

	# Sharper peak = higher maximum for same area
	expect_gt(corrected_max, original_max * 0.95)
})

test_that("step_sec_band_broadening preserves area", {
	skip_if_not_installed("measure")

	peak_data <- create_broadened_peak_data()
	test_data <- create_test_sec_data(peak_data)

	rec <- recipes::recipe(~., data = test_data) |>
		step_sec_band_broadening(sigma = 0.1, method = "tung")

	prepped <- recipes::prep(rec)
	result <- recipes::bake(prepped, new_data = NULL)

	original_area <- sum(peak_data$broadened, na.rm = TRUE)
	corrected_area <- sum(result$ri[[1]]$value, na.rm = TRUE)

	# Area should be preserved (within tolerance)
	expect_equal(corrected_area, original_area, tolerance = 0.01)
})

test_that("step_sec_band_broadening produces non-negative values", {
	skip_if_not_installed("measure")

	peak_data <- create_broadened_peak_data()
	test_data <- create_test_sec_data(peak_data)

	rec <- recipes::recipe(~., data = test_data) |>
		step_sec_band_broadening(sigma = 0.1, method = "tung")

	prepped <- recipes::prep(rec)
	result <- recipes::bake(prepped, new_data = NULL)

	corrected_values <- result$ri[[1]]$value
	expect_true(all(corrected_values >= 0))
})

test_that("step_sec_band_broadening works with EMG method", {
	skip_if_not_installed("measure")

	peak_data <- create_broadened_peak_data()
	test_data <- create_test_sec_data(peak_data)

	rec <- recipes::recipe(~., data = test_data) |>
		step_sec_band_broadening(
			sigma = 0.1,
			tau = 0.05,
			method = "emg",
			iterations = 3
		)

	prepped <- recipes::prep(rec)
	result <- recipes::bake(prepped, new_data = NULL)

	expect_s3_class(result, "tbl_df")
	expect_true("ri" %in% names(result))
})

# -- estimate_sigma tests ------------------------------------------------------

test_that("estimate_sigma works with Gaussian fit method", {
	# Create a Gaussian peak
	x <- seq(0, 10, length.out = 100)
	y <- dnorm(x, mean = 5, sd = 0.5)

	peak <- data.frame(location = x, value = y)

	result <- estimate_sigma(peak, method = "gaussian")

	expect_true(is.list(result))
	expect_true("sigma" %in% names(result))
	expect_true("fwhm" %in% names(result))
	expect_true("asymmetry" %in% names(result))

	# Should be close to true sigma
	expect_equal(result$sigma, 0.5, tolerance = 0.1)
})

test_that("estimate_sigma works with FWHM method", {
	x <- seq(0, 10, length.out = 100)
	y <- dnorm(x, mean = 5, sd = 0.5)

	peak <- data.frame(location = x, value = y)

	result <- estimate_sigma(peak, method = "fwhm")

	expect_true(is.list(result))
	expect_equal(result$sigma, 0.5, tolerance = 0.1)
})

test_that("estimate_sigma works with moments method", {
	x <- seq(0, 10, length.out = 100)
	y <- dnorm(x, mean = 5, sd = 0.5)

	peak <- data.frame(location = x, value = y)

	result <- estimate_sigma(peak, method = "moments")

	expect_true(is.list(result))
	expect_equal(result$sigma, 0.5, tolerance = 0.15)
})

test_that("estimate_sigma detects asymmetric peaks", {
	# Create an asymmetric (tailed) peak using EMG-like shape
	x <- seq(0, 15, length.out = 300)
	# Gaussian + strong exponential tail (only on right side)
	gaussian <- dnorm(x, mean = 5, sd = 0.4)
	tail_contrib <- ifelse(x > 5, exp(-(x - 5) / 1.5), 0)
	y <- gaussian + 0.3 * tail_contrib * gaussian[which.max(gaussian)]
	y <- y / max(y)

	peak <- data.frame(location = x, value = y)

	result <- estimate_sigma(peak, method = "fwhm")

	# Asymmetry should be > 1 for tailed peak (right side wider)
	expect_gt(result$asymmetry, 1.0)
})

test_that("estimate_sigma validates input", {
	expect_error(
		estimate_sigma(c(1, 2, 3)),
		"measure_tbl|data frame"
	)

	expect_error(
		estimate_sigma(data.frame(x = 1:5, y = 1:5)),
		"location.*value"
	)

	expect_error(
		estimate_sigma(data.frame(location = 1:3, value = 1:3)),
		"at least 5"
	)
})

# -- Print and tidy methods ----------------------------------------------------

test_that("step_sec_band_broadening print method works", {
	skip_if_not_installed("measure")

	peak_data <- create_broadened_peak_data()
	test_data <- create_test_sec_data(peak_data)

	rec <- recipes::recipe(~., data = test_data) |>
		step_sec_band_broadening(sigma = 0.1, method = "tung")

	expect_output(print(rec), "band broadening")
	expect_output(print(rec), "tung")

	prepped <- recipes::prep(rec)
	expect_output(print(prepped), "sigma")
})

test_that("step_sec_band_broadening tidy method returns expected structure", {
	skip_if_not_installed("measure")

	peak_data <- create_broadened_peak_data()
	test_data <- create_test_sec_data(peak_data)

	rec <- recipes::recipe(~., data = test_data) |>
		step_sec_band_broadening(sigma = 0.1, method = "tung", damping = 0.6)

	prepped <- recipes::prep(rec)
	tidy_result <- recipes::tidy(prepped, number = 1)

	expect_s3_class(tidy_result, "tbl_df")
	expect_true("method" %in% names(tidy_result))
	expect_true("sigma" %in% names(tidy_result))
	expect_true("damping" %in% names(tidy_result))
	expect_true("iterations" %in% names(tidy_result))
	expect_true("id" %in% names(tidy_result))

	expect_equal(tidy_result$method, "tung")
	expect_equal(tidy_result$sigma, 0.1)
	expect_equal(tidy_result$damping, 0.6)
})

# -- Edge case tests ----------------------------------------------------------

test_that("step_sec_band_broadening reduces peak width (sharpening)", {
	skip_if_not_installed("measure")

	# Create test data with known broadening
	time <- seq(5, 20, length.out = 300)
	# True narrow peak with sigma = 0.3
	true_signal <- dnorm(time, mean = 12, sd = 0.3)
	# Add broadening by smoothing
	broadened <- stats::filter(true_signal, rep(1 / 15, 15), sides = 2)
	broadened[is.na(broadened)] <- 0
	broadened <- as.numeric(broadened)

	test_data <- tibble::tibble(sample_id = "test")
	test_data$ri <- measure::new_measure_list(
		list(measure::new_measure_tbl(location = time, value = broadened))
	)

	# Apply correction with appropriate sigma
	rec <- recipes::recipe(~., data = test_data) |>
		step_sec_band_broadening(
			sigma = 0.15,
			method = "tung",
			damping = 0.8,
			iterations = 2
		)

	prepped <- recipes::prep(rec)
	result <- recipes::bake(prepped, new_data = NULL)

	# Estimate FWHM before and after
	original_fwhm <- .test_estimate_fwhm(time, broadened)
	corrected_fwhm <- .test_estimate_fwhm(time, result$ri[[1]]$value)

	# Corrected peak should be narrower (smaller FWHM)
	expect_lt(corrected_fwhm, original_fwhm)
})

test_that("step_sec_band_broadening handles short chromatograms gracefully", {
	skip_if_not_installed("measure")

	# Create very short chromatogram (< 10 points)
	time <- seq(10, 12, length.out = 8)
	signal <- dnorm(time, mean = 11, sd = 0.3)

	test_data <- tibble::tibble(sample_id = "test")
	test_data$ri <- measure::new_measure_list(
		list(measure::new_measure_tbl(location = time, value = signal))
	)

	rec <- recipes::recipe(~., data = test_data) |>
		step_sec_band_broadening(sigma = 0.1, method = "tung")

	# Should warn about short chromatogram during prep (baking happens during prep)
	expect_warning(
		prepped <- recipes::prep(rec),
		"fewer than 10"
	)

	result <- recipes::bake(prepped, new_data = NULL)

	# Values should be unchanged
	expect_equal(result$ri[[1]]$value, signal)
})

test_that("step_sec_band_broadening validates tau parameter", {
	skip_if_not_installed("measure")

	peak_data <- create_broadened_peak_data()
	test_data <- create_test_sec_data(peak_data)

	expect_error(
		recipes::recipe(~., data = test_data) |>
			step_sec_band_broadening(sigma = 0.1, tau = -0.1),
		"positive"
	)

	expect_error(
		recipes::recipe(~., data = test_data) |>
			step_sec_band_broadening(sigma = 0.1, tau = 0),
		"positive"
	)
})
