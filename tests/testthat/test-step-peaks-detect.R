# ==============================================================================
# Tests for step_sec_peaks_detect
# ==============================================================================

# Helper to create test data with measure columns containing peaks
create_test_sec_peaks_data <- function() {
	# Create SEC-like chromatogram with 2-3 peaks
	time <- seq(10, 25, by = 0.05)
	# Three Gaussian peaks
	peak1 <- 100 * dnorm(time, mean = 13, sd = 0.8)
	peak2 <- 150 * dnorm(time, mean = 16, sd = 1.0)
	peak3 <- 50 * dnorm(time, mean = 20, sd = 0.6)
	signal <- peak1 + peak2 + peak3 + rnorm(length(time), sd = 1)

	test_data <- tibble::tibble(sample_id = "test")
	test_data$ri <- measure::new_measure_list(
		list(measure::new_measure_tbl(location = time, value = signal))
	)
	test_data
}

# Helper to create flat signal (no peaks)
create_flat_signal_data <- function() {
	time <- seq(10, 25, by = 0.05)
	signal <- rep(10, length(time)) + rnorm(length(time), sd = 0.1)

	test_data <- tibble::tibble(sample_id = "test")
	test_data$ri <- measure::new_measure_list(
		list(measure::new_measure_tbl(location = time, value = signal))
	)
	test_data
}

# -- Constructor validation tests ----------------------------------------------

test_that("step_sec_peaks_detect validates algorithm parameter", {
	skip_if_not_installed("measure")

	test_data <- create_test_sec_peaks_data()

	# Only finderskeepers is supported
	expect_error(
		recipes::recipe(~., data = test_data) |>
			step_sec_peaks_detect(algorithm = "prominence"),
		"finderskeepers.*supported"
	)
})

test_that("step_sec_peaks_detect validates min_height parameter", {
	skip_if_not_installed("measure")

	test_data <- create_test_sec_peaks_data()

	expect_error(
		recipes::recipe(~., data = test_data) |>
			step_sec_peaks_detect(min_height = -1),
		"min_height.*non-negative"
	)
})

test_that("step_sec_peaks_detect validates loess_span parameter", {
	skip_if_not_installed("measure")

	test_data <- create_test_sec_peaks_data()

	expect_error(
		recipes::recipe(~., data = test_data) |>
			step_sec_peaks_detect(loess_span = 1.5),
		"loess_span.*between 0 and 1"
	)

	expect_error(
		recipes::recipe(~., data = test_data) |>
			step_sec_peaks_detect(loess_span = 0),
		"loess_span.*between 0 and 1"
	)
})

test_that("step_sec_peaks_detect validates ist_points parameter", {
	skip_if_not_installed("measure")

	test_data <- create_test_sec_peaks_data()

	expect_error(
		recipes::recipe(~., data = test_data) |>
			step_sec_peaks_detect(ist_points = 5),
		"ist_points.*at least 10"
	)
})

# -- Prep and bake tests -------------------------------------------------------

test_that("step_sec_peaks_detect can prep with finderskeepers", {
	skip_if_not_installed("measure")
	skip_if_not_installed("changepoint")

	test_data <- create_test_sec_peaks_data()

	rec <- recipes::recipe(~., data = test_data) |>
		step_sec_peaks_detect(
			algorithm = "finderskeepers",
			min_height = 5,
			loess_span = 0.1 # Larger span for test data
		)

	# Suppress LOESS warnings about near-singularities
	prepped <- suppressWarnings(recipes::prep(rec))
	expect_s3_class(prepped, "recipe")
})

test_that("step_sec_peaks_detect bake creates .peaks column", {
	skip_if_not_installed("measure")
	skip_if_not_installed("changepoint")

	test_data <- create_test_sec_peaks_data()

	rec <- recipes::recipe(~., data = test_data) |>
		step_sec_peaks_detect(
			algorithm = "finderskeepers",
			min_height = 5,
			loess_span = 0.1
		)

	result <- suppressWarnings(
		recipes::prep(rec) |>
			recipes::bake(new_data = NULL)
	)

	expect_true(".peaks" %in% names(result))
	expect_true(measure::is_peaks_list(result$.peaks))
})

test_that("step_sec_peaks_detect detects multiple peaks", {
	skip_if_not_installed("measure")
	skip_if_not_installed("changepoint")

	test_data <- create_test_sec_peaks_data()

	rec <- recipes::recipe(~., data = test_data) |>
		step_sec_peaks_detect(
			algorithm = "finderskeepers",
			min_height = 3,
			loess_span = 0.1
		)

	result <- suppressWarnings(
		recipes::prep(rec) |>
			recipes::bake(new_data = NULL)
	)

	peaks <- result$.peaks[[1]]
	# Should detect 2-3 peaks in our test data
	expect_gte(nrow(peaks), 2)
	expect_lte(nrow(peaks), 5)
})

test_that("step_sec_peaks_detect returns empty peaks_tbl for flat signal", {
	skip_if_not_installed("measure")
	skip_if_not_installed("changepoint")

	test_data <- create_flat_signal_data()

	rec <- recipes::recipe(~., data = test_data) |>
		step_sec_peaks_detect(
			algorithm = "finderskeepers",
			min_height = 5,
			loess_span = 0.1
		)

	result <- suppressWarnings(
		recipes::prep(rec) |>
			recipes::bake(new_data = NULL)
	)

	peaks <- result$.peaks[[1]]
	expect_equal(nrow(peaks), 0)
})

# -- Peak output structure tests -----------------------------------------------

test_that("step_sec_peaks_detect output has correct columns", {
	skip_if_not_installed("measure")
	skip_if_not_installed("changepoint")

	test_data <- create_test_sec_peaks_data()

	rec <- recipes::recipe(~., data = test_data) |>
		step_sec_peaks_detect(
			algorithm = "finderskeepers",
			min_height = 3,
			loess_span = 0.1
		)

	result <- suppressWarnings(
		recipes::prep(rec) |>
			recipes::bake(new_data = NULL)
	)

	peaks <- result$.peaks[[1]]
	expected_cols <- c(
		"peak_id",
		"location",
		"height",
		"left_base",
		"right_base",
		"area"
	)
	expect_true(all(expected_cols %in% names(peaks)))
})

test_that("step_sec_peaks_detect peak locations are within data range", {
	skip_if_not_installed("measure")
	skip_if_not_installed("changepoint")

	test_data <- create_test_sec_peaks_data()
	input_location <- test_data$ri[[1]]$location

	rec <- recipes::recipe(~., data = test_data) |>
		step_sec_peaks_detect(
			algorithm = "finderskeepers",
			min_height = 3,
			loess_span = 0.1
		)

	result <- suppressWarnings(
		recipes::prep(rec) |>
			recipes::bake(new_data = NULL)
	)

	peaks <- result$.peaks[[1]]
	if (nrow(peaks) > 0) {
		expect_true(all(peaks$location >= min(input_location)))
		expect_true(all(peaks$location <= max(input_location)))
	}
})

# -- Tidy method tests ---------------------------------------------------------

test_that("tidy.step_sec_peaks_detect returns expected tibble", {
	skip_if_not_installed("measure")
	skip_if_not_installed("changepoint")

	test_data <- create_test_sec_peaks_data()

	rec <- recipes::recipe(~., data = test_data) |>
		step_sec_peaks_detect(
			algorithm = "finderskeepers",
			min_height = 5,
			loess_span = 0.1
		)

	# Before prep
	tidy_result <- recipes::tidy(rec, number = 1)
	expect_s3_class(tidy_result, "tbl_df")
	expect_true("algorithm" %in% names(tidy_result))
	expect_equal(tidy_result$algorithm, "finderskeepers")

	# After prep
	prepped <- suppressWarnings(recipes::prep(rec))
	tidy_prepped <- recipes::tidy(prepped, number = 1)
	expect_s3_class(tidy_prepped, "tbl_df")
	expect_equal(tidy_prepped$algorithm, "finderskeepers")
})

# -- Print method tests --------------------------------------------------------

test_that("print.step_sec_peaks_detect works", {
	skip_if_not_installed("measure")

	test_data <- create_test_sec_peaks_data()

	rec <- recipes::recipe(~., data = test_data) |>
		step_sec_peaks_detect(algorithm = "finderskeepers")

	expect_output(print(rec), "finderskeepers")
})

# -- required_pkgs tests -------------------------------------------------------

test_that("required_pkgs.step_sec_peaks_detect lists changepoint for finderskeepers", {
	skip_if_not_installed("measure")

	test_data <- create_test_sec_peaks_data()

	rec <- recipes::recipe(~., data = test_data) |>
		step_sec_peaks_detect(algorithm = "finderskeepers")

	pkgs <- recipes::required_pkgs(rec)
	expect_true("changepoint" %in% pkgs)
})
