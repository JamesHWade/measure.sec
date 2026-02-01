# ==============================================================================
# Tests for step_sec_peaks_refine
# ==============================================================================

# Helper to create test data with baseline-corrected signal and detected peaks
create_test_refine_data <- function() {
	# Create SEC-like chromatogram: one clean Gaussian peak (no noise for
	# deterministic tests)
	time <- seq(10, 25, by = 0.05)
	signal <- 100 * dnorm(time, mean = 16, sd = 1.0)

	test_data <- tibble::tibble(sample_id = "test")
	test_data$ri <- measure::new_measure_list(
		list(measure::new_measure_tbl(location = time, value = signal))
	)

	# Manually create peaks with wide boundaries (simulating baseline-return)
	peaks <- measure:::new_peaks_tbl(
		peak_id = 1L,
		location = 16.0,
		height = max(signal),
		left_base = 10.0,
		right_base = 25.0,
		area = NA_real_
	)
	test_data$.peaks <- measure:::new_peaks_list(list(peaks))
	test_data
}

# Helper with multiple peaks
create_test_multi_peak_data <- function() {
	time <- seq(10, 30, by = 0.05)
	peak1 <- 100 * dnorm(time, mean = 14, sd = 0.8)
	peak2 <- 150 * dnorm(time, mean = 20, sd = 1.2)
	signal <- peak1 + peak2

	test_data <- tibble::tibble(sample_id = "test")
	test_data$ri <- measure::new_measure_list(
		list(measure::new_measure_tbl(location = time, value = signal))
	)

	peaks <- measure:::new_peaks_tbl(
		peak_id = c(1L, 2L),
		location = c(14.0, 20.0),
		height = c(max(100 * dnorm(0, sd = 0.8)), max(150 * dnorm(0, sd = 1.2))),
		left_base = c(10.0, 16.5),
		right_base = c(16.5, 30.0),
		area = c(NA_real_, NA_real_)
	)
	test_data$.peaks <- measure:::new_peaks_list(list(peaks))
	test_data
}

# Helper with no peaks detected
create_test_no_peaks_data <- function() {
	time <- seq(10, 25, by = 0.05)
	signal <- rep(10, length(time))

	test_data <- tibble::tibble(sample_id = "test")
	test_data$ri <- measure::new_measure_list(
		list(measure::new_measure_tbl(location = time, value = signal))
	)

	peaks <- measure:::new_peaks_tbl()
	test_data$.peaks <- measure:::new_peaks_list(list(peaks))
	test_data
}

# -- Constructor validation tests ----------------------------------------------

test_that("step_sec_peaks_refine validates method parameter", {
	skip_if_not_installed("measure")

	test_data <- create_test_refine_data()

	expect_error(
		recipes::recipe(~., data = test_data) |>
			step_sec_peaks_refine(method = "invalid"),
		"must be one of"
	)
})

test_that("step_sec_peaks_refine validates cutoff parameter", {
	skip_if_not_installed("measure")

	test_data <- create_test_refine_data()

	expect_error(
		recipes::recipe(~., data = test_data) |>
			step_sec_peaks_refine(cutoff = 0),
		"cutoff.*between 0 and 1"
	)

	expect_error(
		recipes::recipe(~., data = test_data) |>
			step_sec_peaks_refine(cutoff = 1),
		"cutoff.*between 0 and 1"
	)

	expect_error(
		recipes::recipe(~., data = test_data) |>
			step_sec_peaks_refine(cutoff = -0.5),
		"cutoff.*between 0 and 1"
	)

	expect_error(
		recipes::recipe(~., data = test_data) |>
			step_sec_peaks_refine(cutoff = "abc"),
		"cutoff.*between 0 and 1"
	)
})

# -- Prep validation tests -----------------------------------------------------

test_that("step_sec_peaks_refine errors when peaks column is missing", {
	skip_if_not_installed("measure")

	test_data <- tibble::tibble(sample_id = "test")
	test_data$ri <- measure::new_measure_list(
		list(
			measure::new_measure_tbl(
				location = seq(10, 25, by = 0.05),
				value = rnorm(301)
			)
		)
	)

	expect_error(
		recipes::recipe(~., data = test_data) |>
			step_sec_peaks_refine() |>
			recipes::prep(),
		"not found.*step_sec_peaks_detect"
	)
})

test_that("step_sec_peaks_refine auto-detects measure column", {
	skip_if_not_installed("measure")

	test_data <- create_test_refine_data()

	# Default measures_col = ".measures" won't be found, should fall back to "ri"
	rec <- recipes::recipe(~., data = test_data) |>
		step_sec_peaks_refine()

	expect_message(
		prepped <- recipes::prep(rec),
		"Using.*ri.*instead"
	)
	expect_s3_class(prepped, "recipe")
})

test_that("step_sec_peaks_refine errors when no measure columns exist", {
	skip_if_not_installed("measure")

	test_data <- tibble::tibble(sample_id = "test")
	peaks <- measure:::new_peaks_tbl(
		peak_id = 1L,
		location = 16.0,
		height = 100,
		left_base = 10.0,
		right_base = 25.0,
		area = NA_real_
	)
	test_data$.peaks <- measure:::new_peaks_list(list(peaks))

	expect_error(
		recipes::recipe(~., data = test_data) |>
			step_sec_peaks_refine() |>
			recipes::prep(),
		"No measure columns"
	)
})

# -- Bake tests ----------------------------------------------------------------

test_that("step_sec_peaks_refine tightens boundaries on single peak", {
	skip_if_not_installed("measure")

	test_data <- create_test_refine_data()

	rec <- recipes::recipe(~., data = test_data) |>
		step_sec_peaks_refine(measures_col = "ri", cutoff = 0.005)

	result <- recipes::prep(rec) |>
		recipes::bake(new_data = NULL)

	peaks <- result$.peaks[[1]]

	# Refined boundaries should be tighter than originals
	expect_gt(peaks$left_base[1], peaks$original_left_base[1])
	expect_lt(peaks$right_base[1], peaks$original_right_base[1])

	# Original boundaries should match what we set

	expect_equal(peaks$original_left_base[1], 10.0)
	expect_equal(peaks$original_right_base[1], 25.0)
})

test_that("step_sec_peaks_refine handles multiple peaks", {
	skip_if_not_installed("measure")

	test_data <- create_test_multi_peak_data()

	rec <- recipes::recipe(~., data = test_data) |>
		step_sec_peaks_refine(measures_col = "ri", cutoff = 0.005)

	result <- recipes::prep(rec) |>
		recipes::bake(new_data = NULL)

	peaks <- result$.peaks[[1]]

	# Both peaks should have tightened boundaries
	expect_equal(nrow(peaks), 2)
	for (j in seq_len(nrow(peaks))) {
		expect_gte(peaks$left_base[j], peaks$original_left_base[j])
		expect_lte(peaks$right_base[j], peaks$original_right_base[j])
	}
})

test_that("step_sec_peaks_refine preserves empty peaks", {
	skip_if_not_installed("measure")

	test_data <- create_test_no_peaks_data()

	rec <- recipes::recipe(~., data = test_data) |>
		step_sec_peaks_refine(measures_col = "ri")

	result <- recipes::prep(rec) |>
		recipes::bake(new_data = NULL)

	peaks <- result$.peaks[[1]]
	expect_equal(nrow(peaks), 0)
})

test_that("step_sec_peaks_refine preserves peaks_list class", {
	skip_if_not_installed("measure")

	test_data <- create_test_refine_data()

	rec <- recipes::recipe(~., data = test_data) |>
		step_sec_peaks_refine(measures_col = "ri")

	result <- recipes::prep(rec) |>
		recipes::bake(new_data = NULL)

	expect_true(measure::is_peaks_list(result$.peaks))
})

test_that("step_sec_peaks_refine adds original boundary columns", {
	skip_if_not_installed("measure")

	test_data <- create_test_refine_data()

	rec <- recipes::recipe(~., data = test_data) |>
		step_sec_peaks_refine(measures_col = "ri")

	result <- recipes::prep(rec) |>
		recipes::bake(new_data = NULL)

	peaks <- result$.peaks[[1]]
	expect_true("original_left_base" %in% names(peaks))
	expect_true("original_right_base" %in% names(peaks))
})

test_that("higher cutoff produces tighter boundaries", {
	skip_if_not_installed("measure")

	test_data <- create_test_refine_data()

	# Loose cutoff
	rec_loose <- recipes::recipe(~., data = test_data) |>
		step_sec_peaks_refine(measures_col = "ri", cutoff = 0.001)
	result_loose <- recipes::prep(rec_loose) |>
		recipes::bake(new_data = NULL)

	# Tight cutoff
	rec_tight <- recipes::recipe(~., data = test_data) |>
		step_sec_peaks_refine(measures_col = "ri", cutoff = 0.05)
	result_tight <- recipes::prep(rec_tight) |>
		recipes::bake(new_data = NULL)

	peaks_loose <- result_loose$.peaks[[1]]
	peaks_tight <- result_tight$.peaks[[1]]

	# Tight cutoff should have narrower boundaries
	expect_gte(peaks_tight$left_base[1], peaks_loose$left_base[1])
	expect_lte(peaks_tight$right_base[1], peaks_loose$right_base[1])
})

# -- Edge case tests -----------------------------------------------------------

test_that("step_sec_peaks_refine handles peak with negative apex", {
	skip_if_not_installed("measure")

	time <- seq(10, 25, by = 0.05)
	signal <- rep(-5, length(time))

	test_data <- tibble::tibble(sample_id = "test")
	test_data$ri <- measure::new_measure_list(
		list(measure::new_measure_tbl(location = time, value = signal))
	)

	peaks <- measure:::new_peaks_tbl(
		peak_id = 1L,
		location = 16.0,
		height = 0,
		left_base = 12.0,
		right_base = 20.0,
		area = NA_real_
	)
	test_data$.peaks <- measure:::new_peaks_list(list(peaks))

	rec <- recipes::recipe(~., data = test_data) |>
		step_sec_peaks_refine(measures_col = "ri")

	result <- recipes::prep(rec) |>
		recipes::bake(new_data = NULL)

	# Boundaries should be unchanged (apex <= 0)
	peaks_out <- result$.peaks[[1]]
	expect_equal(peaks_out$left_base[1], 12.0)
	expect_equal(peaks_out$right_base[1], 20.0)
})

test_that("step_sec_peaks_refine handles peak with few data points", {
	skip_if_not_installed("measure")

	# Only 4 points in the peak region (< 5 minimum)
	time <- c(15.0, 15.5, 16.0, 16.5)
	signal <- c(10, 50, 100, 50)

	test_data <- tibble::tibble(sample_id = "test")
	test_data$ri <- measure::new_measure_list(
		list(measure::new_measure_tbl(location = time, value = signal))
	)

	peaks <- measure:::new_peaks_tbl(
		peak_id = 1L,
		location = 16.0,
		height = 100,
		left_base = 15.0,
		right_base = 16.5,
		area = NA_real_
	)
	test_data$.peaks <- measure:::new_peaks_list(list(peaks))

	rec <- recipes::recipe(~., data = test_data) |>
		step_sec_peaks_refine(measures_col = "ri")

	result <- recipes::prep(rec) |>
		recipes::bake(new_data = NULL)

	# Boundaries should be unchanged (< 5 points)
	peaks_out <- result$.peaks[[1]]
	expect_equal(peaks_out$left_base[1], 15.0)
	expect_equal(peaks_out$right_base[1], 16.5)
})

# -- Tidy method tests ---------------------------------------------------------

test_that("tidy.step_sec_peaks_refine returns expected tibble", {
	skip_if_not_installed("measure")

	test_data <- create_test_refine_data()

	rec <- recipes::recipe(~., data = test_data) |>
		step_sec_peaks_refine(measures_col = "ri", cutoff = 0.01)

	# Before prep
	tidy_result <- recipes::tidy(rec, number = 1)
	expect_s3_class(tidy_result, "tbl_df")
	expect_true("method" %in% names(tidy_result))
	expect_equal(tidy_result$method, "height_fraction")
	expect_equal(tidy_result$cutoff, 0.01)

	# After prep
	prepped <- recipes::prep(rec)
	tidy_prepped <- recipes::tidy(prepped, number = 1)
	expect_s3_class(tidy_prepped, "tbl_df")
	expect_equal(tidy_prepped$method, "height_fraction")
	expect_equal(tidy_prepped$measures_col, "ri")
})

# -- Print method tests --------------------------------------------------------

test_that("print.step_sec_peaks_refine works", {
	skip_if_not_installed("measure")

	test_data <- create_test_refine_data()

	rec <- recipes::recipe(~., data = test_data) |>
		step_sec_peaks_refine(cutoff = 0.005)

	expect_output(print(rec), "height_fraction")
	expect_output(print(rec), "0.005")
})

# -- required_pkgs tests -------------------------------------------------------

test_that("required_pkgs.step_sec_peaks_refine lists correct packages", {
	skip_if_not_installed("measure")

	test_data <- create_test_refine_data()

	rec <- recipes::recipe(~., data = test_data) |>
		step_sec_peaks_refine()

	pkgs <- recipes::required_pkgs(rec)
	expect_true("measure.sec" %in% pkgs)
	expect_true("measure" %in% pkgs)
})
