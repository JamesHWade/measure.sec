# ==============================================================================
# Edge case tests for SEC steps
# ==============================================================================

# -- step_sec_baseline edge cases ----------------------------------------------

test_that("step_sec_baseline validates left_frac parameter", {
	skip_if_not_installed("measure")

	time <- seq(5, 25, by = 0.1)
	signal <- dnorm(time, mean = 15, sd = 1)

	test_data <- tibble::tibble(sample_id = "test")
	test_data$ri <- measure::new_measure_list(
		list(measure::new_measure_tbl(location = time, value = signal))
	)

	# left_frac <= 0 - validation happens at prep time
	expect_error(
		recipes::recipe(~., data = test_data) |>
			step_sec_baseline(measures = "ri", left_frac = 0) |>
			recipes::prep(),
		"between 0 and 0.5"
	)

	# left_frac >= 0.5
	expect_error(
		recipes::recipe(~., data = test_data) |>
			step_sec_baseline(measures = "ri", left_frac = 0.5) |>
			recipes::prep(),
		"between 0 and 0.5"
	)
})

test_that("step_sec_baseline validates right_frac parameter", {
	skip_if_not_installed("measure")

	time <- seq(5, 25, by = 0.1)
	signal <- dnorm(time, mean = 15, sd = 1)

	test_data <- tibble::tibble(sample_id = "test")
	test_data$ri <- measure::new_measure_list(
		list(measure::new_measure_tbl(location = time, value = signal))
	)

	# right_frac <= 0 - validation at prep time
	expect_error(
		recipes::recipe(~., data = test_data) |>
			step_sec_baseline(measures = "ri", right_frac = 0) |>
			recipes::prep(),
		"between 0 and 0.5"
	)
})

test_that("step_sec_baseline validates overlapping fractions", {
	skip_if_not_installed("measure")

	time <- seq(5, 25, by = 0.1)
	signal <- dnorm(time, mean = 15, sd = 1)

	test_data <- tibble::tibble(sample_id = "test")
	test_data$ri <- measure::new_measure_list(
		list(measure::new_measure_tbl(location = time, value = signal))
	)

	# Sum of fractions >= 0.5 - validation at prep time
	expect_error(
		recipes::recipe(~., data = test_data) |>
			step_sec_baseline(measures = "ri", left_frac = 0.3, right_frac = 0.3) |>
			recipes::prep(),
		"less than 0.5"
	)
})

test_that("step_sec_baseline handles short chromatograms", {
	skip_if_not_installed("measure")

	# Very short chromatogram (< 10 points)
	time <- seq(1, 5, by = 1) # Only 5 points
	signal <- c(0.1, 0.5, 1.0, 0.5, 0.1)

	test_data <- tibble::tibble(sample_id = "test")
	test_data$ri <- measure::new_measure_list(
		list(measure::new_measure_tbl(location = time, value = signal))
	)

	rec <- recipes::recipe(~., data = test_data) |>
		step_sec_baseline(measures = "ri")

	# Warning emitted during prep, data returned unchanged
	expect_warning(prepped <- recipes::prep(rec), "fewer than 10 points")
	result <- recipes::bake(prepped, new_data = NULL)

	# Original values should be unchanged
	expect_equal(result$ri[[1]]$value, signal)
})

test_that("step_sec_baseline handles all-NA values", {
	skip_if_not_installed("measure")

	time <- seq(5, 25, by = 0.1)
	signal <- rep(NA_real_, length(time))

	test_data <- tibble::tibble(sample_id = "test")
	test_data$ri <- measure::new_measure_list(
		list(measure::new_measure_tbl(location = time, value = signal))
	)

	rec <- recipes::recipe(~., data = test_data) |>
		step_sec_baseline(measures = "ri")

	# Warning emitted during prep
	expect_warning(prepped <- recipes::prep(rec), "all NA")
	result <- recipes::bake(prepped, new_data = NULL)

	# Values should still be all NA
	expect_true(all(is.na(result$ri[[1]]$value)))
})

# -- step_sec_mw_averages edge cases -------------------------------------------

test_that("step_sec_mw_averages handles zero-weight chromatograms", {
	skip_if_not_installed("measure")

	log_mw <- seq(3, 6, by = 0.05)
	conc <- rep(0, length(log_mw)) # All zeros

	test_data <- tibble::tibble(sample_id = "test")
	test_data$ri <- measure::new_measure_list(
		list(measure::new_measure_tbl(location = log_mw, value = conc))
	)

	rec <- recipes::recipe(~., data = test_data) |>
		step_sec_mw_averages(measures = "ri")

	prepped <- recipes::prep(rec)
	result <- recipes::bake(prepped, new_data = NULL)

	# Should return NA values, not errors
	expect_true(is.na(result$mw_mn[1]))
	expect_true(is.na(result$mw_mw[1]))
	expect_true(is.na(result$mw_mz[1]))
})

test_that("step_sec_mw_averages handles single valid point", {
	skip_if_not_installed("measure")

	log_mw <- seq(3, 6, by = 0.05)
	conc <- rep(0, length(log_mw))
	conc[30] <- 1.0 # Only one non-zero point

	test_data <- tibble::tibble(sample_id = "test")
	test_data$ri <- measure::new_measure_list(
		list(measure::new_measure_tbl(location = log_mw, value = conc))
	)

	rec <- recipes::recipe(~., data = test_data) |>
		step_sec_mw_averages(measures = "ri")

	prepped <- recipes::prep(rec)
	result <- recipes::bake(prepped, new_data = NULL)

	# Should return NA values for insufficient data
	expect_true(is.na(result$mw_mn[1]))
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

# -- step_sec_mw_fractions edge cases ------------------------------------------

test_that("step_sec_mw_fractions works with single cutoff", {
	skip_if_not_installed("measure")

	log_mw <- seq(3, 6, by = 0.05)
	conc <- dnorm(log_mw, mean = 4.5, sd = 0.5)

	test_data <- tibble::tibble(sample_id = "test")
	test_data$ri <- measure::new_measure_list(
		list(measure::new_measure_tbl(location = log_mw, value = conc))
	)

	rec <- recipes::recipe(~., data = test_data) |>
		step_sec_mw_fractions(measures = "ri", cutoffs = 50000)

	prepped <- recipes::prep(rec)
	result <- recipes::bake(prepped, new_data = NULL)

	expect_s3_class(result, "tbl_df")
	# Should have fraction columns
	frac_cols <- names(result)[grepl("frac", names(result), fixed = TRUE)]
	expect_true(length(frac_cols) > 0)
})

# -- step_sec_mw_distribution edge cases ---------------------------------------

test_that("step_sec_mw_distribution validates n_points", {
	skip_if_not_installed("measure")

	log_mw <- seq(3, 6, by = 0.05)
	conc <- dnorm(log_mw, mean = 4.5, sd = 0.5)

	test_data <- tibble::tibble(sample_id = "test")
	test_data$ri <- measure::new_measure_list(
		list(measure::new_measure_tbl(location = log_mw, value = conc))
	)

	expect_error(
		recipes::recipe(~., data = test_data) |>
			step_sec_mw_distribution(measures = "ri", n_points = 5),
		"n_points"
	)
})

test_that("step_sec_mw_distribution generates cumulative distribution", {
	skip_if_not_installed("measure")

	log_mw <- seq(3, 6, by = 0.05)
	conc <- dnorm(log_mw, mean = 4.5, sd = 0.5)

	test_data <- tibble::tibble(sample_id = "test")
	test_data$ri <- measure::new_measure_list(
		list(measure::new_measure_tbl(location = log_mw, value = conc))
	)

	rec <- recipes::recipe(~., data = test_data) |>
		step_sec_mw_distribution(measures = "ri", type = "cumulative")

	prepped <- recipes::prep(rec)
	result <- recipes::bake(prepped, new_data = NULL)

	expect_s3_class(result, "tbl_df")
	expect_true("ri" %in% names(result))

	# Cumulative distribution should have values between 0 and 1
	values <- result$ri[[1]]$value
	expect_true(all(values >= 0 & values <= 1, na.rm = TRUE))
})

# -- QC function edge cases ----------------------------------------------------

test_that("measure_sec_asymmetry rejects non-positive inputs", {
	expect_error(measure_sec_asymmetry(0, 0.3), "positive")
	expect_error(measure_sec_asymmetry(-0.2, 0.3), "positive")
	expect_error(measure_sec_asymmetry(0.2, -0.3), "positive")
	expect_error(measure_sec_asymmetry(0.2, 0), "positive")
})

test_that("measure_sec_plate_count validates width_type", {
	# Invalid width type - match.arg will throw an error
	expect_error(
		measure_sec_plate_count(10.0, 0.5, width_type = "invalid")
	)
})

# -- step_sec_aggregates edge cases --------------------------------------------

test_that("step_sec_aggregates requires boundaries for manual method", {
	skip_if_not_installed("measure")

	time <- seq(5, 20, by = 0.1)
	signal <- dnorm(time, mean = 12, sd = 0.5)

	test_data <- tibble::tibble(sample_id = "test")
	test_data$uv <- measure::new_measure_list(
		list(measure::new_measure_tbl(location = time, value = signal))
	)

	# Missing monomer_end
	expect_error(
		recipes::recipe(~., data = test_data) |>
			step_sec_aggregates(
				measures = "uv",
				method = "manual",
				monomer_start = 10
			),
		"monomer_start.*monomer_end"
	)

	# Missing monomer_start
	expect_error(
		recipes::recipe(~., data = test_data) |>
			step_sec_aggregates(measures = "uv", method = "manual", monomer_end = 14),
		"monomer_start.*monomer_end"
	)
})

test_that("step_sec_aggregates validates hmws_threshold", {
	skip_if_not_installed("measure")

	time <- seq(5, 20, by = 0.1)
	signal <- dnorm(time, mean = 12, sd = 0.5)

	test_data <- tibble::tibble(sample_id = "test")
	test_data$uv <- measure::new_measure_list(
		list(measure::new_measure_tbl(location = time, value = signal))
	)

	# hmws_threshold < 0
	expect_error(
		recipes::recipe(~., data = test_data) |>
			step_sec_aggregates(measures = "uv", hmws_threshold = -0.1),
		"between 0 and 1"
	)

	# hmws_threshold > 1
	expect_error(
		recipes::recipe(~., data = test_data) |>
			step_sec_aggregates(measures = "uv", hmws_threshold = 1.5),
		"between 0 and 1"
	)
})
