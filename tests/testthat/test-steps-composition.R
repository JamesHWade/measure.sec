# ==============================================================================
# Tests for composition and protein SEC steps
# ==============================================================================

test_that("step_sec_uv_ri_ratio calculates ratio correctly", {
	skip_if_not_installed("measure")

	time <- seq(5, 15, by = 0.1)
	ri_signal <- dnorm(time, mean = 10, sd = 0.5)
	uv_signal <- 0.5 * dnorm(time, mean = 10, sd = 0.5) # Half the RI

	test_data <- tibble::tibble(sample_id = "test")

	test_data$ri <- measure::new_measure_list(
		list(measure::new_measure_tbl(location = time, value = ri_signal))
	)
	test_data$uv <- measure::new_measure_list(
		list(measure::new_measure_tbl(location = time, value = uv_signal))
	)

	rec <- recipes::recipe(~., data = test_data) |>
		step_sec_uv_ri_ratio(uv_col = "uv", ri_col = "ri", smooth = FALSE)

	prepped <- recipes::prep(rec)
	result <- recipes::bake(prepped, new_data = NULL)

	expect_s3_class(result, "tbl_df")
	expect_true("uv_ri_ratio" %in% names(result))

	# Ratio should be approximately 0.5 at the peak
	ratio_vals <- result$uv_ri_ratio[[1]]$value
	peak_idx <- which.max(ri_signal)
	expect_equal(ratio_vals[peak_idx], 0.5, tolerance = 0.1)
})

test_that("step_sec_composition calculates weight fraction", {
	skip_if_not_installed("measure")

	time <- seq(5, 15, by = 0.1)

	# Simulate a copolymer with varying composition
	# Component A: high UV, high RI
	# Component B: low UV, low RI
	ri_signal <- dnorm(time, mean = 10, sd = 0.5)
	uv_signal <- 0.8 * dnorm(time, mean = 10, sd = 0.5)

	test_data <- tibble::tibble(sample_id = "test")

	test_data$ri <- measure::new_measure_list(
		list(measure::new_measure_tbl(location = time, value = ri_signal))
	)
	test_data$uv <- measure::new_measure_list(
		list(measure::new_measure_tbl(location = time, value = uv_signal))
	)

	rec <- recipes::recipe(~., data = test_data) |>
		step_sec_composition(
			uv_col = "uv",
			ri_col = "ri",
			component_a_uv = 1.0,
			component_a_ri = 0.185,
			component_b_uv = 0.1,
			component_b_ri = 0.084
		)

	prepped <- recipes::prep(rec)
	result <- recipes::bake(prepped, new_data = NULL)

	expect_s3_class(result, "tbl_df")
	expect_true("composition_a" %in% names(result))

	# Composition should be between 0 and 1
	comp_vals <- result$composition_a[[1]]$value
	valid_comp <- comp_vals[!is.na(comp_vals)]
	expect_true(all(valid_comp >= 0 & valid_comp <= 1))
})

test_that("step_sec_aggregates quantifies HMWS/monomer/LMWS", {
	skip_if_not_installed("measure")

	time <- seq(5, 20, by = 0.1)

	# Create a typical protein SEC profile:
	# Small aggregate peak at t=7, main monomer at t=12, fragment at t=16
	signal <- 0.05 *
		dnorm(time, mean = 7, sd = 0.3) + # HMWS
		0.90 * dnorm(time, mean = 12, sd = 0.5) + # Monomer
		0.05 * dnorm(time, mean = 16, sd = 0.4) # LMWS

	test_data <- tibble::tibble(sample_id = "test")

	test_data$uv <- measure::new_measure_list(
		list(measure::new_measure_tbl(location = time, value = signal))
	)

	rec <- recipes::recipe(~., data = test_data) |>
		step_sec_aggregates(
			measures = "uv",
			monomer_start = 10,
			monomer_end = 14,
			method = "manual"
		)

	prepped <- recipes::prep(rec)
	result <- recipes::bake(prepped, new_data = NULL)

	expect_s3_class(result, "tbl_df")
	expect_true("purity_hmws" %in% names(result))
	expect_true("purity_monomer" %in% names(result))
	expect_true("purity_lmws" %in% names(result))

	# Monomer should be the largest fraction
	expect_gt(result$purity_monomer[1], result$purity_hmws[1])
	expect_gt(result$purity_monomer[1], result$purity_lmws[1])

	# Sum should be close to 100%
	total <- result$purity_hmws[1] +
		result$purity_monomer[1] +
		result$purity_lmws[1]
	expect_equal(total, 100, tolerance = 1)
})

test_that("step_sec_aggregates auto-detects main peak", {
	skip_if_not_installed("measure")

	time <- seq(5, 20, by = 0.1)
	signal <- dnorm(time, mean = 12, sd = 0.5)

	test_data <- tibble::tibble(sample_id = "test")

	test_data$uv <- measure::new_measure_list(
		list(measure::new_measure_tbl(location = time, value = signal))
	)

	rec <- recipes::recipe(~., data = test_data) |>
		step_sec_aggregates(measures = "uv", method = "tallest")

	prepped <- recipes::prep(rec)
	result <- recipes::bake(prepped, new_data = NULL)

	expect_s3_class(result, "tbl_df")
	expect_true("purity_main_start" %in% names(result))
	expect_true("purity_main_end" %in% names(result))

	# Main peak should be around t=12
	expect_gt(result$purity_main_start[1], 10)
	expect_lt(result$purity_main_end[1], 15)
})
