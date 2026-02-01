# ==============================================================================
# Tests for step_sec_exclude_regions
# ==============================================================================

# Helper to create test data with measure columns
create_test_sec_data <- function() {
	time <- seq(8, 20, by = 0.05)
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

test_that("step_sec_exclude_regions validates regions is a data frame", {
	skip_if_not_installed("measure")

	test_data <- create_test_sec_data()

	expect_error(
		recipes::recipe(~., data = test_data) |>
			step_sec_exclude_regions(regions = c(1, 2, 3)),
		"regions.*data frame"
	)
})

test_that("step_sec_exclude_regions validates regions has start and end", {
	skip_if_not_installed("measure")

	test_data <- create_test_sec_data()

	# Missing end
	expect_error(
		recipes::recipe(~., data = test_data) |>
			step_sec_exclude_regions(regions = tibble::tibble(start = 10)),
		"start.*end"
	)

	# Missing start
	expect_error(
		recipes::recipe(~., data = test_data) |>
			step_sec_exclude_regions(regions = tibble::tibble(end = 18)),
		"start.*end"
	)
})

test_that("step_sec_exclude_regions validates start and end are numeric", {
	skip_if_not_installed("measure")

	test_data <- create_test_sec_data()

	expect_error(
		recipes::recipe(~., data = test_data) |>
			step_sec_exclude_regions(
				regions = tibble::tibble(start = "10", end = "18")
			),
		"numeric"
	)
})

test_that("step_sec_exclude_regions validates start < end", {
	skip_if_not_installed("measure")

	test_data <- create_test_sec_data()

	expect_error(
		recipes::recipe(~., data = test_data) |>
			step_sec_exclude_regions(
				regions = tibble::tibble(start = 18, end = 10)
			),
		"start.*>=.*end"
	)

	expect_error(
		recipes::recipe(~., data = test_data) |>
			step_sec_exclude_regions(
				regions = tibble::tibble(start = 10, end = 10)
			),
		"start.*>=.*end"
	)
})

test_that("step_sec_exclude_regions validates purpose", {
	skip_if_not_installed("measure")

	test_data <- create_test_sec_data()
	regions <- tibble::tibble(start = 18, end = 20)

	expect_error(
		recipes::recipe(~., data = test_data) |>
			step_sec_exclude_regions(regions = regions, purpose = "invalid"),
		"purpose"
	)
})

# -- Prep and bake tests -------------------------------------------------------

test_that("step_sec_exclude_regions can prep with no regions", {
	skip_if_not_installed("measure")

	test_data <- create_test_sec_data()

	rec <- recipes::recipe(~., data = test_data) |>
		step_sec_exclude_regions()

	prepped <- recipes::prep(rec)
	expect_s3_class(prepped, "recipe")
})

test_that("step_sec_exclude_regions can prep with regions", {
	skip_if_not_installed("measure")

	test_data <- create_test_sec_data()
	regions <- tibble::tibble(start = 18.5, end = 20.0, reason = "Solvent peak")

	rec <- recipes::recipe(~., data = test_data) |>
		step_sec_exclude_regions(regions = regions)

	prepped <- recipes::prep(rec)
	expect_s3_class(prepped, "recipe")
})

test_that("step_sec_exclude_regions bake creates .excluded_regions column", {
	skip_if_not_installed("measure")

	test_data <- create_test_sec_data()
	regions <- tibble::tibble(start = 18.5, end = 20.0)

	rec <- recipes::recipe(~., data = test_data) |>
		step_sec_exclude_regions(regions = regions)

	result <- recipes::prep(rec) |>
		recipes::bake(new_data = NULL)

	expect_true(".excluded_regions" %in% names(result))
	expect_equal(length(result$.excluded_regions), 1)
})

test_that("step_sec_exclude_regions output has correct structure", {
	skip_if_not_installed("measure")

	test_data <- create_test_sec_data()
	regions <- tibble::tibble(
		start = 18.5,
		end = 20.0,
		reason = "Solvent peak"
	)

	rec <- recipes::recipe(~., data = test_data) |>
		step_sec_exclude_regions(regions = regions, purpose = "both")

	result <- recipes::prep(rec) |>
		recipes::bake(new_data = NULL)

	excl <- result$.excluded_regions[[1]]
	expect_s3_class(excl, "tbl_df")
	expect_true(all(c("start", "end", "purpose", "reason") %in% names(excl)))
	expect_equal(excl$start, 18.5)
	expect_equal(excl$end, 20.0)
	expect_equal(excl$purpose, "both")
	expect_equal(excl$reason, "Solvent peak")
})

test_that("step_sec_exclude_regions handles multiple regions", {
	skip_if_not_installed("measure")

	test_data <- create_test_sec_data()
	regions <- tibble::tibble(
		start = c(8.0, 18.5),
		end = c(9.0, 20.0),
		reason = c("Void volume", "Solvent peak")
	)

	rec <- recipes::recipe(~., data = test_data) |>
		step_sec_exclude_regions(regions = regions)

	result <- recipes::prep(rec) |>
		recipes::bake(new_data = NULL)

	excl <- result$.excluded_regions[[1]]
	expect_equal(nrow(excl), 2)
	expect_equal(excl$start, c(8.0, 18.5))
	expect_equal(excl$end, c(9.0, 20.0))
})

test_that("step_sec_exclude_regions respects purpose parameter", {
	skip_if_not_installed("measure")

	test_data <- create_test_sec_data()
	regions <- tibble::tibble(start = 18.5, end = 20.0)

	# Test each purpose
	for (p in c("baseline", "integration", "both")) {
		rec <- recipes::recipe(~., data = test_data) |>
			step_sec_exclude_regions(regions = regions, purpose = p)

		result <- recipes::prep(rec) |>
			recipes::bake(new_data = NULL)

		excl <- result$.excluded_regions[[1]]
		expect_equal(excl$purpose, p)
	}
})

test_that("step_sec_exclude_regions handles no regions gracefully", {
	skip_if_not_installed("measure")

	test_data <- create_test_sec_data()

	rec <- recipes::recipe(~., data = test_data) |>
		step_sec_exclude_regions(regions = NULL)

	result <- recipes::prep(rec) |>
		recipes::bake(new_data = NULL)

	excl <- result$.excluded_regions[[1]]
	expect_equal(nrow(excl), 0)
	expect_true(all(c("start", "end", "purpose", "reason") %in% names(excl)))
})

test_that("step_sec_exclude_regions works with multiple samples", {
	skip_if_not_installed("measure")

	test_data <- create_multi_sample_data()
	regions <- tibble::tibble(start = 18.5, end = 20.0)

	rec <- recipes::recipe(~., data = test_data) |>
		step_sec_exclude_regions(regions = regions)

	result <- recipes::prep(rec) |>
		recipes::bake(new_data = NULL)

	expect_equal(nrow(result), 2)
	expect_equal(length(result$.excluded_regions), 2)

	# Both samples should have the same exclusion
	expect_equal(result$.excluded_regions[[1]]$start, 18.5)
	expect_equal(result$.excluded_regions[[2]]$start, 18.5)
})

test_that("step_sec_exclude_regions handles sample-specific exclusions", {
	skip_if_not_installed("measure")

	test_data <- create_multi_sample_data()

	# Global exclusion + sample-specific exclusion
	regions <- tibble::tibble(
		start = c(18.5, 15.0),
		end = c(20.0, 16.0),
		sample_id = c(NA, "sample2"), # NA = global
		reason = c("Solvent", "Artifact")
	)

	rec <- recipes::recipe(~., data = test_data) |>
		step_sec_exclude_regions(regions = regions)

	result <- recipes::prep(rec) |>
		recipes::bake(new_data = NULL)

	# sample1 should only have the global exclusion
	expect_equal(nrow(result$.excluded_regions[[1]]), 1)
	expect_equal(result$.excluded_regions[[1]]$start, 18.5)

	# sample2 should have both
	expect_equal(nrow(result$.excluded_regions[[2]]), 2)
})

# -- Tidy method tests ---------------------------------------------------------

test_that("tidy.step_sec_exclude_regions returns expected tibble", {
	skip_if_not_installed("measure")

	test_data <- create_test_sec_data()
	regions <- tibble::tibble(
		start = c(8.0, 18.5),
		end = c(9.0, 20.0)
	)

	rec <- recipes::recipe(~., data = test_data) |>
		step_sec_exclude_regions(regions = regions, purpose = "integration")

	# Before prep
	tidy_result <- recipes::tidy(rec, number = 1)
	expect_s3_class(tidy_result, "tbl_df")
	expect_true("n_regions" %in% names(tidy_result))
	expect_true("purpose" %in% names(tidy_result))
	expect_equal(tidy_result$n_regions, 2)
	expect_equal(tidy_result$purpose, "integration")

	# After prep
	prepped <- recipes::prep(rec)
	tidy_prepped <- recipes::tidy(prepped, number = 1)
	expect_s3_class(tidy_prepped, "tbl_df")
	expect_equal(tidy_prepped$n_regions, 2)
})

# -- Print method tests --------------------------------------------------------

test_that("print.step_sec_exclude_regions works", {
	skip_if_not_installed("measure")

	test_data <- create_test_sec_data()

	# With regions
	regions <- tibble::tibble(start = c(8.0, 18.5), end = c(9.0, 20.0))
	rec <- recipes::recipe(~., data = test_data) |>
		step_sec_exclude_regions(regions = regions, purpose = "baseline")

	expect_output(print(rec), "2 region")
	expect_output(print(rec), "baseline")

	# Without regions
	rec_empty <- recipes::recipe(~., data = test_data) |>
		step_sec_exclude_regions()

	expect_output(print(rec_empty), "no regions")
})

# -- required_pkgs tests -------------------------------------------------------

test_that("required_pkgs.step_sec_exclude_regions lists measure packages", {
	skip_if_not_installed("measure")

	test_data <- create_test_sec_data()

	rec <- recipes::recipe(~., data = test_data) |>
		step_sec_exclude_regions()

	pkgs <- recipes::required_pkgs(rec)
	expect_true("measure.sec" %in% pkgs)
	expect_true("measure" %in% pkgs)
})
