# ==============================================================================
# Tests for export functions
# ==============================================================================

test_that("measure_sec_slice_table extracts slice data", {
  skip_if_not_installed("measure")

  time <- seq(5, 15, by = 0.5)
  signal <- dnorm(time, mean = 10, sd = 0.5)

  test_data <- tibble::tibble(
    sample_id = c("sample1", "sample2")
  )

  test_data$ri <- measure::new_measure_list(list(
    measure::new_measure_tbl(location = time, value = signal),
    measure::new_measure_tbl(location = time, value = signal * 1.5)
  ))

  slices <- measure_sec_slice_table(test_data, measures = "ri", sample_id = "sample_id")

  expect_s3_class(slices, "tbl_df")
  expect_true("sample_id" %in% names(slices))
  expect_true("slice" %in% names(slices))
  expect_true("location" %in% names(slices))
  expect_true("measure" %in% names(slices))
  expect_true("value" %in% names(slices))

  # Should have correct number of rows
  expect_equal(nrow(slices), 2 * length(time))
})

test_that("measure_sec_slice_table pivots to wide format", {
  skip_if_not_installed("measure")

  time <- seq(5, 15, by = 0.5)
  signal1 <- dnorm(time, mean = 10, sd = 0.5)
  signal2 <- dnorm(time, mean = 10, sd = 0.6)

  test_data <- tibble::tibble(sample_id = "sample1")

  test_data$ri <- measure::new_measure_list(list(
    measure::new_measure_tbl(location = time, value = signal1)
  ))
  test_data$uv <- measure::new_measure_list(list(
    measure::new_measure_tbl(location = time, value = signal2)
  ))

  slices_wide <- measure_sec_slice_table(
    test_data,
    measures = c("ri", "uv"),
    pivot = TRUE
  )

  expect_s3_class(slices_wide, "tbl_df")
  expect_true("ri" %in% names(slices_wide))
  expect_true("uv" %in% names(slices_wide))
  expect_equal(nrow(slices_wide), length(time))
})

test_that("measure_sec_slice_table handles missing sample_id", {
  skip_if_not_installed("measure")

  time <- seq(5, 15, by = 0.5)
  signal <- dnorm(time, mean = 10, sd = 0.5)

  test_data <- tibble::tibble(x = 1)
  test_data$ri <- measure::new_measure_list(list(
    measure::new_measure_tbl(location = time, value = signal)
  ))

  slices <- measure_sec_slice_table(test_data, measures = "ri")

  # Should use row numbers as sample_id
  expect_equal(unique(slices$sample_id), 1)
})

test_that("measure_sec_summary_table creates summary", {
  test_data <- tibble::tibble(
    sample_id = c("A", "B", "C"),
    Mn = c(10000, 15000, 20000),
    Mw = c(15000, 22000, 30000),
    Mz = c(20000, 30000, 40000),
    dispersity = c(1.5, 1.47, 1.5),
    purity_monomer = c(98, 95, 92)
  )

  summary <- measure_sec_summary_table(test_data, sample_id = "sample_id")

  expect_s3_class(summary, "sec_summary_table")
  expect_true("Mn" %in% names(summary))
  expect_true("Mw" %in% names(summary))
  expect_true("purity_monomer" %in% names(summary))
  expect_equal(nrow(summary), 3)
})

test_that("measure_sec_summary_table handles missing columns gracefully", {
  test_data <- tibble::tibble(
    sample_id = c("A", "B"),
    Mn = c(10000, 15000)
  )

  summary <- measure_sec_summary_table(test_data, sample_id = "sample_id")

  expect_s3_class(summary, "sec_summary_table")
  expect_true("Mn" %in% names(summary))
  # Missing columns should not cause errors
  expect_false("purity_monomer" %in% names(summary))
})

test_that("measure_sec_summary_table rounds to specified digits", {
  test_data <- tibble::tibble(
    sample_id = "A",
    Mn = 12345.6789,
    Mw = 23456.7891
  )

  summary <- measure_sec_summary_table(test_data, sample_id = "sample_id", digits = 0)

  expect_equal(summary$Mn[1], 12346)
  expect_equal(summary$Mw[1], 23457)
})

test_that("measure_sec_summary_table print method works", {
  test_data <- tibble::tibble(
    sample_id = "A",
    Mn = 10000,
    Mw = 15000
  )

  summary <- measure_sec_summary_table(test_data, sample_id = "sample_id")

  expect_output(print(summary), "SEC Analysis Summary")
  expect_output(print(summary), "Samples: 1")
})
