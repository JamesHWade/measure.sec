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

  slices <- measure_sec_slice_table(
    test_data,
    measures = "ri",
    sample_id = "sample_id"
  )

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

  summary <- measure_sec_summary_table(
    test_data,
    sample_id = "sample_id",
    digits = 0
  )

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


# ==============================================================================
# measure_sec_compare tests
# ==============================================================================

test_that("measure_sec_compare compares MW averages", {
  sample1 <- tibble::tibble(
    sample_id = "A",
    mw_mn = 10000,
    mw_mw = 15000,
    mw_mz = 20000,
    mw_dispersity = 1.5
  )

  sample2 <- tibble::tibble(
    sample_id = "B",
    mw_mn = 12000,
    mw_mw = 18000,
    mw_mz = 24000,
    mw_dispersity = 1.5
  )

  result <- measure_sec_compare(
    sample1,
    sample2,
    samples = c("Sample 1", "Sample 2"),
    metrics = "mw_averages",
    plot = FALSE
  )

  expect_s3_class(result, "sec_comparison")
  expect_equal(nrow(result$summary), 2)
  expect_true("mw_mn" %in% names(result$summary))
  expect_true("mw_mw" %in% names(result$summary))
  expect_equal(result$samples, c("Sample 1", "Sample 2"))
  expect_equal(result$reference, "Sample 1")
})

test_that("measure_sec_compare calculates differences", {
  sample1 <- tibble::tibble(
    sample_id = "A",
    mw_mn = 10000,
    mw_mw = 20000
  )

  sample2 <- tibble::tibble(
    sample_id = "B",
    mw_mn = 11000,
    mw_mw = 22000
  )

  result <- measure_sec_compare(
    sample1,
    sample2,
    metrics = "mw_averages",
    plot = FALSE
  )

  # Check differences are calculated
  expect_true("mw_mn_diff" %in% names(result$differences))
  expect_true("mw_mn_pct" %in% names(result$differences))

  # Reference sample should have 0 difference
  expect_equal(result$differences$mw_mn_diff[1], 0)
  expect_equal(result$differences$mw_mn_pct[1], 0)

  # Second sample should have 10% increase
  expect_equal(result$differences$mw_mn_diff[2], 1000)
  expect_equal(result$differences$mw_mn_pct[2], 10)
})

test_that("measure_sec_compare accepts character reference", {
  sample1 <- tibble::tibble(mw_mn = 10000)
  sample2 <- tibble::tibble(mw_mn = 12000)

  result <- measure_sec_compare(
    sample1,
    sample2,
    samples = c("Control", "Test"),
    reference = "Test",
    metrics = "mw_averages",
    plot = FALSE
  )

  expect_equal(result$reference, "Test")
})

test_that("measure_sec_compare errors with fewer than 2 samples", {
  sample1 <- tibble::tibble(mw_mn = 10000)

  expect_error(
    measure_sec_compare(sample1, plot = FALSE),
    "At least 2 samples"
  )
})

test_that("measure_sec_compare errors with non-data frame input", {
  sample1 <- tibble::tibble(mw_mn = 10000)
  sample2 <- "not a data frame"

  expect_error(
    measure_sec_compare(sample1, sample2, plot = FALSE),
    "must be data frames"
  )
})

test_that("measure_sec_compare errors with wrong samples length", {
  sample1 <- tibble::tibble(mw_mn = 10000)
  sample2 <- tibble::tibble(mw_mn = 12000)

  expect_error(
    measure_sec_compare(
      sample1,
      sample2,
      samples = c("A", "B", "C"),
      plot = FALSE
    ),
    "must match number of data inputs"
  )
})

test_that("measure_sec_compare errors with invalid reference", {
  sample1 <- tibble::tibble(mw_mn = 10000)
  sample2 <- tibble::tibble(mw_mn = 12000)

  expect_error(
    measure_sec_compare(
      sample1,
      sample2,
      reference = "NonExistent",
      plot = FALSE
    ),
    "not found in sample names"
  )
})

test_that("measure_sec_compare handles branching metrics", {
  sample1 <- tibble::tibble(
    mw_mn = 10000,
    branching_index = 0.85,
    g_ratio = 0.9
  )

  sample2 <- tibble::tibble(
    mw_mn = 12000,
    branching_index = 0.75,
    g_ratio = 0.8
  )

  result <- measure_sec_compare(
    sample1,
    sample2,
    metrics = c("mw_averages", "branching"),
    plot = FALSE
  )

  expect_true("branching_index" %in% names(result$summary))
  expect_true("g_ratio" %in% names(result$summary))
})

test_that("measure_sec_compare print method works", {
  sample1 <- tibble::tibble(mw_mn = 10000, mw_mw = 15000)
  sample2 <- tibble::tibble(mw_mn = 12000, mw_mw = 18000)

  result <- measure_sec_compare(
    sample1,
    sample2,
    samples = c("Control", "Test"),
    metrics = "mw_averages",
    plot = FALSE
  )

  expect_output(print(result), "SEC Multi-Sample Comparison")
  expect_output(print(result), "Samples: 2")
  expect_output(print(result), "Reference: Control")
})

test_that("measure_sec_compare auto-generates sample names", {
  sample1 <- tibble::tibble(mw_mn = 10000)
  sample2 <- tibble::tibble(mw_mn = 12000)

  result <- measure_sec_compare(
    sample1,
    sample2,
    metrics = "mw_averages",
    plot = FALSE
  )

  # Should have auto-generated names
  expect_equal(length(result$samples), 2)
  expect_true(all(nchar(result$samples) > 0))
})

test_that("measure_sec_compare errors with numeric reference out of bounds", {
  sample1 <- tibble::tibble(mw_mn = 10000)
  sample2 <- tibble::tibble(mw_mn = 12000)

  expect_error(
    measure_sec_compare(sample1, sample2, reference = 5, plot = FALSE),
    "must be between 1 and"
  )

  expect_error(
    measure_sec_compare(sample1, sample2, reference = 0, plot = FALSE),
    "must be between 1 and"
  )
})

test_that("measure_sec_compare aggregates multi-row samples by mean", {
  sample1 <- tibble::tibble(
    mw_mn = c(10000, 10200, 9800),
    mw_mw = c(15000, 15300, 14700)
  )

  sample2 <- tibble::tibble(
    mw_mn = 12000,
    mw_mw = 18000
  )

  result <- measure_sec_compare(
    sample1,
    sample2,
    samples = c("Multi", "Single"),
    metrics = "mw_averages",
    plot = FALSE
  )

  # Should use mean of multi-row sample
  expect_equal(result$summary$mw_mn[1], 10000)
  expect_equal(result$summary$mw_mw[1], 15000)
})

test_that("measure_sec_compare warns about invalid metrics", {
  sample1 <- tibble::tibble(mw_mn = 10000)
  sample2 <- tibble::tibble(mw_mn = 12000)

  expect_warning(
    measure_sec_compare(
      sample1,
      sample2,
      metrics = c("mw_averages", "invalid_metric"),
      plot = FALSE
    ),
    "Unrecognized metrics ignored"
  )
})

test_that("measure_sec_compare errors with no valid metrics", {
  sample1 <- tibble::tibble(mw_mn = 10000)
  sample2 <- tibble::tibble(mw_mn = 12000)

  expect_error(
    measure_sec_compare(
      sample1,
      sample2,
      metrics = "not_a_real_metric",
      plot = FALSE
    ),
    "No valid metrics specified"
  )
})

test_that("measure_sec_compare validates plot parameter", {
  sample1 <- tibble::tibble(mw_mn = 10000)
  sample2 <- tibble::tibble(mw_mn = 12000)

  expect_error(
    measure_sec_compare(sample1, sample2, plot = "yes"),
    "must be TRUE or FALSE"
  )
})

test_that("measure_sec_compare validates digits parameter", {
  sample1 <- tibble::tibble(mw_mn = 10000)
  sample2 <- tibble::tibble(mw_mn = 12000)

  expect_error(
    measure_sec_compare(sample1, sample2, digits = -1, plot = FALSE),
    "must be a non-negative"
  )

  expect_error(
    measure_sec_compare(sample1, sample2, digits = "two", plot = FALSE),
    "must be a non-negative"
  )
})

test_that("measure_sec_compare works with 3+ samples", {
  samples <- lapply(1:5, function(i) {
    tibble::tibble(mw_mn = 10000 * i, mw_mw = 15000 * i)
  })

  result <- do.call(
    measure_sec_compare,
    c(
      samples,
      list(
        samples = paste("Batch", 1:5),
        metrics = "mw_averages",
        plot = FALSE
      )
    )
  )

  expect_equal(nrow(result$summary), 5)
  expect_equal(nrow(result$differences), 5)
})

test_that("measure_sec_compare respects digits parameter", {
  sample1 <- tibble::tibble(mw_mn = 10000.555)
  sample2 <- tibble::tibble(mw_mn = 12000.666)

  result_default <- measure_sec_compare(
    sample1,
    sample2,
    metrics = "mw_averages",
    plot = FALSE
  )

  result_0_digits <- measure_sec_compare(
    sample1,
    sample2,
    metrics = "mw_averages",
    plot = FALSE,
    digits = 0
  )

  expect_equal(result_default$summary$mw_mn[1], 10000.56)
  expect_equal(result_0_digits$summary$mw_mn[1], 10001)
})

test_that("measure_sec_compare warns when samples lack mw column for plot", {
  sample1 <- tibble::tibble(mw_mn = 10000)
  sample2 <- tibble::tibble(mw_mn = 12000)

  expect_warning(
    measure_sec_compare(
      sample1,
      sample2,
      samples = c("A", "B"),
      metrics = "mwd",
      plot = TRUE
    ),
    "excluded from plot"
  )
})

test_that("measure_sec_compare recognizes alternative column names", {
  sample1 <- tibble::tibble(Mn = 10000, Mw = 15000, dispersity = 1.5)
  sample2 <- tibble::tibble(Mn = 12000, Mw = 18000, dispersity = 1.5)

  result <- measure_sec_compare(
    sample1,
    sample2,
    metrics = "mw_averages",
    plot = FALSE
  )

  expect_true("Mn" %in% names(result$summary))
  expect_true("Mw" %in% names(result$summary))
})
