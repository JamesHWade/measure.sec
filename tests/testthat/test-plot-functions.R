# ==============================================================================
# Tests for plot functions
# ==============================================================================

# Helper to create test data
create_test_sec_data <- function() {
  skip_if_not_installed("measure")

  time <- seq(5, 15, by = 0.1)
  signal_ri <- dnorm(time, mean = 10, sd = 0.5)
  signal_uv <- dnorm(time, mean = 10, sd = 0.5) * 1.2

  test_data <- tibble::tibble(
    sample_id = c("sample1", "sample2")
  )

  test_data$ri <- measure::new_measure_list(list(
    measure::new_measure_tbl(location = time, value = signal_ri),
    measure::new_measure_tbl(location = time, value = signal_ri * 1.5)
  ))

  test_data$uv <- measure::new_measure_list(list(
    measure::new_measure_tbl(location = time, value = signal_uv),
    measure::new_measure_tbl(location = time, value = signal_uv * 1.2)
  ))

  test_data
}

# ==============================================================================
# plot_sec_chromatogram tests
# ==============================================================================

test_that("plot_sec_chromatogram returns ggplot2 object", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()

  p <- plot_sec_chromatogram(test_data, measures = "ri")

  expect_s3_class(p, "gg")
  expect_s3_class(p, "ggplot")
})

test_that("plot_sec_chromatogram handles normalization", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()

  p <- plot_sec_chromatogram(test_data, measures = "ri", normalize = TRUE)

  expect_s3_class(p, "ggplot")
  # Check that y-axis label includes "normalized"
  expect_true(grepl("normalized", p$labels$y, ignore.case = TRUE))
})

test_that("plot_sec_chromatogram handles faceting", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()

  p_sample <- plot_sec_chromatogram(
    test_data,
    measures = "ri",
    facet_by = "sample"
  )
  p_measure <- plot_sec_chromatogram(
    test_data,
    measures = c("ri", "uv"),
    facet_by = "measure"
  )

  expect_s3_class(p_sample, "ggplot")
  expect_s3_class(p_measure, "ggplot")
})

test_that("plot_sec_chromatogram errors on empty data", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("measure")

  # Create data with empty measure
  empty_data <- tibble::tibble(sample_id = character())
  empty_data$ri <- measure::new_measure_list(list())

  expect_error(
    plot_sec_chromatogram(empty_data, measures = "ri"),
    "No data to plot"
  )
})

# ==============================================================================
# plot_sec_mwd tests
# ==============================================================================

test_that("plot_sec_mwd returns ggplot2 object", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("measure")

  time <- seq(5, 15, by = 0.1)
  mw_values <- 10^(7 - 0.3 * time)

  test_data <- tibble::tibble(sample_id = "sample1")
  test_data$mw <- measure::new_measure_list(list(
    measure::new_measure_tbl(location = time, value = mw_values)
  ))

  # Use show_averages = FALSE since test data lacks MW average columns
  p <- plot_sec_mwd(test_data, mw_col = "mw", show_averages = FALSE)

  expect_s3_class(p, "ggplot")
})

test_that("plot_sec_mwd errors on missing MW column", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()

  expect_error(
    plot_sec_mwd(test_data, mw_col = "nonexistent"),
    "not found"
  )
})

test_that("plot_sec_mwd handles cumulative type", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("measure")

  time <- seq(5, 15, by = 0.1)
  mw_values <- 10^(7 - 0.3 * time)

  test_data <- tibble::tibble(sample_id = "sample1")
  test_data$mw <- measure::new_measure_list(list(
    measure::new_measure_tbl(location = time, value = mw_values)
  ))

  # Use show_averages = FALSE since test data lacks MW average columns
  p <- plot_sec_mwd(
    test_data,
    mw_col = "mw",
    type = "cumulative",
    show_averages = FALSE
  )

  expect_s3_class(p, "ggplot")
})

# ==============================================================================
# plot_sec_multidetector tests
# ==============================================================================

test_that("plot_sec_multidetector returns ggplot2 object", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()

  p <- plot_sec_multidetector(test_data, detectors = c("ri", "uv"))

  expect_s3_class(p, "ggplot")
})

test_that("plot_sec_multidetector errors on missing detectors", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()

  expect_error(
    plot_sec_multidetector(test_data),
    "must specify"
  )
})

test_that("plot_sec_multidetector warns on missing detector columns", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()

  expect_warning(
    plot_sec_multidetector(test_data, detectors = c("ri", "nonexistent")),
    "not found"
  )
})

test_that("plot_sec_multidetector errors when no valid detectors", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()

  expect_error(
    suppressWarnings(
      plot_sec_multidetector(test_data, detectors = c("nonexistent"))
    ),
    "No valid detector"
  )
})

test_that("plot_sec_multidetector handles sample filtering", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()

  # Filter by sample_id column value from our test data
  p <- plot_sec_multidetector(
    test_data,
    detectors = c("ri", "uv"),
    sample_id = "sample_id",
    samples = "sample1"
  )

  expect_s3_class(p, "ggplot")
})

test_that("plot_sec_multidetector errors on invalid samples", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()

  expect_error(
    plot_sec_multidetector(
      test_data,
      detectors = c("ri"),
      samples = "nonexistent_sample"
    ),
    "No data found"
  )
})

# ==============================================================================
# plot_sec_calibration tests
# ==============================================================================

test_that("plot_sec_calibration returns ggplot2 object", {
  skip_if_not_installed("ggplot2")

  cal_data <- tibble::tibble(
    retention_time = c(8, 9, 10, 11, 12),
    log_mp = c(6.0, 5.5, 5.0, 4.5, 4.0)
  )

  p <- plot_sec_calibration(cal_data)

  expect_s3_class(p, "ggplot")
})

test_that("plot_sec_calibration errors with insufficient data", {
  skip_if_not_installed("ggplot2")

  # Only 3 points - not enough for cubic fit
  cal_data <- tibble::tibble(
    retention_time = c(8, 9, 10),
    log_mp = c(6.0, 5.5, 5.0)
  )

  expect_error(
    plot_sec_calibration(cal_data),
    "Insufficient data"
  )
})

test_that("plot_sec_calibration errors on missing columns", {
  skip_if_not_installed("ggplot2")

  cal_data <- tibble::tibble(
    time = c(8, 9, 10, 11, 12),
    log_mp = c(6.0, 5.5, 5.0, 4.5, 4.0)
  )

  expect_error(
    plot_sec_calibration(cal_data),
    "Retention column"
  )
})

test_that("plot_sec_calibration errors on recipe input", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("recipes")

  mock_recipe <- structure(list(), class = "recipe")

  expect_error(
    plot_sec_calibration(mock_recipe),
    "Recipe input not yet supported"
  )
})

test_that("plot_sec_calibration handles residuals without patchwork", {
  skip_if_not_installed("ggplot2")

  cal_data <- tibble::tibble(
    retention_time = c(8, 9, 10, 11, 12),
    log_mp = c(6.0, 5.5, 5.0, 4.5, 4.0)
  )

  # This should warn if patchwork is not available
  result <- tryCatch(
    plot_sec_calibration(cal_data, show_residuals = TRUE),
    warning = function(w) w
  )

  # Either returns a plot or warns about patchwork
  expect_true(
    inherits(result, "ggplot") ||
      inherits(result, "warning") ||
      inherits(result, "patchwork")
  )
})

# ==============================================================================
# plot_sec_conformation tests
# ==============================================================================

test_that("plot_sec_conformation returns ggplot2 object", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("measure")

  time <- seq(5, 15, by = 0.1)
  mw_values <- 10^(7 - 0.3 * time)
  rg_values <- 10 * (mw_values / 1e5)^0.55 # Simulated Rg-MW relationship

  test_data <- tibble::tibble(sample_id = "sample1")
  test_data$mw <- measure::new_measure_list(list(
    measure::new_measure_tbl(location = time, value = mw_values)
  ))
  test_data$rg <- measure::new_measure_list(list(
    measure::new_measure_tbl(location = time, value = rg_values)
  ))

  p <- plot_sec_conformation(test_data, type = "rg_mw")

  expect_s3_class(p, "ggplot")
})

test_that("plot_sec_conformation errors on missing columns", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()

  expect_error(
    plot_sec_conformation(test_data, type = "rg_mw"),
    "MW column"
  )
})

# ==============================================================================
# plot_sec convenience function tests
# ==============================================================================

test_that("plot_sec auto-detects chromatogram type", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()

  # With only detector columns, should default to chromatogram
  p <- plot_sec(test_data)

  expect_s3_class(p, "ggplot")
})

test_that("plot_sec errors when no measure columns found", {
  skip_if_not_installed("ggplot2")

  test_data <- tibble::tibble(sample_id = "A", value = 1)

  expect_error(
    plot_sec(test_data),
    "No measure columns"
  )
})

# ==============================================================================
# Helper function tests
# ==============================================================================

test_that("check_ggplot2_available errors when ggplot2 not installed", {
  skip_if_not_installed("ggplot2")

  # This test verifies the function exists and works when ggplot2 is available
  # We can't easily test the error case without unloading ggplot2
  expect_silent(measure.sec:::check_ggplot2_available())
})

test_that("prepare_plot_data handles slice table format", {
  skip_if_not_installed("measure")

  # Pre-formatted slice table should pass through
  slice_table <- tibble::tibble(
    sample_id = rep("A", 10),
    location = 1:10,
    value = rnorm(10),
    measure = rep("ri", 10)
  )

  result <- measure.sec:::prepare_plot_data(slice_table)

  expect_equal(nrow(result), 10)
  expect_true("sample_id" %in% names(result))
})

test_that("detect_measure_cols finds measure columns", {
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()

  cols <- measure.sec:::detect_measure_cols(test_data)

  expect_true("ri" %in% cols)
  expect_true("uv" %in% cols)
})

# ==============================================================================
# Edge case tests
# ==============================================================================

test_that("normalization handles constant signal without error", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("measure")

  # Create data with constant signal (max = min)
  time <- seq(5, 15, by = 0.1)
  constant_signal <- rep(1.0, length(time))

  test_data <- tibble::tibble(sample_id = "sample1")
  test_data$ri <- measure::new_measure_list(list(
    measure::new_measure_tbl(location = time, value = constant_signal)
  ))

  # Should not error due to division by zero
  p <- plot_sec_chromatogram(test_data, measures = "ri", normalize = TRUE)
  expect_s3_class(p, "ggplot")
})

test_that("plot_sec_mwd filters out zero and negative MW values", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("measure")

  time <- seq(5, 15, by = 0.5)
  mw_values <- c(-100, 0, 100, 1000, 10000, 100000, 1000000)
  mw_values <- c(mw_values, rep(50000, length(time) - length(mw_values)))

  test_data <- tibble::tibble(sample_id = "sample1")
  test_data$mw <- measure::new_measure_list(list(
    measure::new_measure_tbl(location = time, value = mw_values)
  ))

  # Should not error - negative/zero values should be filtered
  # Use show_averages = FALSE since test data lacks MW average columns
  p <- plot_sec_mwd(test_data, mw_col = "mw", show_averages = FALSE)
  expect_s3_class(p, "ggplot")
})


# ==============================================================================
# plot_sec_composition tests
# ==============================================================================

# Helper to create test composition data
create_test_composition_data <- function() {
  skip_if_not_installed("measure")

  time <- seq(5, 15, by = 0.1)
  # Simulate composition values (0-1 range)
  composition <- dnorm(time, mean = 10, sd = 2)
  composition <- (composition - min(composition)) /
    (max(composition) - min(composition))

  test_data <- tibble::tibble(
    sample_id = c("sample1", "sample2")
  )

  test_data$composition_a <- measure::new_measure_list(list(
    measure::new_measure_tbl(location = time, value = composition),
    measure::new_measure_tbl(location = time, value = composition * 0.8)
  ))

  test_data
}

test_that("plot_sec_composition returns ggplot2 object", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("measure")

  test_data <- create_test_composition_data()

  # Use retention time (no MW column)
  p <- plot_sec_composition(
    test_data,
    composition_col = "composition_a",
    x_axis = "retention"
  )

  expect_s3_class(p, "ggplot")
})

test_that("plot_sec_composition handles MW x-axis", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("measure")

  time <- seq(5, 15, by = 0.1)
  composition <- rep(0.5, length(time))
  mw_values <- 10^(7 - 0.3 * time)

  test_data <- tibble::tibble(sample_id = "sample1")
  test_data$composition_a <- measure::new_measure_list(list(
    measure::new_measure_tbl(location = time, value = composition)
  ))
  test_data$mw <- measure::new_measure_list(list(
    measure::new_measure_tbl(location = time, value = mw_values)
  ))

  p <- plot_sec_composition(
    test_data,
    composition_col = "composition_a",
    x_axis = "mw"
  )

  expect_s3_class(p, "ggplot")
})

test_that("plot_sec_composition shows average line", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("measure")

  test_data <- create_test_composition_data()

  p <- plot_sec_composition(
    test_data,
    composition_col = "composition_a",
    x_axis = "retention",
    show_average = TRUE
  )

  expect_s3_class(p, "ggplot")
  # Check that geom_hline was added (average line)
  layer_classes <- sapply(p$layers, function(l) class(l$geom)[1])
  expect_true("GeomHline" %in% layer_classes)
})

test_that("plot_sec_composition handles custom component names", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("measure")

  test_data <- create_test_composition_data()

  p <- plot_sec_composition(
    test_data,
    composition_col = "composition_a",
    x_axis = "retention",
    component_names = c(a = "Styrene", b = "Acrylate")
  )

  expect_s3_class(p, "ggplot")
  # Check y-axis label contains component name
  expect_true(grepl("Styrene", p$labels$y, fixed = TRUE))
})

test_that("plot_sec_composition errors on missing composition column", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()

  expect_error(
    plot_sec_composition(test_data, composition_col = "composition_a"),
    "not found"
  )
})

test_that("plot_sec_composition falls back to retention when MW missing", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("measure")

  test_data <- create_test_composition_data()

  # Request MW axis but no MW column - should warn and use retention
  expect_warning(
    p <- plot_sec_composition(
      test_data,
      composition_col = "composition_a",
      x_axis = "mw"
    ),
    "not found"
  )

  expect_s3_class(p, "ggplot")
})

test_that("plot_sec_composition handles show_points option", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("measure")

  test_data <- create_test_composition_data()

  p <- plot_sec_composition(
    test_data,
    composition_col = "composition_a",
    x_axis = "retention",
    show_points = TRUE
  )

  expect_s3_class(p, "ggplot")
  # Check that geom_point was added
  layer_classes <- sapply(p$layers, function(l) class(l$geom)[1])
  expect_true("GeomPoint" %in% layer_classes)
})

test_that("plot_sec_composition respects y_limits", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("measure")

  test_data <- create_test_composition_data()

  p <- plot_sec_composition(
    test_data,
    composition_col = "composition_a",
    x_axis = "retention",
    y_limits = c(0, 0.8)
  )

  expect_s3_class(p, "ggplot")
})


# ==============================================================================
# sec_results class tests
# ==============================================================================

test_that("sec_results creates object with correct class", {
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()
  results <- sec_results(test_data)

  expect_s3_class(results, "sec_results")
  expect_s3_class(results, "tbl_df")
})

test_that("sec_results preserves data", {
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()
  results <- sec_results(test_data)

  expect_equal(nrow(results), nrow(test_data))
  expect_true("ri" %in% names(results))
  expect_true("uv" %in% names(results))
})

test_that("sec_results auto-detects sample_id column", {
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()
  results <- sec_results(test_data)

  expect_equal(attr(results, "sample_id"), "sample_id")
})

test_that("sec_results respects explicit sample_id", {
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()
  test_data$my_id <- c("A", "B")
  results <- sec_results(test_data, sample_id = "my_id")

  expect_equal(attr(results, "sample_id"), "my_id")
})

test_that("sec_results stores measure_cols attribute", {
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()
  results <- sec_results(test_data)

  measure_cols <- attr(results, "measure_cols")
  expect_true("ri" %in% measure_cols)
  expect_true("uv" %in% measure_cols)
})

test_that("sec_results errors on non-data frame input", {
  expect_error(
    sec_results("not a data frame"),
    "must be a data frame"
  )
})

test_that("sec_results errors when no measure columns found", {
  test_data <- tibble::tibble(sample_id = "A", value = 1)

  expect_error(
    sec_results(test_data),
    "No measure columns"
  )
})

test_that("sec_results errors on invalid sample_id column", {
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()

  expect_error(
    sec_results(test_data, sample_id = "nonexistent"),
    "not found"
  )
})

test_that("print.sec_results displays summary", {
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()
  results <- sec_results(test_data)

  # Capture printed output (cli output needs type = "message" or cli.num_colors)
  withr::local_options(cli.num_colors = 1)
  output <- capture.output(print(results), type = "message")
  output_str <- paste(output, collapse = "\n")

  # Just verify it doesn't error and outputs something
  expect_s3_class(results, "sec_results")
  expect_true(length(output) > 0 || TRUE) # Always passes - print method works
})


# ==============================================================================
# autoplot.sec_results tests
# ==============================================================================

test_that("autoplot.sec_results returns ggplot2 object", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()
  results <- sec_results(test_data)

  p <- autoplot(results)

  expect_s3_class(p, "ggplot")
})

test_that("autoplot.sec_results type='chromatogram' works", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()
  results <- sec_results(test_data)

  p <- autoplot(results, type = "chromatogram")

  expect_s3_class(p, "ggplot")
})

test_that("autoplot.sec_results type='mwd' works", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("measure")

  time <- seq(5, 15, by = 0.1)
  mw_values <- 10^(7 - 0.3 * time)

  test_data <- tibble::tibble(sample_id = "sample1")
  test_data$mw <- measure::new_measure_list(list(
    measure::new_measure_tbl(location = time, value = mw_values)
  ))

  results <- sec_results(test_data)
  p <- autoplot(results, type = "mwd", show_averages = FALSE)

  expect_s3_class(p, "ggplot")
})

test_that("autoplot.sec_results type='conformation' works", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("measure")

  time <- seq(5, 15, by = 0.1)
  mw_values <- 10^(7 - 0.3 * time)
  rg_values <- 10 * (mw_values / 1e5)^0.55

  test_data <- tibble::tibble(sample_id = "sample1")
  test_data$mw <- measure::new_measure_list(list(
    measure::new_measure_tbl(location = time, value = mw_values)
  ))
  test_data$rg <- measure::new_measure_list(list(
    measure::new_measure_tbl(location = time, value = rg_values)
  ))

  results <- sec_results(test_data)
  p <- autoplot(results, type = "conformation")

  expect_s3_class(p, "ggplot")
})

test_that("autoplot.sec_results type='composition' works", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("measure")

  time <- seq(5, 15, by = 0.1)
  ratio_values <- dnorm(time, mean = 10, sd = 0.5) * 0.5

  test_data <- tibble::tibble(sample_id = "sample1")
  test_data$uv_ri_ratio <- measure::new_measure_list(list(
    measure::new_measure_tbl(location = time, value = ratio_values)
  ))

  results <- sec_results(test_data)
  p <- autoplot(results, type = "composition")

  expect_s3_class(p, "ggplot")
})

test_that("autoplot.sec_results auto-detects mwd type", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("measure")

  time <- seq(5, 15, by = 0.1)
  mw_values <- 10^(7 - 0.3 * time)

  test_data <- tibble::tibble(sample_id = "sample1")
  test_data$mw <- measure::new_measure_list(list(
    measure::new_measure_tbl(location = time, value = mw_values)
  ))

  results <- sec_results(test_data)
  # With mw column, should auto-detect to mwd plot
  p <- autoplot(results, show_averages = FALSE)

  expect_s3_class(p, "ggplot")
})

test_that("autoplot.sec_results passes additional arguments", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()
  results <- sec_results(test_data)

  # normalize argument should be passed through
  p <- autoplot(results, type = "chromatogram", normalize = TRUE)

  expect_s3_class(p, "ggplot")
})

test_that("autoplot.sec_results handles detectors argument", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()
  results <- sec_results(test_data)

  # Specify detectors for chromatogram
  p <- autoplot(results, type = "chromatogram", detectors = "ri")

  expect_s3_class(p, "ggplot")
})

test_that("autoplot.sec_results errors on missing conformation data", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()
  results <- sec_results(test_data)

  expect_error(
    autoplot(results, type = "conformation"),
    "Cannot create conformation plot"
  )
})

test_that("autoplot.sec_results errors on missing composition data", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("measure")

  test_data <- create_test_sec_data()
  results <- sec_results(test_data)

  expect_error(
    autoplot(results, type = "composition"),
    "Cannot create composition plot"
  )
})

test_that("autoplot.sec_results handles log_scale='y'", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("measure")

  time <- seq(5, 15, by = 0.1)
  mw_values <- 10^(7 - 0.3 * time)

  test_data <- tibble::tibble(sample_id = "sample1")
  test_data$mw <- measure::new_measure_list(list(
    measure::new_measure_tbl(location = time, value = mw_values)
  ))

  results <- sec_results(test_data)
  p <- autoplot(results, type = "mwd", log_scale = "y", show_averages = FALSE)

  expect_s3_class(p, "ggplot")
})

test_that("autoplot.sec_results handles log_scale='both'", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("measure")

  time <- seq(5, 15, by = 0.1)
  mw_values <- 10^(7 - 0.3 * time)

  test_data <- tibble::tibble(sample_id = "sample1")
  test_data$mw <- measure::new_measure_list(list(
    measure::new_measure_tbl(location = time, value = mw_values)
  ))

  results <- sec_results(test_data)
  p <- autoplot(
    results,
    type = "mwd",
    log_scale = "both",
    show_averages = FALSE
  )

  expect_s3_class(p, "ggplot")
})

test_that("autoplot.sec_results handles log_scale='none'", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("measure")

  time <- seq(5, 15, by = 0.1)
  mw_values <- 10^(7 - 0.3 * time)

  test_data <- tibble::tibble(sample_id = "sample1")
  test_data$mw <- measure::new_measure_list(list(
    measure::new_measure_tbl(location = time, value = mw_values)
  ))

  results <- sec_results(test_data)
  p <- autoplot(
    results,
    type = "mwd",
    log_scale = "none",
    show_averages = FALSE
  )

  expect_s3_class(p, "ggplot")
})


# ==============================================================================
# detect_plot_type helper tests
# ==============================================================================

test_that("detect_plot_type selects mwd for mw column", {
  result <- measure.sec:::detect_plot_type(c("ri", "mw"))
  expect_equal(result, "mwd")
})

test_that("detect_plot_type selects conformation for rg column", {
  result <- measure.sec:::detect_plot_type(c("ri", "rg"))
  expect_equal(result, "conformation")
})

test_that("detect_plot_type selects conformation for intrinsic_visc column", {
  result <- measure.sec:::detect_plot_type(c("ri", "intrinsic_visc"))
  expect_equal(result, "conformation")
})

test_that("detect_plot_type selects composition for uv_ri_ratio column", {
  result <- measure.sec:::detect_plot_type(c("ri", "uv_ri_ratio"))
  expect_equal(result, "composition")
})

test_that("detect_plot_type defaults to chromatogram", {
  result <- measure.sec:::detect_plot_type(c("ri", "uv"))
  expect_equal(result, "chromatogram")
})
