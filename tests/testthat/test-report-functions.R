# ==============================================================================
# Tests for report generation functions
# ==============================================================================

test_that("list_sec_templates returns available templates", {
  templates <- list_sec_templates()

  expect_s3_class(templates, "tbl_df")
  expect_true("template" %in% names(templates))
  expect_true("description" %in% names(templates))
  expect_true("formats" %in% names(templates))
  expect_equal(nrow(templates), 3)
  expect_true("standard" %in% templates$template)
  expect_true("detailed" %in% templates$template)
  expect_true("qc" %in% templates$template)
})

test_that("get_sec_template finds templates", {
  skip_on_cran()

  # Standard template
  path <- get_sec_template("standard")
  expect_true(file.exists(path))
  expect_true(grepl("\\.qmd$", path))

  # Detailed template
  path <- get_sec_template("detailed")
  expect_true(file.exists(path))

  # QC template
  path <- get_sec_template("qc")
  expect_true(file.exists(path))
})

test_that("get_sec_template validates template name", {
  expect_error(
    get_sec_template("nonexistent"),
    "should be one of"
  )
})

test_that("measure_sec_report errors without quarto package", {
  skip_if(rlang::is_installed("quarto"), "quarto is installed")

  test_data <- tibble::tibble(sample_id = "A")

  expect_error(
    measure_sec_report(test_data),
    "Package.*quarto.*is required"
  )
})

test_that("measure_sec_report validates data input", {
  skip_if_not_installed("quarto")
  skip_if(!nzchar(Sys.which("quarto")), "Quarto CLI not found")

  expect_error(
    measure_sec_report("not a data frame"),
    "must be a data frame"
  )
})

test_that("measure_sec_report validates template argument", {
  skip_if_not_installed("quarto")
  skip_if(!nzchar(Sys.which("quarto")), "Quarto CLI not found")

  test_data <- tibble::tibble(sample_id = "A")

  expect_error(
    measure_sec_report(test_data, template = "invalid"),
    "should be one of"
  )
})

test_that("measure_sec_report validates output_format argument", {
  skip_if_not_installed("quarto")
  skip_if(!nzchar(Sys.which("quarto")), "Quarto CLI not found")

  test_data <- tibble::tibble(sample_id = "A")

  expect_error(
    measure_sec_report(test_data, output_format = "rtf"),
    "should be one of"
  )
})


# ==============================================================================
# Integration tests (require quarto CLI)
# These tests are skipped by default as they require a working Quarto
# installation and can be slow. Run manually with:
#   Sys.setenv(RUN_QUARTO_TESTS = "true")
#   devtools::test(filter = "report")
# ==============================================================================

test_that("measure_sec_report generates HTML report", {
  skip_if_not_installed("quarto")
  skip_if_not_installed("measure")
  skip_on_cran()
  skip_on_ci()
  skip_if(!nzchar(Sys.which("quarto")), "Quarto CLI not found")
  skip_if(
    !identical(Sys.getenv("RUN_QUARTO_TESTS"), "true"),
    "Quarto integration tests disabled (set RUN_QUARTO_TESTS=true to enable)"
  )

  # Create minimal test data with measure columns
  time <- seq(5, 15, by = 0.1)
  signal <- dnorm(time, mean = 10, sd = 1)

  test_data <- tibble::tibble(
    sample_id = "Test Sample"
  )

  test_data$ri <- measure::new_measure_list(list(
    measure::new_measure_tbl(location = time, value = signal)
  ))

  # Add MW columns for summary
  test_data$Mn <- 50000
  test_data$Mw <- 100000
  test_data$Mz <- 150000
  test_data$dispersity <- 2.0

  # Generate report in temp directory
  output_file <- tempfile(fileext = ".html")

  withr::with_tempdir({
    result <- measure_sec_report(
      test_data,
      template = "standard",
      output_format = "html",
      output_file = output_file,
      title = "Test Report",
      sample_id = "sample_id",
      include_plots = FALSE,
      open = FALSE,
      quiet = TRUE
    )

    expect_true(file.exists(result))
    expect_true(grepl("\\.html$", result))

    # Check file has content
    content <- readLines(result, warn = FALSE)
    expect_true(length(content) > 0)
    expect_true(any(grepl("Test Report", content, fixed = TRUE)))
  })
})

test_that("measure_sec_report generates detailed report", {
  skip_if_not_installed("quarto")
  skip_if_not_installed("measure")
  skip_on_cran()
  skip_on_ci()
  skip_if(!nzchar(Sys.which("quarto")), "Quarto CLI not found")
  skip_if(
    !identical(Sys.getenv("RUN_QUARTO_TESTS"), "true"),
    "Quarto integration tests disabled"
  )

  time <- seq(5, 15, by = 0.1)
  signal <- dnorm(time, mean = 10, sd = 1)

  test_data <- tibble::tibble(sample_id = "Sample A")

  test_data$ri <- measure::new_measure_list(list(
    measure::new_measure_tbl(location = time, value = signal)
  ))

  test_data$Mn <- 50000
  test_data$Mw <- 100000

  output_file <- tempfile(fileext = ".html")

  withr::with_tempdir({
    result <- measure_sec_report(
      test_data,
      template = "detailed",
      output_format = "html",
      output_file = output_file,
      include_plots = FALSE,
      include_slice_table = FALSE,
      open = FALSE,
      quiet = TRUE
    )

    expect_true(file.exists(result))
  })
})

test_that("measure_sec_report generates QC report", {
  skip_if_not_installed("quarto")
  skip_if_not_installed("measure")
  skip_on_cran()
  skip_on_ci()
  skip_if(!nzchar(Sys.which("quarto")), "Quarto CLI not found")
  skip_if(
    !identical(Sys.getenv("RUN_QUARTO_TESTS"), "true"),
    "Quarto integration tests disabled"
  )

  time <- seq(5, 15, by = 0.1)
  signal <- dnorm(time, mean = 10, sd = 1)

  test_data <- tibble::tibble(sample_id = c("Std 1", "Std 2"))

  test_data$ri <- measure::new_measure_list(list(
    measure::new_measure_tbl(location = time, value = signal),
    measure::new_measure_tbl(location = time, value = signal * 0.9)
  ))

  output_file <- tempfile(fileext = ".html")

  withr::with_tempdir({
    result <- measure_sec_report(
      test_data,
      template = "qc",
      output_format = "html",
      output_file = output_file,
      specs = list(plate_count_min = 5000),
      include_plots = FALSE,
      open = FALSE,
      quiet = TRUE
    )

    expect_true(file.exists(result))
  })
})
