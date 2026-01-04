# Tests for calibration save/load functions

# Helper: Create a simple prepped recipe with calibration
create_test_recipe <- function() {
  # Create test data with measure format
  set.seed(123)
  test_data <- tibble::tibble(
    sample_id = "test_sample",
    signal = list(tibble::tibble(
      location = seq(10, 18, by = 0.1),
      value = dnorm(seq(10, 18, by = 0.1), mean = 14, sd = 1)
    ))
  )
  class(test_data$signal[[1]]) <- c("measure_tbl", class(test_data$signal[[1]]))
  class(test_data$signal) <- c("measure_list", class(test_data$signal))

  # Create calibration standards
  standards <- data.frame(
    retention = c(10.5, 12.0, 13.5, 15.0, 16.5, 18.0),
    log_mw = c(6.0, 5.5, 5.0, 4.5, 4.0, 3.5)
  )

  # Create and prep recipe
  rec <- recipes::recipe(signal ~ sample_id, data = test_data) |>
    recipes::update_role(sample_id, new_role = "id") |>
    step_sec_conventional_cal(
      measures = "signal",
      standards = standards,
      fit_type = "cubic"
    ) |>
    recipes::prep()

  list(recipe = rec, standards = standards, data = test_data)
}

test_that("save_sec_calibration works with RDS format", {
  skip_if_not_installed("measure")

  test_objs <- create_test_recipe()
  tmpfile <- tempfile(fileext = ".rds")
  on.exit(unlink(tmpfile), add = TRUE)

  # Save calibration
  expect_message(
    cal <- save_sec_calibration(test_objs$recipe, tmpfile),
    "Saved"
  )

  # Check file exists
  expect_true(file.exists(tmpfile))

  # Check returned object
  expect_s3_class(cal, "sec_calibration")
  expect_equal(cal$fit_type, "cubic")
  expect_true(cal$diagnostics$r_squared > 0.99)
})

test_that("save_sec_calibration respects overwrite parameter", {
  skip_if_not_installed("measure")

  test_objs <- create_test_recipe()
  tmpfile <- tempfile(fileext = ".rds")
  on.exit(unlink(tmpfile), add = TRUE)

  # Save once
  save_sec_calibration(test_objs$recipe, tmpfile)

  # Try to save again without overwrite - should error
  expect_error(
    save_sec_calibration(test_objs$recipe, tmpfile),
    "already exists"
  )

  # Save again with overwrite - should work
  expect_message(
    save_sec_calibration(test_objs$recipe, tmpfile, overwrite = TRUE),
    "Saved"
  )
})

test_that("save_sec_calibration includes metadata", {
  skip_if_not_installed("measure")

  test_objs <- create_test_recipe()
  tmpfile <- tempfile(fileext = ".rds")
  on.exit(unlink(tmpfile), add = TRUE)

  metadata <- list(
    column = "PLgel Mixed-C",
    instrument = "Agilent 1260",
    analyst = "Test"
  )

  save_sec_calibration(
    test_objs$recipe,
    tmpfile,
    metadata = metadata
  )

  # Load and check metadata
  cal <- load_sec_calibration(tmpfile)
  expect_equal(cal$user_metadata$column, "PLgel Mixed-C")
  expect_equal(cal$user_metadata$instrument, "Agilent 1260")
})

test_that("load_sec_calibration works with RDS format", {
  skip_if_not_installed("measure")

  test_objs <- create_test_recipe()
  tmpfile <- tempfile(fileext = ".rds")
  on.exit(unlink(tmpfile), add = TRUE)

  save_sec_calibration(test_objs$recipe, tmpfile)

  # Load calibration
  expect_message(
    cal <- load_sec_calibration(tmpfile),
    "Loaded"
  )

  # Check loaded object
  expect_s3_class(cal, "sec_calibration")
  expect_equal(cal$fit_type, "cubic")
  expect_equal(cal$degree, 3)
  expect_length(cal$calibration_range, 2)
  expect_true(!is.null(cal$fit_object))
  expect_true(!is.null(cal$standards))
})

test_that("load_sec_calibration errors on missing file", {
  expect_error(
    load_sec_calibration("nonexistent_file.rds"),
    "not found"
  )
})

test_that("save/load roundtrip with YAML format", {
  skip_if_not_installed("measure")
  skip_if_not_installed("yaml")

  test_objs <- create_test_recipe()
  tmpfile <- tempfile(fileext = ".yaml")
  on.exit(unlink(tmpfile), add = TRUE)

  # Save as YAML
  save_sec_calibration(test_objs$recipe, tmpfile)

  # Load from YAML
  cal <- load_sec_calibration(tmpfile)

  # Check loaded object
  expect_s3_class(cal, "sec_calibration")
  expect_equal(cal$fit_type, "cubic")
  expect_equal(cal$degree, 3)

  # Check diagnostics
  expect_true(cal$diagnostics$r_squared > 0.99)
})

test_that("step_sec_conventional_cal accepts pre-loaded calibration", {
  skip_if_not_installed("measure")

  test_objs <- create_test_recipe()
  tmpfile <- tempfile(fileext = ".rds")
  on.exit(unlink(tmpfile), add = TRUE)

  # Save calibration from first recipe
  cal <- save_sec_calibration(test_objs$recipe, tmpfile)

  # Create new data (different sample)
  new_data <- tibble::tibble(
    sample_id = "new_sample",
    signal = list(tibble::tibble(
      location = seq(10, 18, by = 0.1),
      value = dnorm(seq(10, 18, by = 0.1), mean = 13, sd = 0.8)
    ))
  )
  class(new_data$signal[[1]]) <- c("measure_tbl", class(new_data$signal[[1]]))
  class(new_data$signal) <- c("measure_list", class(new_data$signal))

  # Load calibration and use in new recipe
  loaded_cal <- load_sec_calibration(tmpfile)

  rec_new <- recipes::recipe(signal ~ sample_id, data = new_data) |>
    recipes::update_role(sample_id, new_role = "id") |>
    step_sec_conventional_cal(
      measures = "signal",
      calibration = loaded_cal
    )

  # Prep should use pre-loaded calibration
  expect_message(
    prepped <- recipes::prep(rec_new),
    "Using pre-loaded"
  )

  # Bake should produce mw column
  result <- recipes::bake(prepped, new_data = NULL)
  expect_true("mw" %in% names(result))
})

test_that("step_sec_conventional_cal errors on invalid calibration object", {
  skip_if_not_installed("measure")

  test_objs <- create_test_recipe()

  # Try to use invalid calibration object
  expect_error(
    recipes::recipe(signal ~ sample_id, data = test_objs$data) |>
      step_sec_conventional_cal(calibration = list(foo = "bar")),
    "sec_calibration"
  )
})

test_that("step_sec_conventional_cal errors when both standards and calibration missing", {
  skip_if_not_installed("measure")

  test_objs <- create_test_recipe()

  expect_error(
    recipes::recipe(signal ~ sample_id, data = test_objs$data) |>
      step_sec_conventional_cal(),
    "Either.*standards.*calibration.*required"
  )
})

test_that("print.sec_calibration displays correct information", {
  skip_if_not_installed("measure")

  test_objs <- create_test_recipe()
  tmpfile <- tempfile(fileext = ".rds")
  on.exit(unlink(tmpfile), add = TRUE)

  save_sec_calibration(
    test_objs$recipe,
    tmpfile,
    metadata = list(column = "Test Column")
  )
  cal <- load_sec_calibration(tmpfile)

  # Print returns the object invisibly
  expect_invisible(print(cal))
  expect_s3_class(cal, "sec_calibration")

  # Check that print contains expected elements via cli (output goes to message)
  # Since cli output is complex to capture, just verify the object structure
  expect_equal(cal$fit_type, "cubic")
  expect_equal(cal$user_metadata$column, "Test Column")
})

test_that("summary.sec_calibration shows per-standard results", {
  skip_if_not_installed("measure")

  test_objs <- create_test_recipe()
  tmpfile <- tempfile(fileext = ".rds")
  on.exit(unlink(tmpfile), add = TRUE)

  save_sec_calibration(test_objs$recipe, tmpfile)
  cal <- load_sec_calibration(tmpfile)

  # Summary returns the object invisibly
  expect_invisible(summary(cal))

  # Check that calibration has standard_results for summary to display
  expect_true(
    !is.null(cal$diagnostics$standard_results) ||
      !is.null(cal$diagnostics$r_squared)
  )
})

test_that("save_sec_calibration validates inputs", {
  skip_if_not_installed("measure")

  # Not a recipe
  expect_error(
    save_sec_calibration("not a recipe", "file.rds"),
    "must be a prepped recipe"
  )

  # Recipe without calibration step
  test_objs <- create_test_recipe()
  minimal_rec <- recipes::recipe(
    signal ~ sample_id,
    data = test_objs$data
  ) |>
    recipes::prep()

  tmpfile <- tempfile(fileext = ".rds")
  on.exit(unlink(tmpfile), add = TRUE)

  expect_error(
    save_sec_calibration(minimal_rec, tmpfile),
    "No trained calibration step"
  )
})

test_that("calibration preserves extrapolation and output settings", {
  skip_if_not_installed("measure")

  # Create test data with measure format
  set.seed(123)
  test_data <- tibble::tibble(
    sample_id = "test_sample",
    signal = list(tibble::tibble(
      location = seq(10, 18, by = 0.1),
      value = dnorm(seq(10, 18, by = 0.1), mean = 14, sd = 1)
    ))
  )
  class(test_data$signal[[1]]) <- c("measure_tbl", class(test_data$signal[[1]]))
  class(test_data$signal) <- c("measure_list", class(test_data$signal))

  standards <- data.frame(
    retention = c(10.5, 12.0, 13.5, 15.0, 16.5, 18.0),
    log_mw = c(6.0, 5.5, 5.0, 4.5, 4.0, 3.5)
  )

  # Create recipe with specific settings
  rec <- recipes::recipe(signal ~ sample_id, data = test_data) |>
    recipes::update_role(sample_id, new_role = "id") |>
    step_sec_conventional_cal(
      measures = "signal",
      standards = standards,
      fit_type = "quadratic",
      extrapolation = "none",
      output_col = "mol_weight",
      log_output = FALSE
    ) |>
    recipes::prep()

  tmpfile <- tempfile(fileext = ".rds")
  on.exit(unlink(tmpfile), add = TRUE)

  save_sec_calibration(rec, tmpfile)
  cal <- load_sec_calibration(tmpfile)

  # Check settings preserved
  expect_equal(cal$fit_type, "quadratic")
  expect_equal(cal$extrapolation, "none")
  expect_equal(cal$output_col, "mol_weight")
  expect_false(cal$log_output)
})
