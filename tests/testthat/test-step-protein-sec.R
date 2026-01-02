# ==============================================================================
# Tests for step_sec_oligomer and step_sec_protein
# ==============================================================================

# -- Helper functions ----------------------------------------------------------

#' Create simulated protein SEC data with multiple species
create_protein_sec_data <- function(
  monomer_pct = 95,
  dimer_pct = 3,
  hmw_pct = 1,
  lmw_pct = 1,
  n_points = 300
) {
  time <- seq(5, 25, length.out = n_points)

  # Monomer peak (main peak at ~15 min)
  monomer <- monomer_pct * dnorm(time, mean = 15, sd = 0.5)

  # Dimer peak (elutes earlier at ~13 min)
  dimer <- dimer_pct * dnorm(time, mean = 13, sd = 0.4)

  # HMW (even earlier at ~11 min)
  hmw <- hmw_pct * dnorm(time, mean = 11, sd = 0.3)

  # LMW/fragments (later at ~18 min)
  lmw <- lmw_pct * dnorm(time, mean = 18, sd = 0.4)

  # Combined signal
  signal <- monomer + dimer + hmw + lmw

  # Add small noise
  set.seed(42)
  signal <- signal + rnorm(n_points, 0, 0.001 * max(signal))
  signal <- pmax(signal, 0)

  list(
    time = time,
    signal = signal,
    expected_monomer = monomer_pct,
    expected_dimer = dimer_pct,
    expected_hmw = hmw_pct,
    expected_lmw = lmw_pct
  )
}

#' Create test data tibble from simulated data
create_test_protein_tibble <- function(sec_data) {
  test_data <- tibble::tibble(sample_id = "test")
  test_data$uv <- measure::new_measure_list(
    list(measure::new_measure_tbl(
      location = sec_data$time,
      value = sec_data$signal
    ))
  )
  test_data
}

# ==============================================================================
# step_sec_oligomer tests
# ==============================================================================

test_that("step_sec_oligomer requires monomer_mw", {
  skip_if_not_installed("measure")

  sec_data <- create_protein_sec_data()
  test_data <- create_test_protein_tibble(sec_data)

  expect_error(
    recipes::recipe(~., data = test_data) |>
      step_sec_oligomer(),
    "monomer_mw.*required"
  )
})

test_that("step_sec_oligomer validates monomer_mw is positive", {
  skip_if_not_installed("measure")

  sec_data <- create_protein_sec_data()
  test_data <- create_test_protein_tibble(sec_data)

  expect_error(
    recipes::recipe(~., data = test_data) |>
      step_sec_oligomer(monomer_mw = -150000),
    "positive"
  )

  expect_error(
    recipes::recipe(~., data = test_data) |>
      step_sec_oligomer(monomer_mw = 0),
    "positive"
  )
})

test_that("step_sec_oligomer warns for unusual monomer_mw", {
  skip_if_not_installed("measure")

  sec_data <- create_protein_sec_data()
  test_data <- create_test_protein_tibble(sec_data)

  # Too small
  expect_warning(
    recipes::recipe(~., data = test_data) |>
      step_sec_oligomer(monomer_mw = 1000),
    "outside.*protein range"
  )

  # Too large
  expect_warning(
    recipes::recipe(~., data = test_data) |>
      step_sec_oligomer(monomer_mw = 5000000),
    "outside.*protein range"
  )
})

test_that("step_sec_oligomer validates mw_tolerance", {
  skip_if_not_installed("measure")

  sec_data <- create_protein_sec_data()
  test_data <- create_test_protein_tibble(sec_data)

  expect_error(
    recipes::recipe(~., data = test_data) |>
      step_sec_oligomer(monomer_mw = 150000, mw_tolerance = 0),
    "between 0 and 1"
  )

  expect_error(
    recipes::recipe(~., data = test_data) |>
      step_sec_oligomer(monomer_mw = 150000, mw_tolerance = 1.5),
    "between 0 and 1"
  )
})

test_that("step_sec_oligomer detects species", {
  skip_if_not_installed("measure")

  sec_data <- create_protein_sec_data(
    monomer_pct = 90,
    dimer_pct = 5,
    hmw_pct = 3,
    lmw_pct = 2
  )
  test_data <- create_test_protein_tibble(sec_data)

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_oligomer(monomer_mw = 150000)

  prepped <- recipes::prep(rec)
  result <- recipes::bake(prepped, new_data = NULL)

  # Should have output columns
  expect_true("oligo_monomer_pct" %in% names(result))
  expect_true("oligo_hmw_pct" %in% names(result))
  expect_true("oligo_lmw_pct" %in% names(result))
  expect_true("oligo_species_count" %in% names(result))

  # Monomer should be the largest
  expect_gt(result$oligo_monomer_pct[1], 50)

  # Species count should be reasonable
  expect_gt(result$oligo_species_count[1], 0)
})

test_that("step_sec_oligomer handles single peak", {
  skip_if_not_installed("measure")

  # Pure monomer
  sec_data <- create_protein_sec_data(
    monomer_pct = 100,
    dimer_pct = 0,
    hmw_pct = 0,
    lmw_pct = 0
  )
  test_data <- create_test_protein_tibble(sec_data)

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_oligomer(monomer_mw = 150000)

  prepped <- recipes::prep(rec)
  result <- recipes::bake(prepped, new_data = NULL)

  # Should have nearly all monomer
  expect_gt(result$oligo_monomer_pct[1], 90)
})

test_that("step_sec_oligomer print method works", {
  skip_if_not_installed("measure")

  sec_data <- create_protein_sec_data()
  test_data <- create_test_protein_tibble(sec_data)

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_oligomer(monomer_mw = 150000)

  expect_output(print(rec), "oligomer")

  prepped <- recipes::prep(rec)
  expect_output(print(prepped), "150,000")
})

test_that("step_sec_oligomer tidy method works", {
  skip_if_not_installed("measure")

  sec_data <- create_protein_sec_data()
  test_data <- create_test_protein_tibble(sec_data)

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_oligomer(monomer_mw = 150000, mw_tolerance = 0.2)

  prepped <- recipes::prep(rec)
  tidy_result <- recipes::tidy(prepped, number = 1)

  expect_s3_class(tidy_result, "tbl_df")
  expect_equal(tidy_result$monomer_mw, 150000)
  expect_equal(tidy_result$mw_tolerance, 0.2)
})

# ==============================================================================
# step_sec_protein tests
# ==============================================================================

test_that("step_sec_protein works with basic parameters", {
  skip_if_not_installed("measure")

  sec_data <- create_protein_sec_data()
  test_data <- create_test_protein_tibble(sec_data)

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_protein(monomer_mw = 150000)

  prepped <- recipes::prep(rec)
  result <- recipes::bake(prepped, new_data = NULL)

  # Should have aggregate columns
  expect_true("protein_hmws_pct" %in% names(result))
  expect_true("protein_monomer_pct" %in% names(result))
  expect_true("protein_lmws_pct" %in% names(result))
  expect_true("protein_main_start" %in% names(result))
  expect_true("protein_main_end" %in% names(result))

  # Oligomer columns should also be present (include_oligomer defaults to TRUE)
  expect_true("protein_dimer_pct" %in% names(result))
  expect_true("protein_species_count" %in% names(result))
})

test_that("step_sec_protein works without oligomer analysis", {
  skip_if_not_installed("measure")

  sec_data <- create_protein_sec_data()
  test_data <- create_test_protein_tibble(sec_data)

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_protein(include_oligomer = FALSE)

  prepped <- recipes::prep(rec)
  result <- recipes::bake(prepped, new_data = NULL)

  # Should have aggregate columns
  expect_true("protein_hmws_pct" %in% names(result))
  expect_true("protein_monomer_pct" %in% names(result))

  # Should NOT have oligomer columns
  expect_false("protein_dimer_pct" %in% names(result))
  expect_false("protein_species_count" %in% names(result))
})

test_that("step_sec_protein validates oligomer requires monomer_mw", {
  skip_if_not_installed("measure")

  sec_data <- create_protein_sec_data()
  test_data <- create_test_protein_tibble(sec_data)

  expect_error(
    recipes::recipe(~., data = test_data) |>
      step_sec_protein(include_oligomer = TRUE),
    "monomer_mw"
  )
})

test_that("step_sec_protein calculates reasonable percentages", {
  skip_if_not_installed("measure")

  sec_data <- create_protein_sec_data(
    monomer_pct = 95,
    dimer_pct = 3,
    hmw_pct = 1,
    lmw_pct = 1
  )
  test_data <- create_test_protein_tibble(sec_data)

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_protein(monomer_mw = 150000)

  prepped <- recipes::prep(rec)
  result <- recipes::bake(prepped, new_data = NULL)

  # Percentages should sum to ~100
  total <- result$protein_hmws_pct[1] +
    result$protein_monomer_pct[1] +
    result$protein_lmws_pct[1]
  expect_equal(total, 100, tolerance = 5)

  # Monomer should dominate
  expect_gt(result$protein_monomer_pct[1], 80)
})

test_that("step_sec_protein supports different baseline methods", {
  skip_if_not_installed("measure")

  sec_data <- create_protein_sec_data()
  test_data <- create_test_protein_tibble(sec_data)

  # Linear
  rec_linear <- recipes::recipe(~., data = test_data) |>
    step_sec_protein(
      include_oligomer = FALSE,
      baseline_method = "linear"
    )
  result_linear <- recipes::bake(recipes::prep(rec_linear), new_data = NULL)

  # Median
  rec_median <- recipes::recipe(~., data = test_data) |>
    step_sec_protein(
      include_oligomer = FALSE,
      baseline_method = "median"
    )
  result_median <- recipes::bake(recipes::prep(rec_median), new_data = NULL)

  # Both should produce valid results
  expect_true(result_linear$protein_monomer_pct[1] > 0)
  expect_true(result_median$protein_monomer_pct[1] > 0)
})

test_that("step_sec_protein print method works", {
  skip_if_not_installed("measure")

  sec_data <- create_protein_sec_data()
  test_data <- create_test_protein_tibble(sec_data)

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_protein(monomer_mw = 150000, type = "native")

  expect_output(print(rec), "native")

  prepped <- recipes::prep(rec)
  expect_output(print(prepped), "150,000")
})

test_that("step_sec_protein tidy method works", {
  skip_if_not_installed("measure")

  sec_data <- create_protein_sec_data()
  test_data <- create_test_protein_tibble(sec_data)

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_protein(
      monomer_mw = 150000,
      type = "denaturing",
      baseline_method = "median"
    )

  prepped <- recipes::prep(rec)
  tidy_result <- recipes::tidy(prepped, number = 1)

  expect_s3_class(tidy_result, "tbl_df")
  expect_equal(tidy_result$type, "denaturing")
  expect_equal(tidy_result$monomer_mw, 150000)
  expect_equal(tidy_result$baseline_method, "median")
})

test_that("step_sec_protein handles type argument", {
  skip_if_not_installed("measure")

  sec_data <- create_protein_sec_data()
  test_data <- create_test_protein_tibble(sec_data)

  # Native type
  rec_native <- recipes::recipe(~., data = test_data) |>
    step_sec_protein(type = "native", include_oligomer = FALSE)

  prepped_native <- recipes::prep(rec_native)
  expect_equal(prepped_native$steps[[1]]$type, "native")

  # Denaturing type
  rec_denat <- recipes::recipe(~., data = test_data) |>
    step_sec_protein(type = "denaturing", include_oligomer = FALSE)

  prepped_denat <- recipes::prep(rec_denat)
  expect_equal(prepped_denat$steps[[1]]$type, "denaturing")
})

# ==============================================================================
# Additional tests for mw_column, manual peaks, spline baseline, and warnings
# ==============================================================================

#' Create test data with MW column for MW-based species assignment
create_test_protein_tibble_with_mw <- function(sec_data) {
  test_data <- tibble::tibble(sample_id = "test")
  test_data$uv <- measure::new_measure_list(
    list(measure::new_measure_tbl(
      location = sec_data$time,
      value = sec_data$signal
    ))
  )
  # Add simulated MW column (decreasing with time, as in SEC)
  # Monomer at 150 kDa, dimer at 300 kDa, HMW at 450+ kDa
  mw_values <- 1000000 * exp(-0.15 * sec_data$time)
  test_data$mw_mals <- measure::new_measure_list(
    list(measure::new_measure_tbl(
      location = sec_data$time,
      value = mw_values
    ))
  )
  test_data
}

test_that("step_sec_oligomer uses mw_column for species assignment", {
  skip_if_not_installed("measure")

  sec_data <- create_protein_sec_data(
    monomer_pct = 90,
    dimer_pct = 5,
    hmw_pct = 3,
    lmw_pct = 2
  )
  test_data <- create_test_protein_tibble_with_mw(sec_data)

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_oligomer(monomer_mw = 150000, mw_column = "mw_mals")

  prepped <- recipes::prep(rec)

  # Should not warn about RT fallback when mw_column is provided
  expect_no_warning(
    result <- recipes::bake(prepped, new_data = NULL)
  )

  # Should have MW columns in output
  expect_true("oligo_monomer_mw" %in% names(result))
})

test_that("step_sec_oligomer warns when mw_column not found", {
  skip_if_not_installed("measure")

  sec_data <- create_protein_sec_data()
  test_data <- create_test_protein_tibble(sec_data)

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_oligomer(monomer_mw = 150000, mw_column = "nonexistent_col")

  expect_warning(
    recipes::prep(rec),
    "not found"
  )
})

test_that("step_sec_oligomer warns about RT fallback", {
  skip_if_not_installed("measure")

  sec_data <- create_protein_sec_data(
    monomer_pct = 90,
    dimer_pct = 5,
    hmw_pct = 3,
    lmw_pct = 2
  )
  test_data <- create_test_protein_tibble(sec_data)

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_oligomer(monomer_mw = 150000)

  # Warning is emitted during prep (when training data is baked)
  expect_warning(
    recipes::prep(rec),
    "retention time patterns"
  )
})

test_that("step_sec_oligomer works with manual peak detection", {
  skip_if_not_installed("measure")

  sec_data <- create_protein_sec_data()
  test_data <- create_test_protein_tibble(sec_data)

  # Define manual peaks matching the simulated data
  manual_peaks <- data.frame(
    start = c(10, 12, 14, 17),
    end = c(12, 14, 16, 19)
  )

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_oligomer(
      monomer_mw = 150000,
      peak_detection = "manual",
      peaks = manual_peaks
    )

  prepped <- recipes::prep(rec)
  result <- suppressWarnings(recipes::bake(prepped, new_data = NULL))

  # Should detect the manually defined peaks
  expect_gt(result$oligo_species_count[1], 0)
  expect_true("oligo_monomer_pct" %in% names(result))
})

test_that("step_sec_oligomer requires peaks for manual mode", {
  skip_if_not_installed("measure")

  sec_data <- create_protein_sec_data()
  test_data <- create_test_protein_tibble(sec_data)

  expect_error(
    recipes::recipe(~., data = test_data) |>
      step_sec_oligomer(monomer_mw = 150000, peak_detection = "manual"),
    "peaks"
  )
})

test_that("step_sec_protein uses spline baseline", {
  skip_if_not_installed("measure")

  sec_data <- create_protein_sec_data()
  test_data <- create_test_protein_tibble(sec_data)

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_protein(
      include_oligomer = FALSE,
      baseline_method = "spline"
    )

  # Should work without error
  prepped <- recipes::prep(rec)
  result <- recipes::bake(prepped, new_data = NULL)

  expect_true("protein_monomer_pct" %in% names(result))
  expect_gt(result$protein_monomer_pct[1], 0)
})

test_that("step_sec_protein warns for spline with few points", {
  skip_if_not_installed("measure")

  # Create data with only 8 points (< 10)
  sec_data <- create_protein_sec_data(n_points = 8)
  test_data <- create_test_protein_tibble(sec_data)

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_protein(
      include_oligomer = FALSE,
      baseline_method = "spline"
    )

  # Warning is emitted during prep (when training data is baked)
  expect_warning(
    recipes::prep(rec),
    "linear interpolation"
  )
})

test_that("step_sec_protein type='denaturing' works correctly", {
  skip_if_not_installed("measure")

  sec_data <- create_protein_sec_data()
  test_data <- create_test_protein_tibble(sec_data)

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_protein(type = "denaturing", include_oligomer = FALSE)

  # Should run without error and produce valid results
  prepped <- recipes::prep(rec)
  result <- recipes::bake(prepped, new_data = NULL)

  expect_true("protein_monomer_pct" %in% names(result))
  expect_gt(result$protein_monomer_pct[1], 0)

  # Verify type is stored correctly
  expect_equal(prepped$steps[[1]]$type, "denaturing")
})
