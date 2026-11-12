sec_baseline_units_curve <- function(values, location = seq_along(values)) {
  measure::new_measure_tbl(location = location, value = values)
}

sec_baseline_units_data <- function(
  curves,
  unit = NULL,
  cache = NULL,
  ids = seq_along(curves),
  column = "ri"
) {
  data <- tibble::tibble(sample_id = ids)
  data[[column]] <- measure::new_measure_list(curves)
  if (!is.null(unit)) {
    attr(data[[column]], "sec_concentration_unit") <- unit
  }
  if (!is.null(cache)) {
    attr(data[[column]], "old_cached_metadata") <- cache
  }
  data
}

sec_baseline_units_fixture <- function(
  unit_scale = 1,
  peaks = c(2, 4, 2, 1)
) {
  location <- seq_len(20)
  value <- 10 + 0.1 * location
  value[9:12] <- value[9:12] + peaks

  corrected <- numeric(length(location))
  corrected[9:12] <- peaks

  list(
    curve = sec_baseline_units_curve(value * unit_scale, location),
    corrected = corrected * unit_scale
  )
}

sec_baseline_units_prep <- function(training, method = "linear") {
  recipes::recipe(training) |>
    step_sec_baseline(
      measures = "ri",
      left_frac = 0.2,
      right_frac = 0.2,
      method = method
    ) |>
    recipes::prep()
}

test_that("baseline methods correct values and retain mg/mL or g/mL tags", {
  methods <- c("linear", "median", "spline")
  if (requireNamespace("IDPmisc", quietly = TRUE)) {
    methods <- c(methods, "rf")
  }

  training_fixture <- sec_baseline_units_fixture()
  training <- sec_baseline_units_data(
    list(training_fixture$curve),
    unit = "mg/mL"
  )

  for (method in methods) {
    prepped <- sec_baseline_units_prep(training, method)

    for (unit in c("mg/mL", "g/mL")) {
      unit_scale <- if (identical(unit, "mg/mL")) 1 else 0.001
      input_fixture <- sec_baseline_units_fixture(unit_scale)
      input <- sec_baseline_units_data(
        list(input_fixture$curve),
        unit = unit
      )

      result <- recipes::bake(prepped, new_data = input)

      expect_identical(attr(result$ri, "sec_concentration_unit"), unit)
      expect_equal(
        result$ri[[1]]$value,
        input_fixture$corrected,
        tolerance = 1e-8
      )
      expect_false(identical(result$ri[[1]]$value, input$ri[[1]]$value))
    }
  }
})

test_that("bakes after RDS use each current unit and preserve row order only", {
  training_fixture <- sec_baseline_units_fixture()
  training <- sec_baseline_units_data(
    list(training_fixture$curve),
    unit = "mg/mL",
    cache = "from training"
  )
  prepped <- sec_baseline_units_prep(training)
  path <- withr::local_tempfile(fileext = ".rds")
  saveRDS(prepped, path)
  prepped <- readRDS(path)

  first_fixture <- sec_baseline_units_fixture(
    unit_scale = 0.001,
    peaks = c(1, 3, 2, 1)
  )
  second_fixture <- sec_baseline_units_fixture(
    unit_scale = 0.001,
    peaks = c(3, 2, 1, 2)
  )
  input <- sec_baseline_units_data(
    list(second_fixture$curve, first_fixture$curve),
    unit = "g/mL",
    cache = "stale cache",
    ids = c("second", "first")
  )
  input_before <- input
  result <- recipes::bake(prepped, new_data = input)

  expect_identical(result$sample_id, c("second", "first"))
  expect_identical(attr(result$ri, "sec_concentration_unit"), "g/mL")
  expect_equal(result$ri[[1]]$value, second_fixture$corrected, tolerance = 1e-8)
  expect_equal(result$ri[[2]]$value, first_fixture$corrected, tolerance = 1e-8)
  expect_null(attr(result$ri, "old_cached_metadata", exact = TRUE))
  expect_identical(input, input_before)

  mg_input <- sec_baseline_units_data(
    list(training_fixture$curve),
    unit = "mg/mL"
  )
  mg_result <- recipes::bake(prepped, new_data = mg_input)
  expect_identical(attr(mg_result$ri, "sec_concentration_unit"), "mg/mL")
  expect_equal(
    mg_result$ri[[1]]$value,
    training_fixture$corrected,
    tolerance = 1e-8
  )

  untagged <- sec_baseline_units_data(
    list(training_fixture$curve),
    cache = "must not survive"
  )
  untagged_result <- recipes::bake(prepped, new_data = untagged)
  expect_null(attr(untagged_result$ri, "sec_concentration_unit", exact = TRUE))
  expect_null(attr(untagged_result$ri, "old_cached_metadata", exact = TRUE))
})

test_that("baseline rejects invalid current tags and keeps fallback outputs tagged", {
  training_fixture <- sec_baseline_units_fixture()
  training <- sec_baseline_units_data(
    list(training_fixture$curve),
    unit = "mg/mL"
  )
  prepped <- sec_baseline_units_prep(training)

  invalid <- sec_baseline_units_data(
    list(training_fixture$curve),
    unit = "g/L"
  )
  expect_error(
    recipes::bake(prepped, new_data = invalid),
    "recorded concentration unit"
  )

  short_curve <- sec_baseline_units_curve(1:5, location = 1:5)
  short_input <- sec_baseline_units_data(list(short_curve), unit = "g/mL")
  expect_warning(
    short_result <- recipes::bake(prepped, new_data = short_input),
    "fewer than 10 points"
  )
  expect_identical(short_result$ri[[1]]$value, 1:5)
  expect_identical(attr(short_result$ri, "sec_concentration_unit"), "g/mL")

  missing_curve <- sec_baseline_units_curve(
    rep(NA_real_, 20),
    location = 1:20
  )
  missing_input <- sec_baseline_units_data(list(missing_curve), unit = "mg/mL")
  expect_warning(
    missing_result <- recipes::bake(prepped, new_data = missing_input),
    "all NA"
  )
  expect_identical(missing_result$ri[[1]]$value, rep(NA_real_, 20))
  expect_identical(attr(missing_result$ri, "sec_concentration_unit"), "mg/mL")
})

test_that("zero-row and untagged baseline columns retain their unit state", {
  training_fixture <- sec_baseline_units_fixture()
  training <- sec_baseline_units_data(
    list(training_fixture$curve),
    unit = "mg/mL"
  )
  prepped <- sec_baseline_units_prep(training)

  empty_tagged <- sec_baseline_units_data(list(), unit = "g/mL")
  tagged_result <- recipes::bake(prepped, new_data = empty_tagged)
  expect_identical(nrow(tagged_result), 0L)
  expect_length(tagged_result$ri, 0)
  expect_identical(attr(tagged_result$ri, "sec_concentration_unit"), "g/mL")

  empty_untagged <- sec_baseline_units_data(list())
  untagged_result <- recipes::bake(prepped, new_data = empty_untagged)
  expect_identical(nrow(untagged_result), 0L)
  expect_null(attr(untagged_result$ri, "sec_concentration_unit", exact = TRUE))
})

test_that("RI baseline output keeps its unit for downstream LALS calculation", {
  location <- seq_len(20)
  dn_dc <- 0.185
  baseline_concentration <- 0.001 + 0.0001 * location
  peak_concentration <- numeric(length(location))
  peak_concentration[9:12] <- c(0.001, 0.003, 0.002, 0.001)
  ri_concentration <- baseline_concentration + peak_concentration
  wavelength_nm <- 670
  wavelength_cm <- wavelength_nm * 1e-7
  optical_constant <- 4 *
    pi^2 *
    1.333^2 *
    dn_dc^2 /
    (6.022e23 * wavelength_cm^4)
  expected_mw <- 50
  scattering <- optical_constant * peak_concentration * expected_mw

  data <- tibble::tibble(sample_id = "sample")
  data$ri_signal <- measure::new_measure_list(list(
    sec_baseline_units_curve(ri_concentration * dn_dc, location)
  ))
  data$lals <- measure::new_measure_list(list(
    sec_baseline_units_curve(scattering, location)
  ))

  recipe <- recipes::recipe(data) |>
    step_sec_ri(
      measures = "ri_signal",
      output_col = "ri",
      dn_dc = dn_dc,
      concentration_units = "g/mL"
    ) |>
    step_sec_baseline(
      measures = "ri",
      left_frac = 0.2,
      right_frac = 0.2
    ) |>
    step_sec_lals(
      measures = "lals",
      concentration_col = "ri",
      dn_dc = dn_dc,
      calibration_constant = 1,
      min_signal = 0
    )

  result <- recipes::bake(recipes::prep(recipe), new_data = NULL)

  expect_identical(attr(result$ri, "sec_concentration_unit"), "g/mL")
  expect_equal(
    result$ri[[1]]$value[10],
    peak_concentration[10],
    tolerance = 1e-10
  )
  expect_equal(result$mw_lals[[1]]$value[10], expected_mw, tolerance = 1e-8)
})
