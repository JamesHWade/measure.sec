make_ri_test_data <- function(
  dn_dc,
  concentration_g_per_ml,
  instrument_constant
) {
  sample_id <- seq_along(dn_dc)
  data <- tibble::tibble(sample_id = sample_id, dn_dc = dn_dc)
  data$ri <- measure::new_measure_list(lapply(sample_id, function(i) {
    raw_signal <- concentration_g_per_ml[[i]] * dn_dc[[i]] / instrument_constant
    measure::new_measure_tbl(
      location = 1:3,
      value = rep(raw_signal, 3)
    )
  }))
  data
}

test_that("RI normalization accepts a negative scalar dn/dc", {
  concentration_g_per_ml <- 0.002
  dn_dc <- -0.25
  instrument_constant <- 2
  data <- make_ri_test_data(
    dn_dc = dn_dc,
    concentration_g_per_ml = concentration_g_per_ml,
    instrument_constant = instrument_constant
  )

  result <- recipes::recipe(data) |>
    step_sec_ri(
      measures = "ri",
      dn_dc = dn_dc,
      instrument_constant = instrument_constant
    ) |>
    recipes::prep() |>
    recipes::bake(new_data = NULL)

  expect_equal(result$ri[[1]]$value, rep(concentration_g_per_ml, 3))
  expect_identical(attr(result$ri, "sec_concentration_unit"), "g/mL")
})

test_that("RI normalization uses signed per-row dn/dc and selected units", {
  dn_dc <- c(0.2, -0.4, 0.8)
  concentration_g_per_ml <- c(0.001, 0.002, 0.003)
  instrument_constant <- 2
  data <- make_ri_test_data(
    dn_dc = dn_dc,
    concentration_g_per_ml = concentration_g_per_ml,
    instrument_constant = instrument_constant
  )
  row_order <- c(3, 1, 2)

  for (unit in c("g/mL", "mg/mL")) {
    prepped <- recipes::recipe(data) |>
      step_sec_ri(
        measures = "ri",
        dn_dc_column = "dn_dc",
        instrument_constant = instrument_constant,
        concentration_units = unit
      ) |>
      recipes::prep()
    serialized_recipe <- withr::local_tempfile()
    saveRDS(prepped, serialized_recipe)
    restored <- readRDS(serialized_recipe)

    expected <- if (identical(unit, "mg/mL")) {
      concentration_g_per_ml * 1000
    } else {
      concentration_g_per_ml
    }
    result <- recipes::bake(restored, new_data = data)
    reordered_result <- recipes::bake(
      restored,
      new_data = data[row_order, ]
    )

    expect_equal(
      lapply(result$ri, function(measure) measure$value),
      lapply(expected, function(value) rep(value, 3))
    )
    expect_equal(
      lapply(reordered_result$ri, function(measure) measure$value),
      lapply(expected[row_order], function(value) rep(value, 3))
    )
    expect_identical(attr(result$ri, "sec_concentration_unit"), unit)
    expect_identical(
      attr(reordered_result$ri, "sec_concentration_unit"),
      unit
    )
  }
})

test_that("RI without dn/dc preserves detector values and has no unit tag", {
  raw_values <- list(c(-2, 0, 1), c(3, -4, 5))
  data <- tibble::tibble(sample_id = seq_along(raw_values))
  data$ri <- measure::new_measure_list(lapply(raw_values, function(value) {
    measure::new_measure_tbl(location = seq_along(value), value = value)
  }))

  result <- recipes::recipe(data) |>
    step_sec_ri(measures = "ri", dn_dc = NULL) |>
    recipes::prep() |>
    recipes::bake(new_data = NULL)

  expect_identical(
    lapply(result$ri, function(measure) measure$value),
    raw_values
  )
  expect_null(attr(result$ri, "sec_concentration_unit"))
})

test_that("RI rejects zero and nonfinite dn/dc values", {
  recipe <- recipes::recipe(tibble::tibble(sample_id = "sample"))
  for (invalid in list(0, NA_real_, Inf, -Inf, c(0.2, -0.2))) {
    expect_error(step_sec_ri(recipe, dn_dc = invalid), "nonzero finite")
  }

  dn_dc <- c(0.2, -0.4, 0.8)
  concentration_g_per_ml <- c(0.001, 0.002, 0.003)
  instrument_constant <- 2
  training <- make_ri_test_data(
    dn_dc = dn_dc,
    concentration_g_per_ml = concentration_g_per_ml,
    instrument_constant = instrument_constant
  )
  prepped <- recipes::recipe(training) |>
    step_sec_ri(
      measures = "ri",
      dn_dc_column = "dn_dc",
      instrument_constant = instrument_constant
    ) |>
    recipes::prep()

  invalid_columns <- list(
    c(0, -0.4, 0.8),
    c(0.2, NA_real_, 0.8),
    c(0.2, -0.4, Inf),
    c(0.2, -0.4, -Inf),
    c("0.2", "-0.4", "0.8")
  )
  for (invalid in invalid_columns) {
    new_data <- training
    new_data$dn_dc <- invalid
    expect_error(
      recipes::bake(prepped, new_data = new_data),
      "Each dn/dc value must be nonzero and finite"
    )
  }
})
