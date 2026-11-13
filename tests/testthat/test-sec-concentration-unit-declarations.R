sec_unit_declaration_data <- function(concentration, recorded_unit = NULL) {
  location <- 0:2
  data <- tibble::tibble(sample_id = "sample")
  data$ri <- measure::new_measure_list(list(
    measure::new_measure_tbl(location, concentration)
  ))
  if (!is.null(recorded_unit)) {
    attr(data$ri, "sec_concentration_unit") <- recorded_unit
  }
  data$lals <- measure::new_measure_list(list(
    measure::new_measure_tbl(location, c(2.5e-6, 5e-6, 7.5e-6))
  ))
  data
}

sec_unit_declaration_lals <- function(data, concentration_unit = NULL) {
  recipe <- recipes::recipe(data) |>
    step_sec_lals(
      measures = "lals",
      concentration_col = "ri",
      concentration_unit = concentration_unit,
      optical_constant = 1e-7,
      calibration_constant = 2,
      min_signal = 0
    )

  recipes::bake(recipes::prep(recipe), new_data = data)$mw_lals[[1]]$value
}

test_that("concentration unit declarations are plain, canonical scalars", {
  expect_no_error(.sec_validate_concentration_unit(NULL))
  expect_no_error(
    .sec_validate_concentration_unit(structure("mg/mL", names = "unit"))
  )
  expect_identical(
    .sec_validate_concentration_unit(structure("mg/mL", names = "unit")),
    "mg/mL"
  )
  expect_error(
    .sec_validate_concentration_unit(matrix("mg/mL")),
    "concentration_unit"
  )
  expect_error(
    .sec_validate_concentration_unit(structure("mg/mL", class = "unit")),
    "concentration_unit"
  )
  expect_error(
    .sec_validate_concentration_unit(NA_character_),
    "concentration_unit"
  )
  expect_error(
    .sec_validate_concentration_unit(NULL, allow_null = FALSE),
    "concentration_unit"
  )

  expect_equal(
    .sec_concentration_g_per_ml(1000, structure("mg/mL", names = "unit")),
    1
  )
  expect_equal(.sec_concentration_g_per_ml(1, "g/mL"), 1)
  expect_error(.sec_concentration_g_per_ml(1, NULL), "concentration_unit")
  expect_error(
    .sec_concentration_g_per_ml(1, matrix("mg/mL")),
    "concentration_unit"
  )
})

test_that("LALS converts named and recorded concentration units to g/mL", {
  concentration_mg_per_ml <- c(1, 2, 3)
  expected <- rep(50000, 3)
  plain_data <- sec_unit_declaration_data(concentration_mg_per_ml, "mg/mL")
  named_explicit <- sec_unit_declaration_lals(
    plain_data,
    structure("mg/mL", names = "unit")
  )
  plain_explicit <- sec_unit_declaration_lals(plain_data, "mg/mL")

  named_recorded_data <- sec_unit_declaration_data(
    concentration_mg_per_ml,
    structure("mg/mL", names = "recorded")
  )
  named_recorded <- sec_unit_declaration_lals(
    named_recorded_data,
    "mg/mL"
  )

  concentration_g_per_ml <- concentration_mg_per_ml / 1000
  gram_data <- sec_unit_declaration_data(concentration_g_per_ml, "g/mL")
  gram_tagged <- sec_unit_declaration_lals(gram_data, "g/mL")

  expect_equal(plain_explicit, expected, tolerance = 1e-8)
  expect_equal(named_explicit, plain_explicit, tolerance = 1e-8)
  expect_equal(named_recorded, plain_explicit, tolerance = 1e-8)
  expect_equal(gram_tagged, plain_explicit, tolerance = 1e-8)

  named_tag_explicit_canonical <- sec_unit_declaration_lals(
    named_recorded_data,
    structure("mg/mL", names = "explicit")
  )
  expect_equal(named_tag_explicit_canonical, plain_explicit, tolerance = 1e-8)
})

test_that("LALS rejects matrix and classed concentration unit declarations", {
  data <- sec_unit_declaration_data(c(1, 2, 3))
  recipe <- recipes::recipe(data)

  expect_error(
    step_sec_lals(
      recipe,
      measures = "lals",
      concentration_col = "ri",
      concentration_unit = matrix("mg/mL"),
      optical_constant = 1e-7,
      calibration_constant = 2
    ),
    "concentration_unit"
  )
  expect_error(
    step_sec_lals(
      recipe,
      measures = "lals",
      concentration_col = "ri",
      concentration_unit = structure("mg/mL", class = "unit"),
      optical_constant = 1e-7,
      calibration_constant = 2
    ),
    "concentration_unit"
  )

  matrix_tag <- sec_unit_declaration_data(c(1, 2, 3), matrix("mg/mL"))
  expect_error(
    recipes::prep(
      recipes::recipe(matrix_tag) |>
        step_sec_lals(
          measures = "lals",
          concentration_col = "ri",
          optical_constant = 1e-7,
          calibration_constant = 2
        )
    ),
    "recorded concentration unit"
  )

  classed_tag <- sec_unit_declaration_data(
    c(1, 2, 3),
    structure("mg/mL", class = "unit")
  )
  expect_error(
    recipes::prep(
      recipes::recipe(classed_tag) |>
        step_sec_lals(
          measures = "lals",
          concentration_col = "ri",
          optical_constant = 1e-7,
          calibration_constant = 2
        )
    ),
    "recorded concentration unit"
  )
})
