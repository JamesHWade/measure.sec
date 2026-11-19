test_that("a smoothed RI concentration drives LALS on a new injection", {
  concentration <- c(0.4, 1, 2, 1, 0.4) / 1000 # g/mL
  dn_dc <- 0.185
  optical_constant <- 2e-7
  molar_mass <- 66430
  trace <- function(value) {
    measure::new_measure_list(list(measure::new_measure_tbl(8:12, value)))
  }
  injection <- function(amount) {
    tibble::tibble(
      ri = trace(amount * concentration * dn_dc),
      lals = trace(amount * concentration * optical_constant * molar_mass)
    )
  }

  training <- injection(1)
  held_out <- injection(0.7)
  scattering <- function(recipe, unit = NULL) {
    step_sec_lals(
      recipe,
      measures = "lals",
      concentration_col = "concentration",
      concentration_unit = unit,
      optical_constant = optical_constant,
      calibration_constant = 1
    )
  }

  for (unit in c("mg/mL", "g/mL")) {
    smoothed <- recipes::recipe(~., data = training) |>
      step_sec_ri(
        measures = "ri",
        dn_dc = dn_dc,
        output_col = "concentration",
        concentration_units = unit
      ) |>
      measure::step_measure_smooth_ma(
        measures = c("concentration", "lals"),
        window = 3
      )

    result <- recipes::bake(recipes::prep(scattering(smoothed)), held_out)
    expected <- 0.7 * c(0.8, 3.4 / 3, 4 / 3, 3.4 / 3, 0.8)
    if (unit == "g/mL") {
      expected <- expected / 1000
    }

    expect_equal(result$concentration[[1]]$value, expected)
    expect_equal(result$mw_lals[[1]]$value, rep(molar_mass, 5))
    expect_identical(
      measure::measure_value_metadata(result$concentration),
      list(quantity = "mass_concentration", unit = unit)
    )

    normalized <- measure::step_measure_normalize_max(
      smoothed,
      measures = "concentration"
    )
    expect_error(recipes::prep(scattering(normalized)), "normalized_signal")
    expect_error(
      recipes::prep(scattering(normalized, unit)),
      "normalized_signal"
    )

    differentiated <- measure::step_measure_derivative(
      smoothed,
      measures = "concentration"
    )
    expect_error(
      recipes::prep(scattering(differentiated, unit)),
      "signal_derivative"
    )
  }
})
