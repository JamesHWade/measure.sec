test_that("linear extrapolation uses the boundary tangents of nonlinear calibration", {
  x <- 1:5
  y <- 6 - 0.3 * x + 0.02 * x^2
  data <- tibble::tibble(sample_id = "sample")
  data$ri <- measure::new_measure_list(list(
    measure::new_measure_tbl(0:6, rep(1, 7))
  ))
  rec <- recipes::recipe(data) |>
    step_sec_conventional_cal(
      measures = "ri",
      standards = data.frame(location = x, log_mw = y),
      fit_type = "quadratic",
      extrapolation = "linear"
    )
  out <- recipes::bake(recipes::prep(rec), new_data = NULL)$mw[[1]]$value
  expected <- c(
    y[1] + (-0.3 + 0.04 * x[1]) * (0 - x[1]),
    y[5] + (-0.3 + 0.04 * x[5]) * (6 - x[5])
  )
  expect_equal(unname(out[c(1, 7)]), expected, tolerance = 1e-7)
})

test_that("conventional calibration does not silently clip valid fitted values", {
  data <- tibble::tibble(sample_id = "sample")
  data$ri <- measure::new_measure_list(list(measure::new_measure_tbl(
    1:3,
    rep(1, 3)
  )))
  rec <- recipes::recipe(data) |>
    step_sec_conventional_cal(
      measures = "ri",
      standards = data.frame(location = 1:3, log_mw = c(0.3, 0.7, 1.1)),
      fit_type = "linear",
      extrapolation = "none"
    )
  out <- recipes::bake(recipes::prep(rec), new_data = NULL)
  expect_equal(unname(out$mw[[1]]$value), c(0.3, 0.7, 1.1), tolerance = 1e-8)
})

test_that("MALS rejects unsupported multi-angle and radius-of-gyration requests", {
  data <- tibble::tibble(sample_id = "sample")
  data$ri <- measure::new_measure_list(list(measure::new_measure_tbl(
    1:3,
    c(1, 2, 1)
  )))
  data$mals <- measure::new_measure_list(list(measure::new_measure_tbl(
    1:3,
    c(1, 2, 1)
  )))
  expect_error(
    step_sec_mals(recipes::recipe(data), dn_dc = 0.185, angles = c(45, 90)),
    "multi.angle|Multi.angle"
  )
  expect_error(
    step_sec_mals(recipes::recipe(data), dn_dc = 0.185, output_rg = "rg"),
    "Rg|gyration"
  )
})

test_that("MALS requires aligned detector signals and labels uncalibrated output", {
  data <- tibble::tibble(sample_id = "sample")
  data$ri <- measure::new_measure_list(list(measure::new_measure_tbl(
    c(1, 2, 3),
    c(1, 2, 1)
  )))
  data$mals <- measure::new_measure_list(list(measure::new_measure_tbl(
    c(1, 2.1, 3),
    c(1, 2, 1)
  )))
  rec <- recipes::recipe(data) |>
    step_sec_mals(
      mals_col = "mals",
      concentration_col = "ri",
      dn_dc = 0.185,
      calibration_constant = 1,
      concentration_unit = "mg/mL"
    )
  expect_error(recipes::prep(rec), "align|location")
  expect_warning(
    step_sec_mals(recipes::recipe(data), dn_dc = 0.185),
    "relative|uncalibrated"
  )
})

test_that("calibrated single-angle scattering respects concentration units", {
  concentration_mg <- c(1, 2, 1)
  dn_dc <- 0.185
  n0 <- 1.333
  lambda_cm <- 658e-7
  optical_k <- 4 * pi^2 * n0^2 * dn_dc^2 / (6.022e23 * lambda_cm^4)
  signal <- optical_k * (concentration_mg / 1000) * 1e5
  data <- tibble::tibble(sample_id = "sample")
  data$ri <- measure::new_measure_list(list(measure::new_measure_tbl(
    1:3,
    concentration_mg
  )))
  data$mals <- measure::new_measure_list(list(measure::new_measure_tbl(
    1:3,
    signal
  )))
  rec <- recipes::recipe(data) |>
    step_sec_mals(
      mals_col = "mals",
      concentration_col = "ri",
      dn_dc = dn_dc,
      calibration_constant = 1,
      concentration_unit = "mg/mL"
    )
  out <- recipes::bake(recipes::prep(rec), new_data = NULL)
  expect_equal(out$mw_mals[[1]]$value, rep(1e5, 3), tolerance = 1e-8)
  expect_equal(
    attr(out$mw_mals, "sec_scattering_basis"),
    "calibrated_single_angle"
  )
})

test_that("uncalibrated scattering cannot silently feed molecular weight averages", {
  data <- tibble::tibble(sample_id = "sample")
  data$ri <- measure::new_measure_list(list(measure::new_measure_tbl(
    1:3,
    c(1, 2, 1)
  )))
  data$mals <- measure::new_measure_list(list(measure::new_measure_tbl(
    1:3,
    c(1, 2, 1)
  )))
  expect_warning(
    rec <- step_sec_mals(
      recipes::recipe(data),
      mals_col = "mals",
      concentration_col = "ri",
      concentration_unit = "mg/mL",
      dn_dc = 0.185
    ),
    "uncalibrated"
  )
  expect_error(
    recipes::prep(step_sec_mw_averages(rec, measures = "ri")),
    "uncalibrated"
  )
  expect_error(
    recipes::prep(step_sec_mw_averages(
      rec,
      measures = "ri",
      mw_col = "mw_mals",
      mw_scale = "absolute"
    )),
    "uncalibrated"
  )
})

test_that("single-angle output is not silently clipped to arbitrary MW bounds", {
  dn_dc <- 0.185
  k <- 4 * pi^2 * 1.333^2 * dn_dc^2 / (6.022e23 * (658e-7)^4)
  data <- tibble::tibble(sample_id = "sample")
  data$ri <- measure::new_measure_list(list(measure::new_measure_tbl(
    1:3,
    rep(1, 3)
  )))
  data$mals <- measure::new_measure_list(list(measure::new_measure_tbl(
    1:3,
    rep(k * 0.001 * 50, 3)
  )))
  rec <- recipes::recipe(data) |>
    step_sec_mals(
      mals_col = "mals",
      concentration_col = "ri",
      concentration_unit = "mg/mL",
      dn_dc = dn_dc,
      calibration_constant = 1
    )
  out <- recipes::bake(recipes::prep(rec), new_data = NULL)
  expect_equal(out$mw_mals[[1]]$value, rep(50, 3), tolerance = 1e-8)
})

test_that("concentration normalization honors units and actual axis spacing", {
  data <- tibble::tibble(sample_id = "sample")
  data$ri <- measure::new_measure_list(list(measure::new_measure_tbl(
    c(1, 2, 4),
    c(0, 2, 0)
  )))
  make_result <- function(unit) {
    rec <- recipes::recipe(data) |>
      step_sec_concentration(
        measures = "ri",
        injection_mass = 3,
        flow_rate = 1,
        concentration_units = unit
      )
    recipes::bake(recipes::prep(rec), new_data = NULL)$ri
  }
  mg <- make_result("mg/mL")
  g <- make_result("g/mL")
  expect_equal(mg[[1]]$value, c(0, 2, 0))
  expect_equal(g[[1]]$value, mg[[1]]$value / 1000)
  expect_equal(attr(g, "sec_concentration_unit"), "g/mL")
})
