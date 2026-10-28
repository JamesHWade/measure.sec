test_that("RI normalization reports concentration in selected units", {
  concentration_g_per_ml <- c(0.001, 0.002, 0.003)
  dn_dc <- 0.185
  data <- tibble::tibble(sample_id = "sample")
  data$ri <- measure::new_measure_list(list(
    measure::new_measure_tbl(
      location = 1:3,
      value = concentration_g_per_ml * dn_dc
    )
  ))

  normalize_ri <- function(units) {
    rec <- recipes::recipe(data) |>
      step_sec_ri(dn_dc = dn_dc, concentration_units = units)
    recipes::bake(recipes::prep(rec), new_data = NULL)
  }

  result_g <- normalize_ri("g/mL")
  result_mg <- normalize_ri("mg/mL")

  expect_equal(result_g$ri[[1]]$value, concentration_g_per_ml)
  expect_equal(result_mg$ri[[1]]$value, concentration_g_per_ml * 1000)
  expect_identical(attr(result_g$ri, "sec_concentration_unit"), "g/mL")
  expect_identical(attr(result_mg$ri, "sec_concentration_unit"), "mg/mL")
})

test_that("UV normalization reports concentration in selected units", {
  concentration_mg_per_ml <- c(0.5, 1, 2)
  extinction_coef <- 2.5
  path_length <- 2
  absorbance <- concentration_mg_per_ml * extinction_coef * path_length
  data <- tibble::tibble(sample_id = "sample")
  data$uv <- measure::new_measure_list(list(
    measure::new_measure_tbl(location = 1:3, value = absorbance)
  ))

  normalize_uv <- function(units) {
    rec <- recipes::recipe(data) |>
      step_sec_uv(
        extinction_coef = extinction_coef,
        path_length = path_length,
        concentration_units = units
      )
    recipes::bake(recipes::prep(rec), new_data = NULL)
  }

  result_mg <- normalize_uv("mg/mL")
  result_g <- normalize_uv("g/mL")

  expect_equal(result_mg$uv[[1]]$value, concentration_mg_per_ml)
  expect_equal(result_g$uv[[1]]$value, concentration_mg_per_ml / 1000)
  expect_identical(attr(result_mg$uv, "sec_concentration_unit"), "mg/mL")
  expect_identical(attr(result_g$uv, "sec_concentration_unit"), "g/mL")
})

test_that("mass normalization clips negative points and reports their share", {
  data <- tibble::tibble(sample_id = "sample")
  data$ri <- measure::new_measure_list(list(
    measure::new_measure_tbl(
      location = c(0, 1, 2),
      value = c(-1, 2, 1)
    )
  ))
  rec <- recipes::recipe(data) |>
    step_sec_concentration(
      measures = "ri",
      injection_mass = 5,
      flow_rate = 1,
      concentration_units = "mg/mL"
    )

  expect_warning(
    result <- recipes::bake(recipes::prep(rec), new_data = NULL),
    "1 of 3 signal points"
  )

  expect_equal(result$ri[[1]]$value, c(0, 4, 2))
  expect_identical(attr(result$ri, "sec_concentration_unit"), "mg/mL")
})

test_that("LALS and RALS respect concentration units in the single-angle formula", {
  dn_dc <- 0.185
  solvent_ri <- 1.333
  wavelength_nm <- 670
  lambda_cm <- wavelength_nm * 1e-7
  optical_constant <- 4 *
    pi^2 *
    solvent_ri^2 *
    dn_dc^2 /
    (6.022e23 * lambda_cm^4)
  concentration_mg_per_ml <- c(1, 2, 3)
  expected_mw <- 50
  calibration_constant <- 2
  scattering_signal <- optical_constant *
    (concentration_mg_per_ml / 1000) *
    expected_mw /
    calibration_constant
  location <- c(4, 2, 1)

  calculate <- function(detector, unit, step_fn) {
    concentration <- if (unit == "mg/mL") {
      concentration_mg_per_ml
    } else {
      concentration_mg_per_ml / 1000
    }
    data <- tibble::tibble(sample_id = "sample")
    data$conc <- measure::new_measure_list(list(
      measure::new_measure_tbl(location = location, value = concentration)
    ))
    attr(data$conc, "sec_concentration_unit") <- unit
    data[[detector]] <- measure::new_measure_list(list(
      measure::new_measure_tbl(location = location, value = scattering_signal)
    ))

    rec <- recipes::recipe(data) |>
      step_fn(
        measures = detector,
        concentration_col = "conc",
        dn_dc = dn_dc,
        solvent_ri = solvent_ri,
        laser_wavelength = wavelength_nm,
        calibration_constant = calibration_constant
      )
    recipes::bake(recipes::prep(rec), new_data = NULL)
  }

  for (detector_step in list(
    list(name = "lals", fn = step_sec_lals),
    list(name = "rals", fn = step_sec_rals)
  )) {
    for (unit in c("mg/mL", "g/mL")) {
      result <- calculate(detector_step$name, unit, detector_step$fn)
      output <- result[[paste0("mw_", detector_step$name)]]

      expect_equal(output[[1]]$value, rep(expected_mw, 3), tolerance = 1e-8)
      expect_identical(
        attr(output, "sec_scattering_basis"),
        "calibrated_single_angle"
      )
      expect_identical(
        attr(output, "sec_mw_curve"),
        list(source = "conc", scale = "absolute")
      )
    }
  }
})

test_that("scattering steps require units and aligned concentration channels", {
  location <- c(4, 2, 1)
  data <- tibble::tibble(sample_id = "sample")
  data$conc <- measure::new_measure_list(list(
    measure::new_measure_tbl(location = location, value = c(1, 2, 3))
  ))
  data$lals <- measure::new_measure_list(list(
    measure::new_measure_tbl(location = location, value = c(1, 2, 3))
  ))
  data$rals <- measure::new_measure_list(list(
    measure::new_measure_tbl(location = location, value = c(1, 2, 3))
  ))

  expect_error(
    recipes::prep(
      recipes::recipe(data) |>
        step_sec_lals(
          measures = "lals",
          concentration_col = "conc",
          dn_dc = 0.185,
          calibration_constant = 1
        )
    ),
    "Specify.*concentration_unit"
  )

  attr(data$conc, "sec_concentration_unit") <- "mg/mL"
  expect_error(
    recipes::prep(
      recipes::recipe(data) |>
        step_sec_rals(
          measures = "rals",
          concentration_col = "conc",
          concentration_unit = "g/mL",
          dn_dc = 0.185,
          calibration_constant = 1
        )
    ),
    "conflicts"
  )

  data$lals <- measure::new_measure_list(list(
    measure::new_measure_tbl(location = c(4, 2.1, 1), value = c(1, 2, 3))
  ))
  expect_error(
    {
      rec <- recipes::recipe(data) |>
        step_sec_lals(
          measures = "lals",
          concentration_col = "conc",
          dn_dc = 0.185,
          calibration_constant = 1
        )
      recipes::bake(recipes::prep(rec), new_data = NULL)
    },
    "locations must align"
  )
})

test_that("detector calibration arguments and units are validated", {
  data <- tibble::tibble(sample_id = "sample")

  expect_error(
    step_sec_ri(recipes::recipe(data), dn_dc = Inf),
    "dn_dc"
  )
  expect_error(
    step_sec_ri(recipes::recipe(data), instrument_constant = c(1, 2)),
    "instrument_constant"
  )
  expect_error(
    step_sec_uv(recipes::recipe(data), extinction_coef = Inf),
    "extinction_coef"
  )
  expect_error(
    step_sec_uv(recipes::recipe(data), extinction_coef = 0),
    "extinction_coef"
  )
  expect_error(
    step_sec_uv(recipes::recipe(data), wavelength = NA_real_),
    "wavelength"
  )
  expect_error(
    step_sec_lals(
      recipes::recipe(data),
      dn_dc = 0.185,
      calibration_constant = Inf
    ),
    "calibration_constant"
  )
  expect_error(
    step_sec_rals(
      recipes::recipe(data),
      dn_dc = 0.185,
      concentration_unit = "g/L"
    ),
    "concentration_unit"
  )
})

test_that("UV normalization rejects a zero sample-specific extinction value", {
  data <- tibble::tibble(sample_id = "sample", extinction_coef = 0)
  data$uv <- measure::new_measure_list(list(
    measure::new_measure_tbl(location = 1:3, value = c(0, 1e-4, 0))
  ))
  rec <- recipes::recipe(data) |>
    step_sec_uv(extinction_column = "extinction_coef")

  expect_error(
    {
      prepped <- recipes::prep(rec)
      recipes::bake(prepped, new_data = NULL)
    },
    "Each extinction coefficient must be positive and finite"
  )
})

test_that("scattering steps reject concentration unit changes between prep and bake", {
  location <- c(1, 2, 4)
  concentration_mg_per_ml <- c(1, 2, 3)
  training <- tibble::tibble(sample_id = "training")
  training$conc <- measure::new_measure_list(list(
    measure::new_measure_tbl(
      location = location,
      value = concentration_mg_per_ml
    )
  ))
  attr(training$conc, "sec_concentration_unit") <- "mg/mL"
  for (detector in c("lals", "rals", "mals")) {
    training[[detector]] <- measure::new_measure_list(list(
      measure::new_measure_tbl(location = location, value = c(1, 2, 3))
    ))
  }

  bake_data <- tibble::tibble(sample_id = "bake")
  bake_data$conc <- measure::new_measure_list(list(
    measure::new_measure_tbl(
      location = location,
      value = concentration_mg_per_ml / 1000
    )
  ))
  attr(bake_data$conc, "sec_concentration_unit") <- "g/mL"
  for (detector in c("lals", "rals", "mals")) {
    bake_data[[detector]] <- training[[detector]]
  }

  step_builders <- list(
    lals = function(recipe) {
      step_sec_lals(
        recipe,
        measures = "lals",
        concentration_col = "conc",
        dn_dc = 0.185,
        calibration_constant = 1
      )
    },
    rals = function(recipe) {
      step_sec_rals(
        recipe,
        measures = "rals",
        concentration_col = "conc",
        dn_dc = 0.185,
        calibration_constant = 1
      )
    },
    mals = function(recipe) {
      step_sec_mals(
        recipe,
        mals_col = "mals",
        concentration_col = "conc",
        dn_dc = 0.185,
        calibration_constant = 1
      )
    }
  )

  for (builder in step_builders) {
    prepped <- recipes::prep(builder(recipes::recipe(training)))
    expect_error(
      recipes::bake(prepped, new_data = bake_data),
      "conflicts with the concentration curve's recorded units"
    )
  }
})

test_that("non-normalizing concentration conversion preserves recorded units", {
  data <- tibble::tibble(sample_id = "sample")
  data$conc <- measure::new_measure_list(list(
    measure::new_measure_tbl(location = 1:3, value = c(1, 2, 3))
  ))
  attr(data$conc, "sec_concentration_unit") <- "mg/mL"

  out <- recipes::recipe(data) |>
    step_sec_concentration(measures = "conc", normalize_to_mass = FALSE) |>
    recipes::prep() |>
    recipes::bake(new_data = NULL)

  expect_identical(attr(out$conc, "sec_concentration_unit"), "mg/mL")
  expect_equal(out$conc[[1]]$value, c(1, 2, 3))
})

test_that("RI unit metadata survives delay alignment for LALS", {
  location <- c(1, 2, 4)
  dn_dc <- 0.185
  concentration_g_per_ml <- c(0.001, 0.002, 0.003)
  wavelength_nm <- 670
  lambda_cm <- wavelength_nm * 1e-7
  optical_constant <- 4 *
    pi^2 *
    1.333^2 *
    dn_dc^2 /
    (6.022e23 * lambda_cm^4)
  expected_mw <- 50
  data <- tibble::tibble(sample_id = "sample")
  data$ri <- measure::new_measure_list(list(measure::new_measure_tbl(
    location = location,
    value = concentration_g_per_ml * dn_dc
  )))
  data$lals <- measure::new_measure_list(list(measure::new_measure_tbl(
    location = location,
    value = optical_constant * concentration_g_per_ml * expected_mw
  )))

  rec <- recipes::recipe(data) |>
    step_sec_ri(dn_dc = dn_dc) |>
    step_sec_detector_delay(
      reference = "lals",
      targets = "ri",
      delay_volumes = c(ri = 0)
    ) |>
    step_sec_lals(
      measures = "lals",
      concentration_col = "ri",
      dn_dc = dn_dc,
      calibration_constant = 1
    )

  out <- recipes::bake(recipes::prep(rec), new_data = NULL)

  expect_identical(attr(out$ri, "sec_concentration_unit"), "g/mL")
  expect_equal(out$mw_lals[[1]]$value, rep(expected_mw, length(location)))
})

test_that("scattering calculations accept signed nonzero dn/dc", {
  location <- c(1, 2, 4)
  concentration_mg_per_ml <- c(1, 2, 3)
  concentration_g_per_ml <- concentration_mg_per_ml / 1000
  dn_dc <- 0.185
  solvent_ri <- 1.333
  expected_mw <- 50
  wavelength_by_detector <- c(lals = 670, rals = 670, mals = 658)
  optical_constant <- function(wavelength_nm, increment) {
    4 *
      pi^2 *
      solvent_ri^2 *
      increment^2 /
      (6.022e23 * (wavelength_nm * 1e-7)^4)
  }
  data <- tibble::tibble(sample_id = "sample")
  data$conc <- measure::new_measure_list(list(
    measure::new_measure_tbl(
      location = location,
      value = concentration_mg_per_ml
    )
  ))
  attr(data$conc, "sec_concentration_unit") <- "mg/mL"
  for (detector in names(wavelength_by_detector)) {
    data[[detector]] <- measure::new_measure_list(list(
      measure::new_measure_tbl(
        location = location,
        value = optical_constant(wavelength_by_detector[[detector]], dn_dc) *
          concentration_g_per_ml *
          expected_mw
      )
    ))
  }

  build_step <- function(recipe, detector, increment, use_column = FALSE) {
    common <- list(
      recipe = recipe,
      concentration_col = "conc",
      calibration_constant = 1
    )
    if (detector == "lals") {
      common$measures <- "lals"
    } else if (detector == "rals") {
      common$measures <- "rals"
    } else {
      common$mals_col <- "mals"
    }
    if (use_column) {
      common$dn_dc <- NULL
      common$dn_dc_column <- "dn_dc"
    } else {
      common$dn_dc <- increment
    }
    step_fn <- switch(
      detector,
      lals = step_sec_lals,
      rals = step_sec_rals,
      mals = step_sec_mals
    )
    do.call(step_fn, common)
  }

  for (detector in names(wavelength_by_detector)) {
    outputs <- lapply(c(dn_dc, -dn_dc), function(increment) {
      prepped <- recipes::prep(build_step(
        recipes::recipe(data),
        detector,
        increment
      ))
      recipes::bake(prepped, new_data = NULL)[[paste0("mw_", detector)]][[
        1
      ]]$value
    })
    expect_equal(
      outputs[[1]],
      rep(expected_mw, length(location)),
      tolerance = 1e-8
    )
    expect_equal(outputs[[2]], outputs[[1]], tolerance = 1e-8)
  }

  column_data <- data
  column_data$dn_dc <- rep(-dn_dc, nrow(column_data))
  prepped <- recipes::prep(build_step(
    recipes::recipe(column_data),
    "mals",
    increment = NULL,
    use_column = TRUE
  ))
  expect_equal(
    recipes::bake(prepped, new_data = NULL)$mw_mals[[1]]$value,
    rep(expected_mw, length(location)),
    tolerance = 1e-8
  )
})

test_that("scattering rejects zero and nonfinite dn/dc while RI stays positive-only", {
  recipe <- recipes::recipe(tibble::tibble(sample_id = "sample"))
  for (invalid in list(0, NA_real_, Inf, -Inf)) {
    expect_error(step_sec_lals(recipe, dn_dc = invalid), "nonzero finite")
    expect_error(step_sec_rals(recipe, dn_dc = invalid), "nonzero finite")
    expect_error(step_sec_mals(recipe, dn_dc = invalid), "nonzero finite")
  }
  expect_error(step_sec_ri(recipe, dn_dc = -0.185), "positive finite")

  data <- tibble::tibble(sample_id = "sample", dn_dc = 0)
  data$conc <- measure::new_measure_list(list(
    measure::new_measure_tbl(location = 1:2, value = c(1, 2))
  ))
  attr(data$conc, "sec_concentration_unit") <- "mg/mL"
  data$mals <- measure::new_measure_list(list(
    measure::new_measure_tbl(location = 1:2, value = c(1, 2))
  ))
  rec <- recipes::recipe(data) |>
    step_sec_mals(
      mals_col = "mals",
      concentration_col = "conc",
      dn_dc_column = "dn_dc",
      calibration_constant = 1
    )
  expect_error(recipes::prep(rec), "nonzero and finite")
})

test_that("scattering steps preserve inferred unit policy at bake", {
  location <- c(1, 2, 3)
  make_data <- function(concentration, unit = NULL) {
    data <- tibble::tibble(sample_id = "sample")
    data$conc <- measure::new_measure_list(list(
      measure::new_measure_tbl(location, concentration)
    ))
    if (!is.null(unit)) {
      attr(data$conc, "sec_concentration_unit") <- unit
    }
    for (detector in c("lals", "rals", "mals")) {
      data[[detector]] <- measure::new_measure_list(list(
        measure::new_measure_tbl(location, c(1, 2, 3))
      ))
    }
    data
  }

  training <- make_data(c(1, 2, 3), "mg/mL")
  untagged_bake <- make_data(c(0.001, 0.002, 0.003))
  conflicting_bake <- make_data(c(0.001, 0.002, 0.003), "g/mL")
  step_builders <- list(
    lals = function(recipe, concentration_unit = NULL) {
      step_sec_lals(
        recipe,
        measures = "lals",
        concentration_col = "conc",
        concentration_unit = concentration_unit,
        dn_dc = 0.185,
        calibration_constant = 1
      )
    },
    rals = function(recipe, concentration_unit = NULL) {
      step_sec_rals(
        recipe,
        measures = "rals",
        concentration_col = "conc",
        concentration_unit = concentration_unit,
        dn_dc = 0.185,
        calibration_constant = 1
      )
    },
    mals = function(recipe, concentration_unit = NULL) {
      step_sec_mals(
        recipe,
        mals_col = "mals",
        concentration_col = "conc",
        concentration_unit = concentration_unit,
        dn_dc = 0.185,
        calibration_constant = 1
      )
    }
  )

  for (builder in step_builders) {
    prepped <- recipes::prep(builder(recipes::recipe(training)))
    serialized_recipe <- withr::local_tempfile()
    saveRDS(prepped, serialized_recipe)
    prepped <- readRDS(serialized_recipe)

    expect_error(
      recipes::bake(prepped, new_data = untagged_bake),
      "no longer records its unit"
    )

    explicit_prepped <- recipes::prep(
      builder(recipes::recipe(training), concentration_unit = "mg/mL")
    )
    expect_no_error(
      recipes::bake(explicit_prepped, new_data = untagged_bake)
    )
    expect_error(
      recipes::bake(explicit_prepped, new_data = conflicting_bake),
      "conflicts with the concentration curve's recorded units"
    )
  }
})
