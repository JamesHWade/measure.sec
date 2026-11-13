sec_elution_concentration_curve <- function(location, value) {
  measure::new_measure_tbl(location = location, value = value)
}

sec_elution_concentration_data <- function(curves, id = seq_along(curves)) {
  data <- tibble::tibble(sample_id = id)
  data$ri <- measure::new_measure_list(curves)
  data
}

sec_elution_concentration_recipe <- function(
  training,
  elution_units = "min",
  flow_rate = 0.5,
  concentration_units = "mg/mL",
  injection_mass = 2,
  normalize_to_mass = TRUE
) {
  recipes::recipe(training) |>
    step_sec_concentration(
      measures = "ri",
      injection_mass = injection_mass,
      flow_rate = flow_rate,
      concentration_units = concentration_units,
      normalize_to_mass = normalize_to_mass,
      elution_units = elution_units
    )
}

sec_elution_trapezoid <- function(location, value) {
  sum(abs(diff(location)) * (head(value, -1) + tail(value, -1)) / 2)
}

sec_elution_capture_warning <- function(expr) {
  captured <- new.env(parent = emptyenv())
  captured$warnings <- list()
  value <- withCallingHandlers(
    expr,
    measure_sec_negative_signal_clipped = function(warning) {
      captured$warnings[[length(captured$warnings) + 1L]] <- warning
      invokeRestart("muffleWarning")
    }
  )

  list(value = value, warnings = captured$warnings)
}

test_that("concentration normalization handles time and volume axes", {
  time_min <- c(0, 1, 2, 4)
  signal <- c(0, 0.001, 0.002, 0)
  axes <- list(
    min = list(location = time_min, flow = 0.5),
    s = list(location = time_min * 60, flow = 0.5),
    mL = list(location = time_min * 0.5, flow = 8.25),
    uL = list(location = time_min * 500, flow = NULL)
  )

  for (output_units in c("mg/mL", "g/mL")) {
    outputs <- lapply(names(axes), function(units) {
      axis <- axes[[units]]
      training <- sec_elution_concentration_data(list(
        sec_elution_concentration_curve(axis$location, signal)
      ))
      prepped <- recipes::prep(sec_elution_concentration_recipe(
        training,
        elution_units = units,
        flow_rate = axis$flow,
        concentration_units = output_units
      ))
      baked <- recipes::bake(prepped, new_data = NULL)
      curve <- baked$ri[[1]]
      volume_factor <- switch(
        units,
        min = axis$flow,
        s = axis$flow / 60,
        mL = 1,
        uL = 0.001
      )
      mass_mg <- sec_elution_trapezoid(curve$location, curve$value) *
        volume_factor *
        if (identical(output_units, "g/mL")) 1000 else 1

      expect_identical(curve$location, axis$location)
      expect_identical(
        attr(baked$ri, "sec_concentration_unit"),
        output_units
      )
      expect_equal(mass_mg, 2, tolerance = 1e-12)
      expect_equal(
        curve$value,
        if (identical(output_units, "mg/mL")) {
          c(0, 1, 2, 0)
        } else {
          signal
        },
        tolerance = 1e-12
      )
      curve$value
    })

    for (result in outputs[-1]) {
      expect_equal(result, outputs[[1]], tolerance = 1e-12)
    }
  }

  named_unit_training <- sec_elution_concentration_data(list(
    sec_elution_concentration_curve(time_min, signal)
  ))
  named_unit_recipe <- sec_elution_concentration_recipe(
    named_unit_training,
    concentration_units = structure("g/mL", names = "output")
  )
  expect_identical(
    named_unit_recipe$steps[[1]]$concentration_units,
    "g/mL"
  )
  named_unit_output <- recipes::bake(
    recipes::prep(named_unit_recipe),
    new_data = NULL
  )
  expect_equal(named_unit_output$ri[[1]]$value, signal, tolerance = 1e-12)
  expect_identical(
    attr(named_unit_output$ri, "sec_concentration_unit"),
    "g/mL"
  )

  volume_axis <- time_min * 0.5
  volume_outputs <- lapply(c(0.12, 17.5, 1e308), function(flow) {
    training <- sec_elution_concentration_data(list(
      sec_elution_concentration_curve(volume_axis, signal)
    ))
    prepped <- recipes::prep(sec_elution_concentration_recipe(
      training,
      elution_units = "mL",
      flow_rate = flow
    ))
    recipes::bake(prepped, new_data = NULL)$ri[[1]]$value
  })
  expect_identical(volume_outputs[[1]], volume_outputs[[2]])
  expect_identical(volume_outputs[[1]], volume_outputs[[3]])
})

test_that("concentration baking preserves raw coordinates and current row order", {
  training <- sec_elution_concentration_data(list(
    sec_elution_concentration_curve(c(0, 1, 2, 4), c(0, 0.001, 0.002, 0))
  ))
  prepped <- recipes::prep(sec_elution_concentration_recipe(training))

  descending_location <- c(4, 2, 1, 0)
  descending_value <- c(0, 0.002, 0.001, 0)
  irregular_location <- c(0, 0.5, 2, 4)
  irregular_value <- c(0, 0.001, 0.003, 0.0005)
  baked <- recipes::bake(
    prepped,
    new_data = sec_elution_concentration_data(
      list(
        sec_elution_concentration_curve(
          descending_location,
          descending_value
        ),
        sec_elution_concentration_curve(
          irregular_location,
          irregular_value
        )
      ),
      id = c("descending", "irregular")
    )
  )

  irregular_area <- 0.5 *
    (0 + 0.001) /
    2 +
    1.5 * (0.001 + 0.003) / 2 +
    2 * (0.003 + 0.0005) / 2
  expected_irregular <- irregular_value * 2 / (irregular_area * 0.5)

  expect_identical(as.character(baked$sample_id), c("descending", "irregular"))
  expect_identical(baked$ri[[1]]$location, descending_location)
  expect_identical(baked$ri[[2]]$location, irregular_location)
  expect_equal(baked$ri[[1]]$value, c(0, 2, 1, 0), tolerance = 1e-12)
  expect_equal(baked$ri[[2]]$value, expected_irregular, tolerance = 1e-12)
  expect_equal(
    sec_elution_trapezoid(baked$ri[[1]]$location, baked$ri[[1]]$value) * 0.5,
    2,
    tolerance = 1e-12
  )
  expect_equal(
    sec_elution_trapezoid(baked$ri[[2]]$location, baked$ri[[2]]$value) * 0.5,
    2,
    tolerance = 1e-12
  )
})

test_that("minute defaults and historical serialized steps stay compatible", {
  location <- c(0, 1, 3)
  signal <- c(-2, 2, 0)
  training <- sec_elution_concentration_data(list(
    sec_elution_concentration_curve(location, c(1, 2, 1))
  ))
  omitted <- recipes::prep(
    recipes::recipe(training) |>
      step_sec_concentration(
        measures = "ri",
        injection_mass = 6,
        flow_rate = 0.5
      )
  )
  explicit <- recipes::prep(sec_elution_concentration_recipe(
    training,
    elution_units = "min",
    flow_rate = 0.5,
    injection_mass = 6
  ))
  input <- sec_elution_concentration_data(list(
    sec_elution_concentration_curve(location, signal)
  ))
  omitted_result <- sec_elution_capture_warning(
    recipes::bake(omitted, new_data = input)
  )
  explicit_result <- sec_elution_capture_warning(
    recipes::bake(explicit, new_data = input)
  )

  expect_identical(omitted_result$value, explicit_result$value)
  expect_identical(
    omitted_result$warnings[[1]]$diagnostics,
    explicit_result$warnings[[1]]$diagnostics
  )
  expect_identical(omitted_result$warnings[[1]]$elution_units, "min")
  expect_identical(explicit_result$warnings[[1]]$elution_units, "min")

  legacy_recipe <- recipes::recipe(training) |>
    step_sec_concentration(
      measures = "ri",
      injection_mass = 6,
      flow_rate = 0.5
    )
  legacy_recipe$steps[[1]]$elution_units <- NULL
  legacy_recipe$steps[[1]]$elution_units_backup <- "s"
  legacy_prepped <- recipes::prep(legacy_recipe)
  expect_identical(
    recipes::tidy(legacy_prepped, number = 1)$elution_units,
    "min"
  )

  legacy_prepped$steps[[1]]$elution_units <- NULL
  path <- withr::local_tempfile(fileext = ".rds")
  saveRDS(legacy_prepped, path)
  legacy_prepped <- readRDS(path)
  legacy_baked <- sec_elution_capture_warning(
    recipes::bake(legacy_prepped, new_data = input)
  )
  expect_identical(legacy_baked$value, explicit_result$value)
})

test_that("clipping warnings identify raw axis units without changing table shape", {
  base_location <- c(0, 1, 3)
  signal <- c(-2, 2, 0)
  axes <- list(
    min = list(
      location = base_location,
      flow = 0.5,
      positive = 3,
      negative = 1
    ),
    s = list(
      location = base_location * 60,
      flow = 0.5,
      positive = 180,
      negative = 60
    ),
    mL = list(
      location = base_location * 0.5,
      flow = NULL,
      positive = 1.5,
      negative = 0.5
    ),
    uL = list(
      location = base_location * 500,
      flow = 7,
      positive = 1500,
      negative = 500
    )
  )
  training <- sec_elution_concentration_data(list(
    sec_elution_concentration_curve(base_location, c(1, 2, 1))
  ))
  warnings <- lapply(names(axes), function(units) {
    axis <- axes[[units]]
    prepped <- recipes::prep(sec_elution_concentration_recipe(
      training,
      elution_units = units,
      flow_rate = axis$flow,
      injection_mass = 6
    ))
    result <- sec_elution_capture_warning(recipes::bake(
      prepped,
      new_data = sec_elution_concentration_data(list(
        sec_elution_concentration_curve(axis$location, signal)
      ))
    ))
    warning <- result$warnings[[1]]

    expect_length(result$warnings, 1)
    expect_s3_class(warning, "measure_sec_negative_signal_clipped")
    expect_identical(warning$column, "ri")
    expect_identical(warning$elution_units, units)
    expect_named(
      warning$diagnostics,
      c(
        "row",
        "negative_points",
        "points",
        "positive_area",
        "negative_area",
        "excluded_fraction"
      )
    )
    expect_identical(warning$diagnostics$row, 1L)
    expect_identical(warning$diagnostics$negative_points, 1L)
    expect_identical(warning$diagnostics$points, 3L)
    expect_identical(warning$diagnostics$positive_area, axis$positive)
    expect_identical(warning$diagnostics$negative_area, axis$negative)
    expect_identical(warning$diagnostics$excluded_fraction, 0.25)
    expect_equal(result$value$ri[[1]]$value, c(0, 8, 0))

    warning
  })

  expect_snapshot({
    for (warning in warnings) {
      cat(warning$elution_units, ": ", sep = "")
      cat(
        warning$diagnostics$positive_area,
        "/",
        warning$diagnostics$negative_area,
        "\n",
        sep = ""
      )
      cat(conditionMessage(warning), "\n", sep = "")
    }
  })

  empty_training <- sec_elution_concentration_data(list(
    sec_elution_concentration_curve(base_location, c(1, 2, 1))
  ))
  prepped <- recipes::prep(sec_elution_concentration_recipe(empty_training))
  empty_data <- sec_elution_concentration_data(list())
  empty_result <- recipes::bake(prepped, new_data = empty_data)

  expect_identical(nrow(empty_result), 0L)
  expect_identical(attr(empty_result$ri, "sec_concentration_unit"), "mg/mL")
  empty_diagnostics <- .sec_concentration_diagnostics(list())
  expect_identical(empty_diagnostics$row, integer())
  expect_identical(empty_diagnostics$negative_points, integer())
  expect_identical(empty_diagnostics$points, integer())
  expect_identical(empty_diagnostics$positive_area, double())
  expect_identical(empty_diagnostics$negative_area, double())
  expect_identical(empty_diagnostics$excluded_fraction, double())
})

test_that("elution units and flow policies validate across step methods", {
  training <- sec_elution_concentration_data(list(
    sec_elution_concentration_curve(c(0, 1), c(1, 1))
  ))
  make_step <- function(...) {
    recipes::recipe(training) |>
      step_sec_concentration(
        measures = "ri",
        injection_mass = 1,
        ...
      )
  }

  expect_identical(.sec_elution_units(structure("mL", names = "unit")), "mL")
  named_units <- make_step(
    elution_units = structure("mL", names = "unit"),
    flow_rate = NULL
  )
  expect_identical(named_units$steps[[1]]$elution_units, "mL")
  expect_error(make_step(elution_units = NULL), "elution_units")
  expect_error(make_step(elution_units = NA_character_), "elution_units")
  expect_error(make_step(elution_units = c("min", "s")), "elution_units")
  expect_error(make_step(elution_units = factor("min")), "elution_units")
  expect_error(make_step(elution_units = "minutes"), "elution_units")
  expect_error(
    make_step(elution_units = structure("min", class = "axis")),
    "elution_units"
  )
  expect_error(make_step(elution_units = matrix("min")), "elution_units")

  expect_error(make_step(elution_units = "min", flow_rate = NULL), "flow_rate")
  expect_error(make_step(elution_units = "s", flow_rate = NULL), "flow_rate")
  expect_no_error(make_step(elution_units = "mL", flow_rate = NULL))
  expect_no_error(make_step(elution_units = "uL", flow_rate = NULL))
  expect_error(make_step(elution_units = "mL", flow_rate = 0), "flow_rate")
  expect_error(make_step(elution_units = "mL", flow_rate = -1), "flow_rate")
  expect_error(make_step(elution_units = "mL", flow_rate = Inf), "flow_rate")
  expect_error(
    make_step(elution_units = "mL", flow_rate = c(1, 2)),
    "flow_rate"
  )
  expect_error(
    make_step(elution_units = "mL", flow_rate = factor(1)),
    "flow_rate"
  )
  expect_error(
    recipes::recipe(training) |>
      step_sec_concentration(
        measures = "ri",
        injection_mass = matrix(1),
        elution_units = "mL"
      ),
    "positive finite mass"
  )
  expect_error(
    recipes::recipe(training) |>
      step_sec_concentration(
        measures = "ri",
        injection_volume = matrix(100),
        sample_concentration = 2,
        elution_units = "mL"
      ),
    "injection_volume.*plain numeric scalar"
  )
  expect_error(
    recipes::recipe(training) |>
      step_sec_concentration(
        measures = "ri",
        injection_volume = 100,
        sample_concentration = matrix(2),
        elution_units = "mL"
      ),
    "sample_concentration.*plain numeric scalar"
  )
  named_mass <- recipes::prep(
    recipes::recipe(training) |>
      step_sec_concentration(
        measures = "ri",
        injection_mass = structure(1, names = "mass"),
        flow_rate = NULL,
        elution_units = "mL"
      )
  )
  named_mass_output <- recipes::bake(named_mass, new_data = NULL)$ri[[1]]$value
  expect_identical(length(named_mass_output), 2L)
  expect_false(is.matrix(named_mass_output))
  expect_error(
    make_step(elution_units = "s", flow_rate = 5e-324),
    "conversion factor"
  )

  tagged <- sec_elution_concentration_data(list(
    sec_elution_concentration_curve(c(0, 1, 2), c(NA_real_, -1, 2))
  ))
  attr(tagged$ri, "sec_concentration_unit") <- "g/mL"
  no_op <- recipes::prep(
    recipes::recipe(tagged) |>
      step_sec_concentration(
        measures = "ri",
        normalize_to_mass = FALSE,
        flow_rate = NULL,
        elution_units = "mL"
      )
  )
  expect_identical(recipes::bake(no_op, new_data = tagged), tagged)
  expect_error(
    recipes::recipe(tagged) |>
      step_sec_concentration(
        measures = "ri",
        normalize_to_mass = FALSE,
        flow_rate = 0
      ),
    "flow_rate"
  )

  invalid_recipe <- make_step(elution_units = "mL", flow_rate = NULL)
  invalid_recipe$steps[[1]]["elution_units"] <- list(factor("mL"))
  expect_error(recipes::prep(invalid_recipe), "elution_units")

  invalid_flow_recipe <- make_step(elution_units = "mL", flow_rate = NULL)
  invalid_flow_recipe$steps[[1]]$flow_rate <- 0
  expect_error(recipes::prep(invalid_flow_recipe), "flow_rate")

  invalid_prepped <- recipes::prep(make_step(
    elution_units = "mL",
    flow_rate = NULL
  ))
  invalid_prepped$steps[[1]]["elution_units"] <- list(factor("mL"))
  expect_error(
    recipes::bake(invalid_prepped, new_data = training),
    "elution_units"
  )
  expect_error(
    recipes::tidy(invalid_prepped, number = 1),
    "elution_units"
  )
  invalid_prepped$steps[[1]]["elution_units"] <- list("mL")
  invalid_prepped$steps[[1]]$flow_rate <- 0
  expect_error(recipes::bake(invalid_prepped, new_data = training), "flow_rate")
  expect_error(recipes::tidy(invalid_prepped, number = 1), "flow_rate")

  other_training <- tibble::tibble(sample_id = "training")
  other_training$uv <- measure::new_measure_list(list(
    sec_elution_concentration_curve(c(0, 1), c(1, 1))
  ))
  missing_prep <- recipes::recipe(other_training) |>
    step_sec_concentration(
      measures = "ri",
      injection_mass = 1,
      elution_units = "mL"
    )
  expect_error(recipes::prep(missing_prep), "must exist as measure columns")

  normalized_prepped <- recipes::prep(make_step(elution_units = "mL"))
  missing_bake_data <- tibble::tibble(sample_id = "new")
  expect_error(
    bake.step_sec_concentration(
      normalized_prepped$steps[[1]],
      new_data = missing_bake_data
    ),
    "must exist as measure columns"
  )
  pass_through <- recipes::prep(
    recipes::recipe(training) |>
      step_sec_concentration(
        measures = "ri",
        normalize_to_mass = FALSE,
        flow_rate = NULL,
        elution_units = "mL"
      )
  )
  expect_identical(
    bake.step_sec_concentration(
      pass_through$steps[[1]],
      new_data = missing_bake_data
    ),
    missing_bake_data
  )
})

test_that("output units and resolved mass stay valid through step lifecycles", {
  training <- sec_elution_concentration_data(list(
    sec_elution_concentration_curve(c(0, 1, 2, 4), c(0, 0.001, 0.002, 0))
  ))
  make_step <- function(injection_mass = 2, ...) {
    recipes::recipe(training) |>
      step_sec_concentration(
        measures = "ri",
        injection_mass = injection_mass,
        flow_rate = 0.5,
        elution_units = "min",
        ...
      )
  }

  named_unit_recipe <- make_step(
    concentration_units = structure("g/mL", names = "output")
  )
  expect_identical(
    named_unit_recipe$steps[[1]]$concentration_units,
    "g/mL"
  )
  named_unit_output <- recipes::bake(
    recipes::prep(named_unit_recipe),
    new_data = NULL
  )
  expect_equal(
    named_unit_output$ri[[1]]$value,
    c(0, 0.001, 0.002, 0),
    tolerance = 1e-12
  )
  expect_identical(
    attr(named_unit_output$ri, "sec_concentration_unit"),
    "g/mL"
  )

  expect_error(make_step(concentration_units = NULL), "concentration_units")
  expect_error(
    make_step(concentration_units = NA_character_),
    "concentration_units"
  )
  expect_error(
    make_step(concentration_units = c("g/mL", "mg/mL")),
    "concentration_units"
  )
  expect_error(
    make_step(concentration_units = factor("g/mL")),
    "concentration_units"
  )
  expect_error(
    make_step(concentration_units = matrix("g/mL")),
    "concentration_units"
  )
  expect_error(
    make_step(concentration_units = structure("g/mL", class = "unit")),
    "concentration_units"
  )

  invalid_prep_unit <- make_step()
  invalid_prep_unit$steps[[1]]["concentration_units"] <- list(NULL)
  expect_error(recipes::prep(invalid_prep_unit), "concentration_units")

  invalid_prepped_unit <- recipes::prep(make_step())
  invalid_prepped_unit$steps[[1]]$concentration_units <- NA_character_
  expect_error(
    recipes::bake(invalid_prepped_unit, new_data = training),
    "concentration_units"
  )
  invalid_prepped_unit$steps[[1]]["concentration_units"] <- list(NULL)
  expect_error(
    recipes::tidy(invalid_prepped_unit, number = 1),
    "concentration_units"
  )

  expect_error(
    make_step(injection_mass = matrix(2)),
    "positive finite mass"
  )
  expect_error(
    make_step(injection_mass = structure(2, class = "mass")),
    "positive finite mass"
  )

  invalid_prep_mass <- make_step()
  invalid_prep_mass$steps[[1]]$injection_mass <- NULL
  expect_error(recipes::prep(invalid_prep_mass), "positive finite mass")

  invalid_prepped_mass <- recipes::prep(make_step())
  invalid_prepped_mass$steps[[1]]$injection_mass <- NULL
  expect_error(
    recipes::bake(invalid_prepped_mass, new_data = training),
    "positive finite mass"
  )
  expect_error(
    recipes::tidy(invalid_prepped_mass, number = 1),
    "positive finite mass"
  )

  pass_through <- recipes::prep(
    recipes::recipe(training) |>
      step_sec_concentration(
        measures = "ri",
        injection_mass = matrix(2),
        normalize_to_mass = FALSE,
        flow_rate = NULL,
        concentration_units = structure("g/mL", names = "output"),
        elution_units = "mL"
      )
  )
  expect_identical(recipes::bake(pass_through, new_data = training), training)
})

test_that("concentration normalization rejects malformed and unrepresentable curves", {
  malformed_curve <- function(location, value) {
    structure(
      list(location = location, value = value),
      names = c("location", "value"),
      row.names = c(NA_integer_, -length(location)),
      class = c("measure_tbl", "tbl_df", "tbl", "data.frame")
    )
  }
  normalize_curve <- function(curve) {
    .sec_normalize_concentration_measure(
      m = curve,
      injection_mass = 2,
      volume_factor = 0.5,
      concentration_units = "mg/mL",
      row = 1
    )
  }

  expect_error(
    normalize_curve(malformed_curve(c(0, 1), c(1, 2, 3))),
    "equal length"
  )
  expect_error(
    normalize_curve(malformed_curve(c(0, 1), matrix(c(1, 2), ncol = 1))),
    "plain numeric one-dimensional"
  )
  expect_error(
    normalize_curve(malformed_curve(
      c(0, 1),
      structure(c(1, 2), class = "signal")
    )),
    "plain numeric one-dimensional"
  )
  expect_error(
    normalize_curve(malformed_curve(c(0, 1, 1), c(1, 2, 3))),
    "monotonic"
  )
  expect_error(
    normalize_curve(malformed_curve(c(0, Inf), c(1, 2))),
    "monotonic"
  )
  expect_error(
    normalize_curve(malformed_curve(c(0, 1), c(1, NaN))),
    "monotonic"
  )
  expect_error(
    normalize_curve(malformed_curve(c(-1e308, 1e308), c(1, 1))),
    "axis widths"
  )
})
