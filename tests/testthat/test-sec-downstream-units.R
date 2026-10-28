test_that("intrinsic viscosity converts tagged concentration to physical units", {
  location <- c(1, 2, 3)
  specific_viscosity <- c(0.05, 0.05, 0.05)

  make_result <- function(concentration, unit, output_units = "dL/g") {
    data <- tibble::tibble(sample_id = "sample")
    data$visc <- measure::new_measure_list(
      list(measure::new_measure_tbl(
        location = location,
        value = specific_viscosity
      ))
    )
    data$conc <- measure::new_measure_list(
      list(measure::new_measure_tbl(
        location = location,
        value = rep(concentration, length(location))
      ))
    )
    attr(data$conc, "sec_concentration_unit") <- unit

    recipe <- recipes::recipe(~., data = data) |>
      step_sec_intrinsic_visc(
        specific_visc_col = "visc",
        concentration_col = "conc",
        units = output_units
      )

    recipes::bake(recipes::prep(recipe), new_data = NULL)
  }

  result_mg <- make_result(2, "mg/mL")
  result_g <- make_result(0.002, "g/mL")
  result_ml_per_g <- make_result(0.002, "g/mL", "mL/g")

  expect_equal(result_mg$intrinsic_visc[[1]]$value, rep(0.25, 3))
  expect_equal(result_g$intrinsic_visc[[1]]$value, rep(0.25, 3))
  expect_equal(result_ml_per_g$intrinsic_visc[[1]]$value, rep(25, 3))
  expect_identical(
    attr(result_g$intrinsic_visc, "sec_intrinsic_viscosity_unit"),
    "dL/g"
  )
  expect_identical(
    attr(result_ml_per_g$intrinsic_visc, "sec_intrinsic_viscosity_unit"),
    "mL/g"
  )
})

test_that("intrinsic viscosity requires or records an untagged unit contract", {
  location <- c(1, 2, 3)
  data <- tibble::tibble(sample_id = "sample")
  data$visc <- measure::new_measure_list(
    list(measure::new_measure_tbl(location, c(0.05, 0.05, 0.05)))
  )
  data$conc <- measure::new_measure_list(
    list(measure::new_measure_tbl(location, c(2, 2, 2)))
  )

  recipe <- recipes::recipe(~., data = data) |>
    step_sec_intrinsic_visc(
      specific_visc_col = "visc",
      concentration_col = "conc"
    )

  expect_error(
    recipes::prep(recipe),
    "Specify.*concentration_unit.*untagged"
  )

  explicit_recipe <- recipes::recipe(~., data = data) |>
    step_sec_intrinsic_visc(
      specific_visc_col = "visc",
      concentration_col = "conc",
      concentration_unit = "mg/mL"
    )
  result <- recipes::bake(recipes::prep(explicit_recipe), new_data = NULL)

  expect_equal(result$intrinsic_visc[[1]]$value, rep(0.25, 3))
})

test_that("intrinsic viscosity rejects concentration unit drift at bake", {
  location <- c(1, 2, 3)
  make_data <- function(concentration, unit) {
    data <- tibble::tibble(sample_id = "sample")
    data$visc <- measure::new_measure_list(
      list(measure::new_measure_tbl(location, c(0.05, 0.05, 0.05)))
    )
    data$conc <- measure::new_measure_list(
      list(measure::new_measure_tbl(
        location,
        rep(concentration, length(location))
      ))
    )
    attr(data$conc, "sec_concentration_unit") <- unit
    data
  }

  recipe <- recipes::recipe(~., data = make_data(2, "mg/mL")) |>
    step_sec_intrinsic_visc(
      specific_visc_col = "visc",
      concentration_col = "conc"
    )

  prepped <- recipes::prep(recipe)

  expect_error(
    recipes::bake(
      prepped,
      new_data = make_data(0.002, "g/mL")
    ),
    "conflicts with the concentration curve's recorded units"
  )
})

test_that("intrinsic viscosity requires aligned curves and positive inputs", {
  location <- c(1, 2, 3)
  data <- tibble::tibble(sample_id = "sample")
  data$visc <- measure::new_measure_list(
    list(measure::new_measure_tbl(location, c(0.05, 0, -0.05)))
  )
  data$conc <- measure::new_measure_list(
    list(measure::new_measure_tbl(location, c(2, 2, -2)))
  )
  attr(data$conc, "sec_concentration_unit") <- "mg/mL"

  recipe <- recipes::recipe(~., data = data) |>
    step_sec_intrinsic_visc(
      specific_visc_col = "visc",
      concentration_col = "conc"
    )
  result <- recipes::bake(recipes::prep(recipe), new_data = NULL)

  expect_equal(result$intrinsic_visc[[1]]$value, c(0.25, NA_real_, NA_real_))

  misaligned <- data
  misaligned$conc <- measure::new_measure_list(
    list(measure::new_measure_tbl(location + 0.1, c(2, 2, 2)))
  )
  attr(misaligned$conc, "sec_concentration_unit") <- "mg/mL"
  misaligned_recipe <- recipes::recipe(~., data = misaligned) |>
    step_sec_intrinsic_visc(
      specific_visc_col = "visc",
      concentration_col = "conc"
    )

  expect_error(
    recipes::bake(recipes::prep(misaligned_recipe), new_data = NULL),
    "locations must align"
  )
})

test_that("composition solves the calibrated two-component response matrix", {
  location <- c(1, 2, 3)
  uv_value <- rep(0.0035, 3)
  ri_value <- rep(0.01, 3)
  data <- tibble::tibble(sample_id = "sample")
  data$uv <- measure::new_measure_list(
    list(measure::new_measure_tbl(location, uv_value))
  )
  data$ri <- measure::new_measure_list(
    list(measure::new_measure_tbl(location, ri_value))
  )

  recipe <- recipes::recipe(~., data = data) |>
    step_sec_composition(
      uv_col = "uv",
      ri_col = "ri",
      component_a_uv = 2,
      component_a_ri = 1,
      component_b_uv = 0.5,
      component_b_ri = 3
    )
  result <- recipes::bake(recipes::prep(recipe), new_data = NULL)

  expect_equal(result$composition_a[[1]]$value, rep(0.25, 3))
})

test_that("composition normalizes tagged or explicitly declared channels", {
  location <- c(1, 2, 3)
  make_data <- function(uv, ri, uv_unit = NULL, ri_unit = NULL) {
    data <- tibble::tibble(sample_id = "sample")
    data$uv <- measure::new_measure_list(
      list(measure::new_measure_tbl(location, uv))
    )
    data$ri <- measure::new_measure_list(
      list(measure::new_measure_tbl(location, ri))
    )
    if (!is.null(uv_unit)) {
      attr(data$uv, "sec_concentration_unit") <- uv_unit
    }
    if (!is.null(ri_unit)) {
      attr(data$ri, "sec_concentration_unit") <- ri_unit
    }
    data
  }
  make_recipe <- function(data, uv_unit = NULL, ri_unit = NULL) {
    recipes::recipe(~., data = data) |>
      step_sec_composition(
        uv_col = "uv",
        ri_col = "ri",
        component_a_uv = 2,
        component_a_ri = 1,
        component_b_uv = 0.5,
        component_b_ri = 3,
        input_basis = "normalized_concentration",
        uv_concentration_unit = uv_unit,
        ri_concentration_unit = ri_unit
      )
  }

  result_mixed_units <- recipes::bake(
    recipes::prep(make_recipe(
      make_data(rep(3.5, 3), rep(0.01, 3), "mg/mL", "g/mL")
    )),
    new_data = NULL
  )
  result_g_per_ml <- recipes::bake(
    recipes::prep(make_recipe(
      make_data(rep(0.0035, 3), rep(0.01, 3), "g/mL", "g/mL")
    )),
    new_data = NULL
  )
  result_explicit_units <- recipes::bake(
    recipes::prep(make_recipe(
      make_data(rep(3.5, 3), rep(0.01, 3)),
      uv_unit = "mg/mL",
      ri_unit = "g/mL"
    )),
    new_data = NULL
  )

  expect_equal(result_mixed_units$composition_a[[1]]$value, rep(0.25, 3))
  expect_equal(
    result_g_per_ml$composition_a[[1]]$value,
    result_mixed_units$composition_a[[1]]$value
  )
  expect_equal(
    result_explicit_units$composition_a[[1]]$value,
    result_mixed_units$composition_a[[1]]$value
  )

  training_recipe <- make_recipe(
    make_data(rep(3.5, 3), rep(0.01, 3), "mg/mL", "g/mL")
  )
  prepped <- recipes::prep(training_recipe)
  changed_units <- make_data(
    rep(0.0035, 3),
    rep(0.01, 3),
    "g/mL",
    "g/mL"
  )
  expect_error(
    recipes::bake(prepped, new_data = changed_units),
    "conflicts with the concentration curve's recorded units"
  )

  untagged <- make_data(rep(3.5, 3), rep(0.01, 3))
  expect_error(
    recipes::prep(make_recipe(untagged)),
    "Specify.*uv_concentration_unit.*untagged"
  )
})

test_that("composition accepts zero and signed factors when identifiable", {
  location <- c(1, 2, 3)
  data <- tibble::tibble(sample_id = "sample")
  data$uv <- measure::new_measure_list(
    list(measure::new_measure_tbl(location, rep(0.003, 3)))
  )
  data$ri <- measure::new_measure_list(
    list(measure::new_measure_tbl(location, rep(0.007, 3)))
  )

  zero_factor_recipe <- recipes::recipe(~., data = data) |>
    step_sec_composition(
      uv_col = "uv",
      ri_col = "ri",
      component_a_uv = 0,
      component_a_ri = 1,
      component_b_uv = 1,
      component_b_ri = 2
    )
  zero_factor_result <- recipes::bake(
    recipes::prep(zero_factor_recipe),
    new_data = NULL
  )

  signed_data <- data
  signed_data$uv <- measure::new_measure_list(
    list(measure::new_measure_tbl(location, rep(0.0035, 3)))
  )
  signed_data$ri <- measure::new_measure_list(
    list(measure::new_measure_tbl(location, rep(-0.008, 3)))
  )
  signed_factor_recipe <- recipes::recipe(~., data = signed_data) |>
    step_sec_composition(
      uv_col = "uv",
      ri_col = "ri",
      component_a_uv = 2,
      component_a_ri = 1,
      component_b_uv = 0.5,
      component_b_ri = -3
    )
  signed_factor_result <- recipes::bake(
    recipes::prep(signed_factor_recipe),
    new_data = NULL
  )

  expect_equal(zero_factor_result$composition_a[[1]]$value, rep(0.25, 3))
  expect_equal(signed_factor_result$composition_a[[1]]$value, rep(0.25, 3))
  expect_error(
    step_sec_composition(
      recipes::recipe(~., data = data),
      component_a_uv = 1,
      component_a_ri = 2,
      component_b_uv = 2,
      component_b_ri = 4
    ),
    "cannot identify composition"
  )
})

test_that("composition enforces its input basis and curve alignment", {
  location <- c(1, 2, 3)
  data <- tibble::tibble(sample_id = "sample")
  data$uv <- measure::new_measure_list(
    list(measure::new_measure_tbl(location, rep(0.0035, 3)))
  )
  data$ri <- measure::new_measure_list(
    list(measure::new_measure_tbl(location, rep(0.01, 3)))
  )
  tagged_data <- data
  attr(tagged_data$uv, "sec_concentration_unit") <- "g/mL"
  expect_error(
    recipes::prep(
      recipes::recipe(~., data = tagged_data) |>
        step_sec_composition(
          uv_col = "uv",
          ri_col = "ri",
          component_a_uv = 2,
          component_a_ri = 1,
          component_b_uv = 0.5,
          component_b_ri = 3
        )
    ),
    "cannot use concentration-tagged curves"
  )

  data$uv <- measure::new_measure_list(
    list(measure::new_measure_tbl(location + 0.1, rep(0.0035, 3)))
  )
  misaligned_recipe <- recipes::recipe(~., data = data) |>
    step_sec_composition(
      uv_col = "uv",
      ri_col = "ri",
      component_a_uv = 2,
      component_a_ri = 1,
      component_b_uv = 0.5,
      component_b_ri = 3
    )

  expect_error(
    recipes::bake(recipes::prep(misaligned_recipe), new_data = NULL),
    "locations must align"
  )
})

test_that("UV/RI ratio converts concentration channels to common units", {
  location <- c(1, 2, 3)
  make_data <- function(uv, ri, uv_unit = NULL, ri_unit = NULL) {
    data <- tibble::tibble(sample_id = "sample")
    data$uv <- measure::new_measure_list(
      list(measure::new_measure_tbl(location, uv))
    )
    data$ri <- measure::new_measure_list(
      list(measure::new_measure_tbl(location, ri))
    )
    if (!is.null(uv_unit)) {
      attr(data$uv, "sec_concentration_unit") <- uv_unit
    }
    if (!is.null(ri_unit)) {
      attr(data$ri, "sec_concentration_unit") <- ri_unit
    }
    data
  }
  make_result <- function(
    data,
    input_basis = "auto",
    uv_unit = NULL,
    ri_unit = NULL
  ) {
    recipe <- recipes::recipe(~., data = data) |>
      step_sec_uv_ri_ratio(
        uv_col = "uv",
        ri_col = "ri",
        smooth = FALSE,
        input_basis = input_basis,
        uv_concentration_unit = uv_unit,
        ri_concentration_unit = ri_unit
      )
    recipes::bake(recipes::prep(recipe), new_data = NULL)
  }

  result_mixed_units <- make_result(
    make_data(rep(4, 3), rep(0.002, 3), "mg/mL", "g/mL")
  )
  result_g_per_ml <- make_result(
    make_data(rep(0.004, 3), rep(0.002, 3), "g/mL", "g/mL")
  )
  result_explicit_units <- make_result(
    make_data(rep(4, 3), rep(0.002, 3)),
    input_basis = "normalized_concentration",
    uv_unit = "mg/mL",
    ri_unit = "g/mL"
  )

  expect_equal(result_mixed_units$uv_ri_ratio[[1]]$value, rep(2, 3))
  expect_equal(
    result_g_per_ml$uv_ri_ratio[[1]]$value,
    result_mixed_units$uv_ri_ratio[[1]]$value
  )
  expect_equal(
    result_explicit_units$uv_ri_ratio[[1]]$value,
    result_mixed_units$uv_ri_ratio[[1]]$value
  )
  expect_identical(
    attr(result_mixed_units$uv_ri_ratio, "sec_uv_ri_ratio_basis"),
    "normalized_concentration"
  )

  expect_error(
    make_result(
      make_data(rep(4, 3), rep(0.002, 3)),
      input_basis = "normalized_concentration"
    ),
    "Specify.*uv_concentration_unit.*untagged"
  )
})

test_that("UV/RI ratio rejects mixed bases, drift, and misaligned curves", {
  location <- c(1, 2, 3)
  make_data <- function(uv, ri, uv_unit = NULL, ri_unit = NULL) {
    data <- tibble::tibble(sample_id = "sample")
    data$uv <- measure::new_measure_list(
      list(measure::new_measure_tbl(location, uv))
    )
    data$ri <- measure::new_measure_list(
      list(measure::new_measure_tbl(location, ri))
    )
    if (!is.null(uv_unit)) {
      attr(data$uv, "sec_concentration_unit") <- uv_unit
    }
    if (!is.null(ri_unit)) {
      attr(data$ri, "sec_concentration_unit") <- ri_unit
    }
    data
  }
  make_recipe <- function(data) {
    recipes::recipe(~., data = data) |>
      step_sec_uv_ri_ratio(uv_col = "uv", ri_col = "ri", smooth = FALSE)
  }

  mixed_basis <- make_data(rep(1, 3), rep(1, 3), ri_unit = "g/mL")
  expect_error(
    recipes::prep(make_recipe(mixed_basis)),
    "cannot mix a raw detector"
  )

  training <- make_data(rep(4, 3), rep(0.002, 3), "mg/mL", "g/mL")
  prepped <- recipes::prep(make_recipe(training))
  changed_unit <- make_data(rep(0.004, 3), rep(0.002, 3), "g/mL", "g/mL")
  expect_error(
    recipes::bake(prepped, new_data = changed_unit),
    "conflicts with the concentration curve's recorded units"
  )

  misaligned <- make_data(rep(1, 3), rep(1, 3))
  misaligned$uv <- measure::new_measure_list(
    list(measure::new_measure_tbl(location + 0.1, rep(1, 3)))
  )
  expect_error(
    recipes::bake(recipes::prep(make_recipe(misaligned)), new_data = NULL),
    "locations must align"
  )
})
