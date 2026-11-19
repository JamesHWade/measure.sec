# Tests for the SEC bridge to measure value metadata.

test_that("calibrated SEC concentration writers declare their value quantity", {
  data <- tibble::tibble(sample_id = "sample")
  data$ri <- measure::new_measure_list(list(
    measure::new_measure_tbl(1:3, c(1, 2, 1))
  ))
  data$uv <- measure::new_measure_list(list(
    measure::new_measure_tbl(1:3, c(1, 2, 1))
  ))

  ri_result <- recipes::bake(
    recipes::prep(
      recipes::recipe(data) |>
        step_sec_ri(measures = "ri", dn_dc = 0.2)
    ),
    new_data = NULL
  )
  uv_result <- recipes::bake(
    recipes::prep(
      recipes::recipe(data) |>
        step_sec_uv(measures = "uv", extinction_coef = 2)
    ),
    new_data = NULL
  )
  concentration_result <- recipes::bake(
    recipes::prep(
      recipes::recipe(data) |>
        step_sec_concentration(
          measures = "ri",
          injection_mass = 3,
          flow_rate = 1,
          concentration_units = "mg/mL"
        )
    ),
    new_data = NULL
  )

  expect_identical(
    measure::measure_value_metadata(ri_result$ri),
    list(quantity = "mass_concentration", unit = "g/mL")
  )
  expect_identical(attr(ri_result$ri, "sec_concentration_unit"), "g/mL")
  expect_identical(
    measure::measure_value_metadata(uv_result$uv),
    list(quantity = "mass_concentration", unit = "mg/mL")
  )
  expect_identical(attr(uv_result$uv, "sec_concentration_unit"), "mg/mL")
  expect_identical(
    measure::measure_value_metadata(concentration_result$ri),
    list(quantity = "mass_concentration", unit = "mg/mL")
  )
  expect_identical(
    attr(concentration_result$ri, "sec_concentration_unit"),
    "mg/mL"
  )
})

test_that("baseline and detector delay preserve value metadata", {
  location <- seq_len(20)
  data <- tibble::tibble(sample_id = "sample")
  data$ri <- measure::new_measure_list(list(
    measure::new_measure_tbl(location, 10 + 0.1 * location)
  ))
  data$ri <- measure::set_measure_value_metadata(
    data$ri,
    quantity = "mass_concentration",
    unit = "mg/mL"
  )
  attr(data$ri, "sec_concentration_unit") <- "mg/mL"

  baseline_result <- recipes::bake(
    recipes::prep(
      recipes::recipe(data) |>
        step_sec_baseline(
          measures = "ri",
          left_frac = 0.2,
          right_frac = 0.2
        )
    ),
    new_data = NULL
  )
  expect_identical(
    measure::measure_value_metadata(baseline_result$ri),
    list(quantity = "mass_concentration", unit = "mg/mL")
  )
  expect_identical(
    attr(baseline_result$ri, "sec_concentration_unit"),
    "mg/mL"
  )

  delay_data <- tibble::tibble(sample_id = "sample")
  delay_data$reference <- measure::new_measure_list(list(
    measure::new_measure_tbl(1:3, c(0, 0, 0))
  ))
  delay_data$target <- measure::new_measure_list(list(
    measure::new_measure_tbl(1:3, c(1, 2, 3))
  ))
  delay_data$target <- measure::set_measure_value_metadata(
    delay_data$target,
    quantity = "mass_concentration",
    unit = "g/mL"
  )
  attr(delay_data$target, "sec_concentration_unit") <- "g/mL"

  delay_result <- recipes::bake(
    recipes::prep(
      recipes::recipe(delay_data) |>
        step_sec_detector_delay(
          reference = "reference",
          targets = "target",
          delay_times = c(target = 0),
          method = "interpolate"
        )
    ),
    new_data = NULL
  )
  expect_identical(
    measure::measure_value_metadata(delay_result$target),
    list(quantity = "mass_concentration", unit = "g/mL")
  )
  expect_identical(
    attr(delay_result$target, "sec_concentration_unit"),
    "g/mL"
  )
})

test_that("legacy SEC concentration attributes are upgraded by preserving steps", {
  location <- seq_len(20)
  baseline_data <- tibble::tibble(sample_id = "sample")
  baseline_data$ri <- measure::new_measure_list(list(
    measure::new_measure_tbl(location, 10 + 0.1 * location)
  ))
  attr(baseline_data$ri, "sec_concentration_unit") <- "mg/mL"

  baseline_result <- recipes::bake(
    recipes::prep(
      recipes::recipe(baseline_data) |>
        step_sec_baseline(
          measures = "ri",
          left_frac = 0.2,
          right_frac = 0.2
        )
    ),
    new_data = NULL
  )
  expect_identical(
    measure::measure_value_metadata(baseline_result$ri),
    list(quantity = "mass_concentration", unit = "mg/mL")
  )

  delay_data <- tibble::tibble(sample_id = "sample")
  delay_data$reference <- measure::new_measure_list(list(
    measure::new_measure_tbl(1:3, c(0, 0, 0))
  ))
  delay_data$target <- measure::new_measure_list(list(
    measure::new_measure_tbl(1:3, c(1, 2, 3))
  ))
  attr(delay_data$target, "sec_concentration_unit") <- "g/mL"

  delay_result <- recipes::bake(
    recipes::prep(
      recipes::recipe(delay_data) |>
        step_sec_detector_delay(
          reference = "reference",
          targets = "target",
          delay_times = c(target = 0),
          method = "interpolate"
        )
    ),
    new_data = NULL
  )
  expect_identical(
    measure::measure_value_metadata(delay_result$target),
    list(quantity = "mass_concentration", unit = "g/mL")
  )
})

test_that("new metadata drives normalized UV/RI ratio consumers", {
  location <- 1:3
  data <- tibble::tibble(sample_id = "sample")
  data$ri <- measure::new_measure_list(list(
    measure::new_measure_tbl(location, c(1, 2, 1))
  ))
  data$uv <- measure::new_measure_list(list(
    measure::new_measure_tbl(location, c(0.5, 1, 0.5))
  ))
  data$ri <- measure::set_measure_value_metadata(
    data$ri,
    quantity = "mass_concentration",
    unit = "mg/mL"
  )
  data$uv <- measure::set_measure_value_metadata(
    data$uv,
    quantity = "mass_concentration",
    unit = "mg/mL"
  )

  result <- recipes::bake(
    recipes::prep(
      recipes::recipe(data) |>
        step_sec_uv_ri_ratio(
          uv_col = "uv",
          ri_col = "ri",
          input_basis = "normalized_concentration",
          smooth = FALSE
        )
    ),
    new_data = NULL
  )

  expect_equal(result$uv_ri_ratio[[1]]$value, rep(0.5, 3))
})

test_that("raw SEC consumers accept native detector-response metadata", {
  data <- tibble::tibble(sample_id = "sample")
  data$ri <- measure::new_measure_list(list(
    measure::new_measure_tbl(1:3, c(1, 2, 1))
  ))
  data$uv <- measure::new_measure_list(list(
    measure::new_measure_tbl(1:3, c(0.5, 1, 0.5))
  ))
  data$ri <- measure::set_measure_value_metadata(
    data$ri,
    quantity = "detector_response",
    unit = "RI-AU"
  )
  data$uv <- measure::set_measure_value_metadata(
    data$uv,
    quantity = "detector_response",
    unit = "UV-AU"
  )

  ratio <- recipes::bake(
    recipes::prep(
      recipes::recipe(data) |>
        step_sec_uv_ri_ratio(
          uv_col = "uv",
          ri_col = "ri",
          input_basis = "detector_response",
          smooth = FALSE
        )
    ),
    new_data = NULL
  )
  auto_ratio <- recipes::bake(
    recipes::prep(
      recipes::recipe(data) |>
        step_sec_uv_ri_ratio(
          uv_col = "uv",
          ri_col = "ri",
          input_basis = "auto",
          smooth = FALSE
        )
    ),
    new_data = NULL
  )
  composition <- recipes::bake(
    recipes::prep(
      recipes::recipe(data) |>
        step_sec_composition(
          uv_col = "uv",
          ri_col = "ri",
          component_a_uv = 1,
          component_a_ri = 0.2,
          component_b_uv = 0.1,
          component_b_ri = 0.8,
          input_basis = "detector_response"
        )
    ),
    new_data = NULL
  )

  expect_equal(ratio$uv_ri_ratio[[1]]$value, rep(0.5, 3))
  expect_equal(auto_ratio$uv_ri_ratio[[1]]$value, rep(0.5, 3))
  expect_true(is.list(composition$composition_a[[1]]))
})

test_that("raw SEC fits enforce declared native units per detector", {
  make_data <- function(
    uv_unit = "UV-AU",
    ri_unit = "RI-AU",
    uv_scale = 1,
    ri_scale = 1
  ) {
    data <- tibble::tibble(sample_id = "sample")
    data$ri <- measure::new_measure_list(list(
      measure::new_measure_tbl(1:3, ri_scale * c(1, 2, 1))
    ))
    data$uv <- measure::new_measure_list(list(
      measure::new_measure_tbl(1:3, uv_scale * c(0.5, 1, 0.5))
    ))
    if (!is.null(ri_unit)) {
      data$ri <- measure::set_measure_value_metadata(
        data$ri,
        quantity = "detector_response",
        unit = ri_unit
      )
    }
    if (!is.null(uv_unit)) {
      data$uv <- measure::set_measure_value_metadata(
        data$uv,
        quantity = "detector_response",
        unit = uv_unit
      )
    }
    data
  }

  training <- make_data()
  ratio_fit <- recipes::prep(
    recipes::recipe(training) |>
      step_sec_uv_ri_ratio(
        uv_col = "uv",
        ri_col = "ri",
        input_basis = "detector_response",
        smooth = FALSE
      )
  )
  composition_fit <- recipes::prep(
    recipes::recipe(training) |>
      step_sec_composition(
        uv_col = "uv",
        ri_col = "ri",
        component_a_uv = 1,
        component_a_ri = 0.2,
        component_b_uv = 0.1,
        component_b_ri = 0.8,
        input_basis = "detector_response"
      )
  )

  expect_identical(ratio_fit$steps[[1]]$uv_raw_unit, "UV-AU")
  expect_identical(ratio_fit$steps[[1]]$ri_raw_unit, "RI-AU")
  expect_identical(composition_fit$steps[[1]]$uv_raw_unit, "UV-AU")
  expect_identical(composition_fit$steps[[1]]$ri_raw_unit, "RI-AU")
  expect_silent(recipes::bake(ratio_fit, new_data = training))
  expect_silent(recipes::bake(composition_fit, new_data = training))

  changed_uv <- make_data(uv_unit = "UV-mAU", uv_scale = 1000)
  expect_error(
    recipes::bake(ratio_fit, new_data = changed_uv),
    "unit"
  )
  expect_error(
    recipes::bake(composition_fit, new_data = changed_uv),
    "unit"
  )

  lost_ri <- make_data(ri_unit = NULL)
  expect_error(
    recipes::bake(ratio_fit, new_data = lost_ri),
    "retain raw detector-response unit"
  )
  expect_error(
    recipes::bake(composition_fit, new_data = lost_ri),
    "retain raw detector-response unit"
  )

  persisted_ratio <- tempfile(fileext = ".rds")
  persisted_composition <- tempfile(fileext = ".rds")
  legacy_ratio <- tempfile(fileext = ".rds")
  legacy_composition <- tempfile(fileext = ".rds")
  on.exit(
    unlink(
      c(
        persisted_ratio,
        persisted_composition,
        legacy_ratio,
        legacy_composition
      )
    ),
    add = TRUE
  )
  saveRDS(ratio_fit, persisted_ratio)
  saveRDS(composition_fit, persisted_composition)
  expect_silent(
    recipes::bake(readRDS(persisted_ratio), new_data = training)
  )
  expect_silent(
    recipes::bake(readRDS(persisted_composition), new_data = training)
  )

  old_ratio_fit <- ratio_fit
  old_ratio_fit$steps[[1]]$uv_raw_unit <- NULL
  old_ratio_fit$steps[[1]]$ri_raw_unit <- NULL
  old_composition_fit <- composition_fit
  old_composition_fit$steps[[1]]$uv_raw_unit <- NULL
  old_composition_fit$steps[[1]]$ri_raw_unit <- NULL
  saveRDS(old_ratio_fit, legacy_ratio)
  saveRDS(old_composition_fit, legacy_composition)
  expect_silent(
    recipes::bake(readRDS(legacy_ratio), new_data = changed_uv)
  )
  expect_silent(
    recipes::bake(readRDS(legacy_composition), new_data = changed_uv)
  )

  untagged <- make_data(uv_unit = NULL, ri_unit = NULL)
  newly_declared <- make_data(uv_unit = "UV-mAU", ri_unit = "RI-mAU")
  untagged_ratio <- recipes::prep(
    recipes::recipe(untagged) |>
      step_sec_uv_ri_ratio(
        uv_col = "uv",
        ri_col = "ri",
        input_basis = "detector_response",
        smooth = FALSE
      )
  )
  untagged_composition <- recipes::prep(
    recipes::recipe(untagged) |>
      step_sec_composition(
        uv_col = "uv",
        ri_col = "ri",
        component_a_uv = 1,
        component_a_ri = 0.2,
        component_b_uv = 0.1,
        component_b_ri = 0.8,
        input_basis = "detector_response"
      )
  )
  expect_silent(
    recipes::bake(untagged_ratio, new_data = newly_declared)
  )
  expect_silent(
    recipes::bake(untagged_composition, new_data = newly_declared)
  )
})

test_that("raw SEC readers reject unsupported quantities before unit checks", {
  data <- tibble::tibble(sample_id = "sample")
  data$ri <- measure::new_measure_list(list(
    measure::new_measure_tbl(1:3, c(1, 2, 1))
  ))
  data$uv <- measure::new_measure_list(list(
    measure::new_measure_tbl(1:3, c(0.5, 1, 0.5))
  ))
  data$ri <- measure::set_measure_value_metadata(
    data$ri,
    quantity = "detector_response",
    unit = "RI-AU"
  )
  data$uv <- measure::set_measure_value_metadata(
    data$uv,
    quantity = "detector_response",
    unit = "UV-AU"
  )
  ratio_fit <- recipes::prep(
    recipes::recipe(data) |>
      step_sec_uv_ri_ratio(
        uv_col = "uv",
        ri_col = "ri",
        input_basis = "detector_response",
        smooth = FALSE
      )
  )
  composition_fit <- recipes::prep(
    recipes::recipe(data) |>
      step_sec_composition(
        uv_col = "uv",
        ri_col = "ri",
        component_a_uv = 1,
        component_a_ri = 0.2,
        component_b_uv = 0.1,
        component_b_ri = 0.8,
        input_basis = "detector_response"
      )
  )
  invalid <- data
  invalid$uv <- measure::set_measure_value_metadata(
    invalid$uv,
    quantity = "normalized_signal",
    unit = "UV-AU"
  )

  expect_error(
    recipes::bake(ratio_fit, new_data = invalid),
    "unsupported value quantity"
  )
  expect_error(
    recipes::bake(composition_fit, new_data = invalid),
    "unsupported value quantity"
  )
})

test_that("auto raw mode accepts mixed tagged and untagged response curves", {
  data <- tibble::tibble(sample_id = "sample")
  data$ri <- measure::new_measure_list(list(
    measure::new_measure_tbl(1:3, c(1, 2, 1))
  ))
  data$uv <- measure::new_measure_list(list(
    measure::new_measure_tbl(1:3, c(0.5, 1, 0.5))
  ))
  data$ri <- measure::set_measure_value_metadata(
    data$ri,
    quantity = "detector_response",
    unit = "RI-AU"
  )

  result <- recipes::bake(
    recipes::prep(
      recipes::recipe(data) |>
        step_sec_uv_ri_ratio(
          uv_col = "uv",
          ri_col = "ri",
          input_basis = "auto",
          smooth = FALSE
        )
    ),
    new_data = NULL
  )

  expect_equal(result$uv_ri_ratio[[1]]$value, rep(0.5, 3))

  mixed_data <- data
  mixed_data$uv <- measure::set_measure_value_metadata(
    mixed_data$uv,
    quantity = "mass_concentration",
    unit = "mg/mL"
  )
  expect_error(
    recipes::prep(
      recipes::recipe(mixed_data) |>
        step_sec_uv_ri_ratio(
          uv_col = "uv",
          ri_col = "ri",
          input_basis = "auto",
          smooth = FALSE
        )
    ),
    "cannot mix"
  )
})

test_that("SEC concentration consumers reject wrong or conflicting declarations", {
  data <- tibble::tibble(sample_id = "sample")
  data$ri <- measure::new_measure_list(list(
    measure::new_measure_tbl(1:3, c(1, 2, 1))
  ))
  data$mals <- measure::new_measure_list(list(
    measure::new_measure_tbl(1:3, c(1, 2, 1))
  ))

  wrong_data <- data
  wrong_data$ri <- measure::set_measure_value_metadata(
    wrong_data$ri,
    quantity = "normalized_signal",
    unit = "mg/mL"
  )
  wrong_recipe <- recipes::recipe(wrong_data) |>
    step_sec_lals(
      measures = "mals",
      concentration_col = "ri",
      concentration_unit = "mg/mL",
      dn_dc = 0.185,
      calibration_constant = 1
    )
  expect_error(
    recipes::prep(wrong_recipe),
    "quantity|mass_concentration"
  )

  conflicting_data <- data
  conflicting_data$ri <- measure::set_measure_value_metadata(
    conflicting_data$ri,
    quantity = "mass_concentration",
    unit = "mg/mL"
  )
  attr(conflicting_data$ri, "sec_concentration_unit") <- "g/mL"
  conflicting_recipe <- recipes::recipe(conflicting_data) |>
    step_sec_lals(
      measures = "mals",
      concentration_col = "ri",
      concentration_unit = "mg/mL",
      dn_dc = 0.185,
      calibration_constant = 1
    )
  expect_error(recipes::prep(conflicting_recipe), "conflicts")
})

test_that("raw ratio and composition reject declared concentration inputs", {
  data <- tibble::tibble(sample_id = "sample")
  data$ri <- measure::new_measure_list(list(
    measure::new_measure_tbl(1:3, c(1, 2, 1))
  ))
  data$uv <- measure::new_measure_list(list(
    measure::new_measure_tbl(1:3, c(0.5, 1, 0.5))
  ))
  data$ri <- measure::set_measure_value_metadata(
    data$ri,
    quantity = "mass_concentration",
    unit = "mg/mL"
  )
  data$uv <- measure::set_measure_value_metadata(
    data$uv,
    quantity = "mass_concentration",
    unit = "mg/mL"
  )

  ratio_recipe <- recipes::recipe(data) |>
    step_sec_uv_ri_ratio(
      uv_col = "uv",
      ri_col = "ri",
      input_basis = "detector_response",
      smooth = FALSE
    )
  expect_error(recipes::prep(ratio_recipe), "Raw detector-response")

  composition_recipe <- recipes::recipe(data) |>
    step_sec_composition(
      uv_col = "uv",
      ri_col = "ri",
      component_a_uv = 1,
      component_a_ri = 0.2,
      component_b_uv = 0.1,
      component_b_ri = 0.8,
      input_basis = "detector_response"
    )
  expect_error(recipes::prep(composition_recipe), "Raw detector-response")

  raw_data <- tibble::tibble(sample_id = "sample")
  raw_data$ri <- measure::new_measure_list(list(
    measure::new_measure_tbl(1:3, c(1, 2, 1))
  ))
  raw_data$uv <- measure::new_measure_list(list(
    measure::new_measure_tbl(1:3, c(0.5, 1, 0.5))
  ))
  trained_ratio <- recipes::prep(
    recipes::recipe(raw_data) |>
      step_sec_uv_ri_ratio(
        uv_col = "uv",
        ri_col = "ri",
        input_basis = "detector_response",
        smooth = FALSE
      )
  )
  trained_composition <- recipes::prep(
    recipes::recipe(raw_data) |>
      step_sec_composition(
        uv_col = "uv",
        ri_col = "ri",
        component_a_uv = 1,
        component_a_ri = 0.2,
        component_b_uv = 0.1,
        component_b_ri = 0.8,
        input_basis = "detector_response"
      )
  )
  expect_error(
    recipes::bake(trained_ratio, new_data = data),
    "cannot bake"
  )
  expect_error(
    recipes::bake(trained_composition, new_data = data),
    "cannot bake"
  )
})

test_that("ratio and composition reject unsupported value quantities at prep", {
  data <- tibble::tibble(sample_id = "sample")
  data$ri <- measure::new_measure_list(list(
    measure::new_measure_tbl(1:3, c(1, 2, 1))
  ))
  data$uv <- measure::new_measure_list(list(
    measure::new_measure_tbl(1:3, c(0.5, 1, 0.5))
  ))
  data$ri <- measure::set_measure_value_metadata(
    data$ri,
    quantity = "normalized_signal",
    unit = "RI-AU"
  )

  ratio_recipe <- recipes::recipe(data) |>
    step_sec_uv_ri_ratio(
      uv_col = "uv",
      ri_col = "ri",
      input_basis = "auto",
      smooth = FALSE
    )
  expect_error(recipes::prep(ratio_recipe), "unsupported value quantity")

  composition_recipe <- recipes::recipe(data) |>
    step_sec_composition(
      uv_col = "uv",
      ri_col = "ri",
      component_a_uv = 1,
      component_a_ri = 0.2,
      component_b_uv = 0.1,
      component_b_ri = 0.8,
      input_basis = "detector_response"
    )
  expect_error(
    recipes::prep(composition_recipe),
    "unsupported value quantity"
  )
})

test_that("ratio and composition reject unsupported value quantities at bake", {
  data <- tibble::tibble(sample_id = "sample")
  data$ri <- measure::new_measure_list(list(
    measure::new_measure_tbl(1:3, c(1, 2, 1))
  ))
  data$uv <- measure::new_measure_list(list(
    measure::new_measure_tbl(1:3, c(0.5, 1, 0.5))
  ))

  ratio_recipe <- recipes::prep(
    recipes::recipe(data) |>
      step_sec_uv_ri_ratio(
        uv_col = "uv",
        ri_col = "ri",
        input_basis = "detector_response",
        smooth = FALSE
      )
  )
  composition_recipe <- recipes::prep(
    recipes::recipe(data) |>
      step_sec_composition(
        uv_col = "uv",
        ri_col = "ri",
        component_a_uv = 1,
        component_a_ri = 0.2,
        component_b_uv = 0.1,
        component_b_ri = 0.8,
        input_basis = "detector_response"
      )
  )
  wrong_data <- data
  wrong_data$ri <- measure::set_measure_value_metadata(
    wrong_data$ri,
    quantity = "signal_derivative",
    unit = "RI-AU"
  )

  expect_error(
    recipes::bake(ratio_recipe, new_data = wrong_data),
    "unsupported value quantity"
  )
  expect_error(
    recipes::bake(composition_recipe, new_data = wrong_data),
    "unsupported value quantity"
  )
})
