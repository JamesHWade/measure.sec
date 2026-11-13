sec_delay_data <- function(
  location,
  target_value = seq_along(location) - 1,
  target_location = location,
  reference_location = location,
  reference_value = rep(0, length(reference_location))
) {
  data <- tibble::tibble(sample_id = "sample")
  data$reference <- measure::new_measure_list(list(
    measure::new_measure_tbl(
      location = reference_location,
      value = reference_value
    )
  ))
  data$target <- measure::new_measure_list(list(
    measure::new_measure_tbl(
      location = target_location,
      value = target_value
    )
  ))
  data
}

sec_delay_result <- function(data, args) {
  rec <- do.call(
    step_sec_detector_delay,
    c(
      list(
        recipe = recipes::recipe(data),
        reference = "reference",
        targets = "target"
      ),
      args
    )
  )
  recipes::bake(recipes::prep(rec), new_data = data)
}

test_that("positive and negative delays align impulses and ramps on both orientations", {
  for (method in c("shift", "interpolate")) {
    ascending <- sec_delay_data(
      location = 0:6,
      target_value = 0:6
    )
    positive <- sec_delay_result(
      ascending,
      list(delay_times = c(target = 1), method = method)
    )
    negative <- sec_delay_result(
      ascending,
      list(delay_times = c(target = -1), method = method)
    )
    expect_equal(positive$target[[1]]$value, c(1:6, NA_real_))
    expect_equal(negative$target[[1]]$value, c(NA_real_, 0:5))
    expect_identical(positive$target[[1]]$location, 0:6)
    expect_identical(positive$reference, ascending$reference)

    impulse <- sec_delay_data(
      location = 0:6,
      target_value = c(0, 0, 0, 5, 0, 0, 0)
    )
    positive_impulse <- sec_delay_result(
      impulse,
      list(delay_times = c(target = 1), method = method)
    )
    negative_impulse <- sec_delay_result(
      impulse,
      list(delay_times = c(target = -1), method = method)
    )
    expect_equal(
      positive_impulse$target[[1]]$value,
      c(0, 0, 5, 0, 0, 0, NA_real_)
    )
    expect_equal(
      negative_impulse$target[[1]]$value,
      c(NA_real_, 0, 0, 0, 5, 0, 0)
    )

    descending <- sec_delay_data(
      location = 6:0,
      target_value = 6:0
    )
    descending_result <- sec_delay_result(
      descending,
      list(delay_times = c(target = 1), method = method)
    )
    expect_equal(
      descending_result$target[[1]]$value,
      c(NA_real_, 6, 5, 4, 3, 2, 1)
    )
    expect_identical(descending_result$target[[1]]$location, 6:0)
  }
})

test_that("interpolation follows corrected x plus delay on irregular axes", {
  ascending <- sec_delay_data(
    location = c(0, 1, 3, 6),
    target_value = c(0, 2, 6, 12)
  )
  result <- sec_delay_result(
    ascending,
    list(
      delay_times = c(target = 0.5),
      method = "interpolate"
    )
  )
  expect_equal(result$target[[1]]$value, c(1, 3, 7, NA_real_))

  descending <- sec_delay_data(
    location = c(6, 3, 1, 0),
    target_value = c(12, 6, 2, 0)
  )
  descending_result <- sec_delay_result(
    descending,
    list(
      delay_times = c(target = 0.5),
      method = "interpolate"
    )
  )
  expect_equal(
    descending_result$target[[1]]$value,
    c(NA_real_, 7, 3, 1)
  )
  expect_error(
    sec_delay_result(
      ascending,
      list(delay_times = c(target = 0.5), method = "shift")
    ),
    "uniformly spaced"
  )
})

test_that("interpolation propagates internal missing values without bridging", {
  data <- sec_delay_data(
    location = 0:4,
    target_value = c(0, 1, NA_real_, 3, 4)
  )
  result <- sec_delay_result(
    data,
    list(
      delay_times = c(target = 0.5),
      method = "interpolate"
    )
  )
  expect_equal(
    result$target[[1]]$value,
    c(0.5, NA_real_, NA_real_, 3.5, NA_real_)
  )
})

test_that("interpolation rejects nonfinite results from finite inputs", {
  for (target_value in list(
    c(-1e308, 1e308, 0),
    c(1e308, -1e308, 0)
  )) {
    data <- sec_delay_data(
      location = 0:2,
      target_value = target_value
    )
    expect_error(
      sec_delay_result(
        data,
        list(
          delay_times = c(target = 0.5),
          method = "interpolate"
        )
      ),
      "finite or ordinary NA"
    )
  }
})

test_that("zero delay preserves irregular curves and concentration units", {
  data <- sec_delay_data(
    location = c(0, 1, 3, 6),
    target_value = c(2, NA_real_, 4, 5)
  )
  attr(data$target, "sec_concentration_unit") <- "g/mL"

  for (method in c("shift", "interpolate")) {
    result <- sec_delay_result(
      data,
      list(
        delay_volumes = c(target = 0),
        method = method,
        flow_rate = 1
      )
    )
    expect_identical(result$target[[1]], data$target[[1]])
    expect_identical(
      attr(result$target, "sec_concentration_unit"),
      "g/mL"
    )
    expect_identical(result$reference, data$reference)
  }
})

test_that("both raw delay dimensions convert across all four axis units", {
  for (flow in c(0.5, 2)) {
    for (units in c("min", "s", "mL", "uL")) {
      axis_scale <- switch(
        units,
        min = 1,
        s = 60,
        mL = flow,
        uL = flow * 1000
      )
      data <- sec_delay_data(
        location = (0:3) * axis_scale,
        target_value = 0:3
      )
      time_flow <- if (units %in% c("min", "s")) NULL else flow
      volume_flow <- if (units %in% c("min", "s")) flow else NULL

      from_time <- sec_delay_result(
        data,
        list(
          delay_times = c(target = 1),
          flow_rate = time_flow,
          elution_units = units
        )
      )
      from_volume <- sec_delay_result(
        data,
        list(
          delay_volumes = c(target = flow),
          flow_rate = volume_flow,
          elution_units = units
        )
      )
      expect_equal(from_time$target[[1]]$value, c(1, 2, 3, NA_real_))
      expect_equal(from_volume$target[[1]]$value, c(1, 2, 3, NA_real_))
    }
  }

  data <- sec_delay_data(location = 0:3, target_value = 0:3)
  huge_irrelevant_flow <- sec_delay_result(
    data,
    list(
      delay_times = c(target = 1),
      flow_rate = 1e308,
      elution_units = "min"
    )
  )
  expect_equal(
    huge_irrelevant_flow$target[[1]]$value,
    c(1, 2, 3, NA_real_)
  )
})

test_that("flow is required only for time-volume conversions", {
  time_min <- sec_delay_data(location = 0:3)
  volume_ml <- sec_delay_data(location = 0:3)

  expect_no_error(sec_delay_result(
    time_min,
    list(
      delay_times = c(target = 1),
      flow_rate = NULL,
      elution_units = "min"
    )
  ))
  expect_no_error(sec_delay_result(
    time_min,
    list(
      delay_times = c(target = 60),
      flow_rate = NULL,
      elution_units = "s"
    )
  ))
  expect_no_error(sec_delay_result(
    volume_ml,
    list(
      delay_volumes = c(target = 1),
      flow_rate = NULL,
      elution_units = "mL"
    )
  ))
  expect_no_error(sec_delay_result(
    volume_ml,
    list(
      delay_volumes = c(target = 1000),
      flow_rate = NULL,
      elution_units = "uL"
    )
  ))

  expect_error(
    step_sec_detector_delay(
      recipes::recipe(time_min),
      reference = "reference",
      targets = "target",
      delay_volumes = c(target = 1),
      flow_rate = NULL,
      elution_units = "min"
    ),
    "flow_rate"
  )
  expect_error(
    step_sec_detector_delay(
      recipes::recipe(volume_ml),
      reference = "reference",
      targets = "target",
      delay_times = c(target = 1),
      flow_rate = NULL,
      elution_units = "mL"
    ),
    "flow_rate"
  )
  expect_error(
    step_sec_detector_delay(
      recipes::recipe(time_min),
      reference = "reference",
      targets = "target",
      delay_times = c(target = 1),
      flow_rate = 0,
      elution_units = "min"
    ),
    "flow_rate"
  )
  expect_error(
    step_sec_detector_delay(
      recipes::recipe(volume_ml),
      reference = "reference",
      targets = "target",
      delay_volumes = c(target = 1),
      flow_rate = Inf,
      elution_units = "mL"
    ),
    "flow_rate"
  )
})

test_that("new steps retain mutually exclusive raw delay declarations", {
  data <- sec_delay_data(location = 0:3)
  time_recipe <- recipes::recipe(data) |>
    step_sec_detector_delay(
      reference = "reference",
      targets = "target",
      delay_times = c(target = 1),
      flow_rate = NULL
    )
  volume_recipe <- recipes::recipe(data) |>
    step_sec_detector_delay(
      reference = "reference",
      targets = "target",
      delay_volumes = c(target = 1)
    )

  time_step <- time_recipe$steps[[1]]
  volume_step <- volume_recipe$steps[[1]]
  expect_identical(time_step$delay_times, c(target = 1))
  expect_null(time_step$delay_volumes)
  expect_identical(volume_step$delay_volumes, c(target = 1))
  expect_null(volume_step$delay_times)

  prepared <- recipes::prep(time_recipe)
  details <- recipes::tidy(prepared, number = 1)
  expect_identical(details$elution_units, "min")
  expect_identical(details$delay_times[[1]], c(target = 1))
  expect_null(details$delay_volumes[[1]])

  saved <- withr::local_tempfile(fileext = ".rds")
  saveRDS(prepared, saved)
  restored <- readRDS(saved)
  expect_equal(
    recipes::bake(restored, new_data = data)$target[[1]]$value,
    c(1, 2, 3, NA_real_)
  )

  old_untrained <- volume_recipe
  old_untrained$steps[[1]]$delay_times <- NULL
  old_untrained$steps[[1]]$elution_units <- NULL
  rebuilt <- recipes::prep(old_untrained)
  expect_null(rebuilt$steps[[1]]$delay_times)
  expect_identical(rebuilt$steps[[1]]$elution_units, "min")
})

test_that("historical steps without new fields retain minute volume semantics", {
  data <- sec_delay_data(location = 0:3)
  prepared <- recipes::recipe(data) |>
    step_sec_detector_delay(
      reference = "reference",
      targets = "target",
      delay_volumes = c(target = 2),
      flow_rate = 2
    ) |>
    recipes::prep()
  historical <- prepared
  historical$steps[[1]]$delay_times <- NULL
  historical$steps[[1]]$elution_units <- NULL

  output <- recipes::bake(historical, new_data = data)
  expect_equal(output$target[[1]]$value, c(1, 2, 3, NA_real_))
  details <- recipes::tidy(historical, number = 1)
  expect_identical(details$elution_units, "min")
  expect_null(details$delay_times[[1]])

  explicit_null <- prepared
  explicit_null$steps[[1]]["elution_units"] <- list(NULL)
  expect_error(
    recipes::bake(explicit_null, new_data = data),
    "elution_units"
  )

  invalid_untrained <- recipes::recipe(data) |>
    step_sec_detector_delay(
      reference = "reference",
      targets = "target",
      delay_volumes = c(target = 1)
    )
  invalid_untrained$steps[[1]]$elution_units <- "minutes"
  expect_error(recipes::prep(invalid_untrained), "elution_units")
  expect_error(
    recipes::tidy(invalid_untrained, number = 1),
    "elution_units"
  )
})

test_that("duplicate serialized physics fields are rejected", {
  data <- sec_delay_data(location = 0:3)
  prepared <- recipes::recipe(data) |>
    step_sec_detector_delay(
      reference = "reference",
      targets = "target",
      delay_times = c(target = 1)
    ) |>
    recipes::prep()
  original_step <- prepared$steps[[1]]
  duplicate_values <- list(
    reference = "other_reference",
    targets = "other_target",
    delay_times = c(target = 2),
    delay_volumes = c(target = 2),
    flow_rate = 2,
    method = "shift"
  )

  for (field in names(duplicate_values)) {
    ambiguous_step <- c(
      unclass(original_step),
      stats::setNames(list(duplicate_values[[field]]), field)
    )
    class(ambiguous_step) <- class(original_step)
    ambiguous_recipe <- prepared
    ambiguous_recipe$steps[[1]] <- ambiguous_step
    expect_error(
      recipes::bake(ambiguous_recipe, new_data = data),
      "duplicate configuration fields"
    )
  }

  unrelated <- c(
    unclass(original_step),
    stats::setNames(
      list("first", "second"),
      c("vendor_metadata", "vendor_metadata")
    )
  )
  class(unrelated) <- class(original_step)
  expect_no_error(.sec_delay_resolve_fields(unrelated))
})

test_that("empty targets and unused valid delays are allowed", {
  data <- sec_delay_data(location = 0:3)
  rec <- recipes::recipe(data) |>
    step_sec_detector_delay(
      reference = "reference",
      targets = character(),
      delay_times = c(unused = 100),
      flow_rate = NULL
    ) |>
    recipes::prep()
  output <- recipes::bake(rec, new_data = data)
  expect_identical(output$target, data$target)
  expect_identical(output$reference, data$reference)
})

test_that("current bake rows validate selected columns and exact grids", {
  training <- sec_delay_data(location = 0:3)
  prepared <- recipes::recipe(training) |>
    step_sec_detector_delay(
      reference = "reference",
      targets = "target",
      delay_times = c(target = 1)
    ) |>
    recipes::prep()

  mismatch <- sec_delay_data(
    location = 0:3,
    target_location = c(0, 1, 2, 4)
  )
  expect_snapshot(
    error = TRUE,
    recipes::bake(prepared, new_data = mismatch)
  )

  missing_target <- training
  missing_target$target <- NULL
  expect_error(
    recipes::bake(prepared, new_data = missing_target),
    "target|measure column|variables"
  )

  missing_reference <- training
  missing_reference$reference <- NULL
  expect_error(
    recipes::bake(prepared, new_data = missing_reference),
    "reference|measure column|variables"
  )

  missing_delay <- prepared
  missing_delay$steps[[1]]$delay_times <- c(other = 1)
  expect_error(
    recipes::bake(missing_delay, new_data = training),
    "Delays not specified"
  )

  integer_reference <- sec_delay_data(
    location = 0:3,
    target_location = as.double(0:3)
  )
  expect_no_error(recipes::bake(prepared, new_data = integer_reference))
})

test_that("each reordered current row uses its own reference grid", {
  training <- tibble::tibble(sample_id = c("train-a", "train-b"))
  training$reference <- measure::new_measure_list(list(
    measure::new_measure_tbl(0:2, c(0, 0, 0)),
    measure::new_measure_tbl(0:2, c(0, 0, 0))
  ))
  training$target <- measure::new_measure_list(list(
    measure::new_measure_tbl(0:2, c(0, 1, 2)),
    measure::new_measure_tbl(0:2, c(0, 10, 20))
  ))
  prepared <- recipes::recipe(training) |>
    step_sec_detector_delay(
      reference = "reference",
      targets = "target",
      delay_times = c(target = 1)
    ) |>
    recipes::prep()

  current <- tibble::tibble(sample_id = c("b", "a"))
  current$reference <- measure::new_measure_list(list(
    measure::new_measure_tbl(10:12, c(1, 1, 1)),
    measure::new_measure_tbl(0:2, c(2, 2, 2))
  ))
  current$target <- measure::new_measure_list(list(
    measure::new_measure_tbl(10:12, c(100, 110, 120)),
    measure::new_measure_tbl(0:2, c(1, 2, 3))
  ))

  output <- recipes::bake(prepared, new_data = current)
  expect_identical(output$sample_id, c("b", "a"))
  expect_equal(output$target[[1]]$value, c(110, 120, NA_real_))
  expect_equal(output$target[[2]]$value, c(2, 3, NA_real_))
  expect_identical(output$reference, current$reference)
})

test_that("invalid declarations and malformed curves fail clearly", {
  data <- sec_delay_data(location = 0:3)
  recipe <- recipes::recipe(data)

  for (units in list(
    NULL,
    NA_character_,
    c("min", "s"),
    factor("min"),
    "minutes"
  )) {
    expect_error(
      step_sec_detector_delay(
        recipe,
        reference = "reference",
        targets = "target",
        delay_times = c(target = 1),
        elution_units = units
      ),
      "elution_units"
    )
  }
  expect_error(
    step_sec_detector_delay(
      recipe,
      reference = "reference",
      targets = "target",
      delay_volumes = numeric()
    ),
    "delay_volumes"
  )
  expect_error(
    step_sec_detector_delay(
      recipe,
      reference = "reference",
      targets = "target",
      delay_volumes = c(target = 1, target = 2)
    ),
    "unique"
  )
  expect_error(
    step_sec_detector_delay(
      recipe,
      reference = "reference",
      targets = "target",
      delay_volumes = c(" " = 1)
    ),
    "nonblank"
  )
  expect_error(
    step_sec_detector_delay(
      recipe,
      reference = "reference",
      targets = "reference",
      delay_volumes = c(reference = 1)
    ),
    "cannot also be a target"
  )
  expect_error(
    step_sec_detector_delay(
      recipe,
      reference = "reference",
      targets = "target",
      delay_volumes = c(target = 1),
      delay_times = c(target = 1)
    ),
    "exactly one"
  )

  for (bad in list(
    sec_delay_data(c(0, 1, 1, 2)),
    sec_delay_data(c(0, 2, 1, 3)),
    sec_delay_data(c(0, 1, Inf, 3)),
    sec_delay_data(c(-1e308, 0, 1e308), target_value = c(0, 1, 2))
  )) {
    expect_error(
      sec_delay_result(
        bad,
        list(delay_times = c(target = 1), method = "interpolate")
      ),
      "finite|monotonic|representable"
    )
  }

  bad_values <- sec_delay_data(0:3, target_value = c(0, NaN, 2, 3))
  expect_error(
    sec_delay_result(
      bad_values,
      list(delay_times = c(target = 1), method = "interpolate")
    ),
    "ordinary NA"
  )
  bad_infinite_values <- sec_delay_data(0:3, target_value = c(0, Inf, 2, 3))
  expect_error(
    sec_delay_result(
      bad_infinite_values,
      list(delay_times = c(target = 1), method = "interpolate")
    ),
    "ordinary NA"
  )

  malformed_curve <- structure(
    list(location = 0:1, value = 0:2),
    class = "data.frame",
    row.names = 1:3
  )
  expect_error(.sec_delay_curve(malformed_curve), "equal-length numeric")
  expect_error(
    .sec_delay_curve(data.frame(location = c("0", "1"), value = c(0, 1))),
    "equal-length numeric"
  )

  bad_reference <- sec_delay_data(
    location = 0:3,
    reference_location = c(0, 1, 1, 3)
  )
  expect_error(
    sec_delay_result(
      bad_reference,
      list(delay_times = c(target = 0), method = "shift")
    ),
    "monotonic"
  )
})

test_that("zero and singleton curves only accept zero correction", {
  for (location in list(numeric(), 2)) {
    data <- sec_delay_data(
      location = location,
      target_value = numeric(length(location))
    )
    zero <- sec_delay_result(
      data,
      list(delay_times = c(target = 0), method = "interpolate")
    )
    expect_identical(zero$target[[1]], data$target[[1]])

    expect_error(
      sec_delay_result(
        data,
        list(delay_times = c(target = 1), method = "interpolate")
      ),
      "at least two"
    )
  }
})

test_that("shift bounds huge offsets and preserves round-to-even quantization", {
  data <- sec_delay_data(location = 0:3, target_value = 0:3)
  huge <- sec_delay_result(
    data,
    list(delay_times = c(target = 1e308), method = "shift")
  )
  expect_true(all(is.na(huge$target[[1]]$value)))

  quantized_zero <- sec_delay_result(
    data,
    list(delay_times = c(target = 0.5), method = "shift")
  )
  expect_identical(quantized_zero$target[[1]], data$target[[1]])
  quantized_two <- sec_delay_result(
    data,
    list(delay_times = c(target = 1.5), method = "shift")
  )
  expect_equal(quantized_two$target[[1]]$value, c(2, 3, NA_real_, NA_real_))
})

test_that("conversion underflow, overflow, and interpolation cancellation are rejected", {
  data <- sec_delay_data(location = 0:3)
  expect_error(
    step_sec_detector_delay(
      recipes::recipe(data),
      reference = "reference",
      targets = "target",
      delay_volumes = c(target = 1e308),
      flow_rate = 1e-308
    ),
    "nonfinite delay"
  )
  expect_error(
    step_sec_detector_delay(
      recipes::recipe(data),
      reference = "reference",
      targets = "target",
      delay_volumes = c(target = 5e-324),
      flow_rate = 1e308
    ),
    "underflowed to zero"
  )

  large_axis <- sec_delay_data(
    location = c(1e16, 1e16 + 2, 1e16 + 4),
    target_value = c(1, 2, 3)
  )
  expect_error(
    sec_delay_result(
      large_axis,
      list(
        delay_times = c(target = 1),
        method = "interpolate"
      )
    ),
    "not representable"
  )
})

test_that("corrected concentration tags and grids remain usable by LALS", {
  location <- 0:4
  dn_dc <- 0.185
  concentration_g_per_ml <- c(0.001, 0.002, 0.003, 0.004, 0.005)
  target_concentration <- c(0.002, 0.003, 0.004, 0.005, 0.006)
  reference_concentration <- c(0.003, 0.004, 0.005, 0.006, 0.006)
  wavelength_nm <- 670
  lambda_cm <- wavelength_nm * 1e-7
  optical_constant <- 4 *
    pi^2 *
    1.333^2 *
    dn_dc^2 /
    (6.022e23 * lambda_cm^4)
  expected_mw <- 50

  data <- tibble::tibble(sample_id = "sample")
  data$ri <- measure::new_measure_list(list(
    measure::new_measure_tbl(
      location = location,
      value = target_concentration * dn_dc
    )
  ))
  data$lals <- measure::new_measure_list(list(
    measure::new_measure_tbl(
      location = location,
      value = optical_constant * reference_concentration * expected_mw
    )
  ))

  rec <- recipes::recipe(data) |>
    step_sec_ri(dn_dc = dn_dc) |>
    step_sec_detector_delay(
      reference = "lals",
      targets = "ri",
      delay_times = c(ri = 1),
      method = "interpolate"
    ) |>
    step_sec_lals(
      measures = "lals",
      concentration_col = "ri",
      dn_dc = dn_dc,
      calibration_constant = 1
    )
  output <- recipes::bake(recipes::prep(rec), new_data = NULL)

  expect_identical(attr(output$ri, "sec_concentration_unit"), "g/mL")
  expect_identical(output$ri[[1]]$location, location)
  expect_equal(
    output$mw_lals[[1]]$value,
    c(rep(expected_mw, 4), NA_real_),
    tolerance = 1e-8
  )
  expect_identical(output$mw_lals[[1]]$location, location)
})
