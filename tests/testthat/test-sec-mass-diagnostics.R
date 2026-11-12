sec_diagnostics_curve <- function(location, value) {
  measure::new_measure_tbl(location = location, value = value)
}

sec_diagnostics_data <- function(...) {
  curves <- list(...)
  n <- length(curves[[1]])
  data <- tibble::tibble(sample_id = seq_len(n))

  for (column in names(curves)) {
    data[[column]] <- measure::new_measure_list(curves[[column]])
  }

  data
}

sec_diagnostics_recipe <- function(
  training,
  measures = names(training)[-1],
  injection_mass = 6,
  flow_rate = 1,
  concentration_units = "mg/mL",
  normalize_to_mass = TRUE
) {
  recipes::recipe(training) |>
    step_sec_concentration(
      measures = measures,
      injection_mass = injection_mass,
      flow_rate = flow_rate,
      concentration_units = concentration_units,
      normalize_to_mass = normalize_to_mass
    ) |>
    recipes::prep()
}

sec_capture_clipping_warnings <- function(expr) {
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

test_that("mass normalization reports sampled positive and negative areas", {
  axis <- c(0, 1, 3)
  curves <- list(
    sec_diagnostics_curve(axis, c(-2, 2, 0)),
    sec_diagnostics_curve(axis, c(-20, 2, 0)),
    sec_diagnostics_curve(axis, c(2, 2, 2))
  )
  training <- sec_diagnostics_data(
    ri = list(sec_diagnostics_curve(axis, c(1, 2, 1)))
  )
  prepped <- sec_diagnostics_recipe(training, measures = "ri")
  input <- sec_diagnostics_data(ri = curves)
  input_before <- input

  captured <- sec_capture_clipping_warnings(
    recipes::bake(prepped, new_data = input)
  )

  expect_length(captured$warnings, 1)
  warning <- captured$warnings[[1]]
  expect_s3_class(warning, "measure_sec_negative_signal_clipped")
  expect_identical(warning$column, "ri")
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
  expect_s3_class(warning$diagnostics, "tbl_df")
  expect_identical(warning$diagnostics$row, 1:3)
  expect_identical(warning$diagnostics$negative_points, c(1L, 1L, 0L))
  expect_identical(warning$diagnostics$points, c(3L, 3L, 3L))
  expect_identical(warning$diagnostics$positive_area, c(3, 3, 6))
  expect_identical(warning$diagnostics$negative_area, c(1, 10, 0))
  expect_equal(
    warning$diagnostics$excluded_fraction,
    c(1 / 4, 10 / 13, 0),
    tolerance = 2e-16
  )

  expect_snapshot(cat(conditionMessage(warning), "\n", sep = ""))
  expect_identical(
    lapply(captured$value$ri, function(m) m$value),
    list(c(0, 4, 0), c(0, 4, 0), c(2, 2, 2))
  )
  expect_identical(names(captured$value), names(input))
  expect_identical(class(captured$value$ri), class(input$ri))
  output_attributes <- attributes(captured$value$ri)
  input_attributes <- attributes(input$ri)
  output_attributes$sec_concentration_unit <- NULL
  expect_identical(output_attributes, input_attributes)
  expect_identical(attr(captured$value$ri, "sec_concentration_unit"), "mg/mL")
  expect_identical(input, input_before)

  prepped_g <- sec_diagnostics_recipe(
    training,
    measures = "ri",
    concentration_units = "g/mL"
  )
  result_g <- sec_capture_clipping_warnings(
    recipes::bake(prepped_g, new_data = input)
  )
  expect_identical(
    lapply(result_g$value$ri, function(m) m$value),
    list(c(0, 0.004, 0), c(0, 0.004, 0), c(0.002, 0.002, 0.002))
  )
})

test_that("diagnostics follow current row order, reversal, and signal scale", {
  axis <- c(0, 1, 3)
  curves <- list(
    sec_diagnostics_curve(axis, c(-2, 2, 0)),
    sec_diagnostics_curve(axis, c(-20, 2, 0)),
    sec_diagnostics_curve(axis, c(2, 2, 2))
  )
  training <- sec_diagnostics_data(
    ri = list(sec_diagnostics_curve(axis, c(1, 2, 1)))
  )
  prepped <- sec_diagnostics_recipe(training, measures = "ri")
  transformed <- lapply(
    curves[c(3, 2, 1)],
    function(m) {
      sec_diagnostics_curve(
        location = rev(m$location),
        value = rev(m$value) * 10
      )
    }
  )
  input <- sec_diagnostics_data(ri = transformed)

  captured <- sec_capture_clipping_warnings(
    recipes::bake(prepped, new_data = input)
  )

  expect_length(captured$warnings, 1)
  diagnostics <- captured$warnings[[1]]$diagnostics
  expect_identical(diagnostics$row, 1:3)
  expect_identical(diagnostics$negative_points, c(0L, 1L, 1L))
  expect_equal(diagnostics$positive_area, c(60, 30, 30))
  expect_equal(diagnostics$negative_area, c(0, 100, 10))
  expect_equal(diagnostics$excluded_fraction, c(0, 10 / 13, 1 / 4))
  expect_identical(
    lapply(captured$value$ri, function(m) m$value),
    list(c(2, 2, 2), c(0, 4, 0), c(0, 4, 0))
  )

  tied <- sec_diagnostics_data(ri = list(curves[[1]], curves[[1]]))
  tied_result <- sec_capture_clipping_warnings(
    recipes::bake(prepped, new_data = tied)
  )
  expect_length(tied_result$warnings, 1)
  expect_match(conditionMessage(tied_result$warnings[[1]]), "row 1")

  many_rows <- sec_diagnostics_data(ri = rep(curves[1], 50))
  many_result <- sec_capture_clipping_warnings(
    recipes::bake(prepped, new_data = many_rows)
  )
  expect_length(many_result$warnings, 1)
  expect_identical(nrow(many_result$warnings[[1]]$diagnostics), 50L)
  expect_length(
    strsplit(conditionMessage(many_result$warnings[[1]]), "\n")[[1]],
    3
  )
})

test_that("bakes after serialization diagnose only the current columns and rows", {
  axis <- c(0, 1, 3)
  training <- sec_diagnostics_data(
    ri = list(sec_diagnostics_curve(axis, c(1, 2, 1))),
    uv = list(sec_diagnostics_curve(axis, c(2, 1, 2)))
  )
  prepped <- sec_diagnostics_recipe(training, measures = c("ri", "uv"))
  path <- withr::local_tempfile(fileext = ".rds")
  saveRDS(prepped, path)
  prepped <- readRDS(path)

  input <- sec_diagnostics_data(
    ri = list(
      sec_diagnostics_curve(axis, c(-2, 2, 0)),
      sec_diagnostics_curve(axis, c(2, 2, 2))
    ),
    uv = list(
      sec_diagnostics_curve(axis, c(1, 1, 1)),
      sec_diagnostics_curve(axis, c(-20, 2, 0))
    )
  )
  captured <- sec_capture_clipping_warnings(
    recipes::bake(prepped, new_data = input)
  )

  expect_length(captured$warnings, 2)
  expect_identical(
    vapply(captured$warnings, `[[`, character(1), "column"),
    c("ri", "uv")
  )
  expect_identical(captured$warnings[[1]]$diagnostics$row, 1:2)
  expect_identical(captured$warnings[[2]]$diagnostics$row, 1:2)
  expect_identical(
    captured$warnings[[1]]$diagnostics$negative_points,
    c(1L, 0L)
  )
  expect_identical(
    captured$warnings[[2]]$diagnostics$negative_points,
    c(0L, 1L)
  )

  positive_only <- sec_diagnostics_data(
    ri = list(sec_diagnostics_curve(axis, c(1, 2, 1))),
    uv = list(sec_diagnostics_curve(axis, c(2, 2, 2)))
  )
  second <- sec_capture_clipping_warnings(
    recipes::bake(prepped, new_data = positive_only)
  )
  expect_length(second$warnings, 0)
  expect_identical(nrow(second$value), 1L)
  expect_identical(attr(second$value$ri, "sec_concentration_unit"), "mg/mL")

  empty <- sec_diagnostics_data(ri = list(), uv = list())
  empty_result <- sec_capture_clipping_warnings(
    recipes::bake(prepped, new_data = empty)
  )
  expect_length(empty_result$warnings, 0)
  expect_identical(nrow(empty_result$value), 0L)
  expect_identical(
    attr(empty_result$value$ri, "sec_concentration_unit"),
    "mg/mL"
  )
  empty_diagnostics <- .sec_concentration_diagnostics(list())
  expect_named(
    empty_diagnostics,
    c(
      "row",
      "negative_points",
      "points",
      "positive_area",
      "negative_area",
      "excluded_fraction"
    )
  )
  expect_identical(empty_diagnostics$row, integer())
  expect_identical(empty_diagnostics$negative_points, integer())
  expect_identical(empty_diagnostics$points, integer())
  expect_identical(empty_diagnostics$positive_area, double())
  expect_identical(empty_diagnostics$negative_area, double())
  expect_identical(empty_diagnostics$excluded_fraction, double())
})

test_that("mass diagnostics use a scaled fraction for extreme finite areas", {
  axis <- c(0, 1, 2)
  training <- sec_diagnostics_data(
    ri = list(sec_diagnostics_curve(axis, c(1, 2, 1)))
  )
  prepped <- sec_diagnostics_recipe(
    training,
    measures = "ri",
    injection_mass = 1
  )
  input <- sec_diagnostics_data(
    ri = list(sec_diagnostics_curve(axis, c(-1.2e308, 1.2e308, -1.2e308)))
  )

  captured <- sec_capture_clipping_warnings(
    recipes::bake(prepped, new_data = input)
  )
  diagnostics <- captured$warnings[[1]]$diagnostics

  expect_identical(diagnostics$positive_area, 1.2e308)
  expect_identical(diagnostics$negative_area, 1.2e308)
  expect_identical(diagnostics$excluded_fraction, 0.5)
  expect_identical(captured$value$ri[[1]]$value, c(0, 1, 0))
})

test_that("normalization failures are reported before clipping warnings", {
  axis <- c(0, 1, 2)
  training <- sec_diagnostics_data(
    ri = list(sec_diagnostics_curve(axis, c(1, 2, 1)))
  )
  prepped <- sec_diagnostics_recipe(training, measures = "ri")

  invalid <- list(
    all_negative = list(
      curves = list(sec_diagnostics_curve(axis, c(-2, -1, -3))),
      mass = 6,
      flow = 1
    ),
    all_zero = list(
      curves = list(sec_diagnostics_curve(axis, c(0, 0, 0))),
      mass = 6,
      flow = 1
    ),
    late_row_all_negative = list(
      curves = list(
        sec_diagnostics_curve(axis, c(-2, 2, 0)),
        sec_diagnostics_curve(axis, c(-2, -1, -3))
      ),
      mass = 6,
      flow = 1
    ),
    positive_area_overflow = list(
      curves = list(sec_diagnostics_curve(c(0, 1), c(1e308, 1e308))),
      mass = 6,
      flow = 1
    ),
    negative_area_overflow = list(
      curves = list(sec_diagnostics_curve(axis, c(-1e308, -1e308, 1))),
      mass = 6,
      flow = 1
    ),
    axis_width_overflow = list(
      curves = list(
        sec_diagnostics_curve(c(-1e308, 1e308), c(1, 1))
      ),
      mass = 6,
      flow = 1
    ),
    denominator_overflow = list(
      curves = list(sec_diagnostics_curve(c(0, 1), c(2, 2))),
      mass = 6,
      flow = 1e308
    ),
    denominator_underflow = list(
      curves = list(sec_diagnostics_curve(c(0, 1), c(1e-200, 1e-200))),
      mass = 6,
      flow = 1e-200
    ),
    factor_underflow = list(
      curves = list(sec_diagnostics_curve(c(0, 1), c(0, 1e308))),
      mass = 1e-308,
      flow = 1
    ),
    factor_overflow = list(
      curves = list(sec_diagnostics_curve(c(0, 1), c(0, 1e-308))),
      mass = 1e308,
      flow = 1
    ),
    output_overflow = list(
      curves = list(sec_diagnostics_curve(c(0, 1), c(1e308, 0))),
      mass = 1e308,
      flow = 1
    ),
    output_underflow = list(
      curves = list(sec_diagnostics_curve(c(0, 1), c(1e-308, 1e-308))),
      mass = 1e-308,
      flow = 1
    )
  )

  errors <- lapply(names(invalid), function(name) {
    case <- invalid[[name]]
    step <- prepped
    step$steps[[1]]$injection_mass <- case$mass
    step$steps[[1]]$flow_rate <- case$flow
    captured <- new.env(parent = emptyenv())
    captured$warnings <- list()
    error <- tryCatch(
      withCallingHandlers(
        recipes::bake(
          step,
          new_data = sec_diagnostics_data(ri = case$curves)
        ),
        measure_sec_negative_signal_clipped = function(warning) {
          captured$warnings[[length(captured$warnings) + 1L]] <- warning
          invokeRestart("muffleWarning")
        }
      ),
      error = identity
    )

    expect_length(captured$warnings, 0)
    conditionMessage(error)
  })
  names(errors) <- names(invalid)
  expect_snapshot({
    for (name in names(errors)) {
      cat(name, ":\n", sep = "")
      cat(errors[[name]], "\n", sep = "")
    }
  })
})

test_that("positive tails may underflow individually and no-op mode stays quiet", {
  axis <- c(0, 1)
  training <- sec_diagnostics_data(
    ri = list(sec_diagnostics_curve(axis, c(1, 1)))
  )
  prepped <- sec_diagnostics_recipe(
    training,
    measures = "ri",
    concentration_units = "g/mL"
  )
  tiny_tail <- sec_diagnostics_data(
    ri = list(sec_diagnostics_curve(axis, c(1, 1e-323)))
  )

  tail_result <- sec_capture_clipping_warnings(
    recipes::bake(prepped, new_data = tiny_tail)
  )
  expect_length(tail_result$warnings, 0)
  expect_gt(tail_result$value$ri[[1]]$value[[1]], 0)
  expect_identical(tail_result$value$ri[[1]]$value[[2]], 0)

  no_op <- sec_diagnostics_recipe(
    training,
    measures = "ri",
    normalize_to_mass = FALSE
  )
  negative <- sec_diagnostics_data(
    ri = list(sec_diagnostics_curve(axis, c(-1, 2)))
  )
  no_op_result <- sec_capture_clipping_warnings(
    recipes::bake(no_op, new_data = negative)
  )
  expect_length(no_op_result$warnings, 0)
  expect_identical(no_op_result$value, negative)
})
