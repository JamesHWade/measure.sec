make_sec_ri_output_test_data <- function(
  values = list(c(-2, 1), c(3, -4)),
  dn_dc = c(-0.25, 0.5)
) {
  sample_id <- seq_along(values)
  data <- tibble::tibble(
    sample_id = sample_id,
    dn_dc = dn_dc[sample_id]
  )
  data$ri <- measure::new_measure_list(lapply(seq_along(values), function(i) {
    measure::new_measure_tbl(
      location = seq_along(values[[i]]),
      value = values[[i]]
    )
  }))
  attr(data$ri, "source_marker") <- list(source = "raw detector", id = 75L)
  data
}

test_that("RI output columns preserve the source and survive fresh recipe bakes", {
  values <- list(c(-2, 1), c(3, -4), c(-1, 2))
  dn_dc <- c(-0.25, 0.5, -0.2)
  instrument_constant <- 2
  data <- make_sec_ri_output_test_data(values, dn_dc)

  prepped <- recipes::recipe(data) |>
    recipes::update_role(ri, new_role = "outcome") |>
    step_sec_ri(
      measures = "ri",
      dn_dc_column = "dn_dc",
      instrument_constant = instrument_constant,
      concentration_units = "mg/mL",
      output_col = "converted RI",
      role = "predictor"
    ) |>
    recipes::prep()

  step_info <- recipes::tidy(prepped, number = 1)
  expect_identical(step_info$output_col, "converted RI")
  variable_info <- prepped$last_term_info
  source_role <- variable_info$role[[match("ri", variable_info$variable)]]
  output_role <- variable_info$role[[match(
    "converted RI",
    variable_info$variable
  )]]
  expect_identical(source_role, "outcome")
  expect_identical(output_role, "predictor")

  serialized_recipe <- withr::local_tempfile()
  saveRDS(prepped, serialized_recipe)
  restored <- readRDS(serialized_recipe)

  row_order <- c(3, 1, 2)
  bake_inputs <- list(data, data[row_order, ])
  for (input in bake_inputs) {
    baked <- recipes::bake(restored, new_data = input)
    indices <- if (identical(input$sample_id, data$sample_id)) {
      seq_len(nrow(data))
    } else {
      row_order
    }
    expected <- Map(
      function(value, dndc) value * instrument_constant / dndc * 1000,
      values[indices],
      dn_dc[indices]
    )
    expect_equal(
      lapply(baked[["converted RI"]], function(curve) curve$value),
      expected
    )
    expect_identical(
      attr(baked[["converted RI"]], "sec_concentration_unit"),
      "mg/mL"
    )
    expect_identical(baked$ri, input$ri)
  }
  expect_identical(
    attr(data$ri, "source_marker"),
    list(source = "raw detector", id = 75L)
  )
})

test_that("RI raw output has no concentration unit and a missing role stays missing", {
  values <- list(c(-2, 1), c(3, -4))
  data <- make_sec_ri_output_test_data(values)
  prepped <- recipes::recipe(data) |>
    step_sec_ri(
      measures = "ri",
      instrument_constant = 3,
      output_col = "raw scaled"
    ) |>
    recipes::prep()

  output_role <- prepped$last_term_info$role[[
    match("raw scaled", prepped$last_term_info$variable)
  ]]
  expect_true(is.na(output_role))

  result <- recipes::bake(prepped, new_data = data)
  expect_equal(
    lapply(result[["raw scaled"]], function(curve) curve$value),
    lapply(values, `*`, 3)
  )
  expect_null(attr(result[["raw scaled"]], "sec_concentration_unit"))
  expect_identical(result$ri, data$ri)

  empty_result <- recipes::bake(prepped, new_data = data[0, ])
  expect_identical(
    names(empty_result),
    c("sample_id", "dn_dc", "ri", "raw scaled")
  )
  expect_length(empty_result[["raw scaled"]], 0)
  expect_null(attr(empty_result[["raw scaled"]], "sec_concentration_unit"))
})

test_that("RI output defaults and an explicit source name both process in place", {
  data <- make_sec_ri_output_test_data()
  default_output <- recipes::recipe(data) |>
    step_sec_ri(measures = "ri", dn_dc = -0.5) |>
    recipes::prep()
  named_output <- recipes::recipe(data) |>
    step_sec_ri(measures = "ri", dn_dc = -0.5, output_col = "ri") |>
    recipes::prep()

  default_result <- recipes::bake(default_output, new_data = data)
  named_result <- recipes::bake(named_output, new_data = data)
  expect_named(default_result, names(data))
  expect_named(named_result, names(data))
  expect_identical(default_result$ri, named_result$ri)
  expect_identical(
    attr(default_result$ri, "sec_concentration_unit"),
    "g/mL"
  )
  expect_identical(
    recipes::tidy(named_output, number = 1)$output_col,
    "ri"
  )
  expect_true(is.na(recipes::tidy(default_output, number = 1)$output_col))
})

test_that("RI output names and collisions are validated before writing", {
  recipe <- recipes::recipe(make_sec_ri_output_test_data())
  invalid_names <- list(
    1,
    character(),
    c("one", "two"),
    NA_character_,
    "",
    "   "
  )
  for (output_col in invalid_names) {
    expect_error(step_sec_ri(recipe, output_col = output_col), "output_col")
  }

  data <- make_sec_ri_output_test_data()
  data$occupied <- 1
  expect_error(
    recipes::recipe(data) |>
      step_sec_ri(measures = "ri", output_col = "occupied") |>
      recipes::prep(),
    "already exists"
  )

  prepped <- recipes::recipe(data[c("sample_id", "dn_dc", "ri")]) |>
    step_sec_ri(measures = "ri", output_col = "occupied") |>
    recipes::prep()
  fresh_data <- data[c("sample_id", "dn_dc", "ri", "occupied")]
  expect_error(
    recipes::bake(
      prepped$steps[[1]],
      new_data = fresh_data
    ),
    "already exists"
  )
})

test_that("in-place RI processing keeps multiple measures while named output does not", {
  data <- make_sec_ri_output_test_data()
  data$ri_extra <- data$ri
  values <- lapply(data$ri, function(curve) curve$value)

  prepped_in_place <- recipes::recipe(data) |>
    step_sec_ri(measures = c("ri", "ri_extra"), dn_dc = -0.5) |>
    recipes::prep()
  result_in_place <- recipes::bake(prepped_in_place, new_data = data)
  expected <- lapply(values, function(value) value / -0.5)
  expect_equal(
    lapply(result_in_place$ri, function(curve) curve$value),
    expected
  )
  expect_equal(
    lapply(result_in_place$ri_extra, function(curve) curve$value),
    expected
  )

  inferred_in_place <- recipes::recipe(data) |>
    step_sec_ri(dn_dc = -0.5) |>
    recipes::prep()
  expect_identical(
    inferred_in_place$steps[[1]]$measures,
    c("ri", "ri_extra")
  )

  expect_error(
    recipes::recipe(data) |>
      step_sec_ri(
        measures = c("ri", "ri_extra"),
        dn_dc = -0.5,
        output_col = "converted RI"
      ) |>
      recipes::prep(),
    "exactly one RI measure"
  )

  expect_error(
    recipes::recipe(data) |>
      step_sec_ri(output_col = "converted RI") |>
      recipes::prep(),
    "exactly one RI measure"
  )

  prepped <- recipes::recipe(data[c("sample_id", "dn_dc", "ri")]) |>
    step_sec_ri(measures = "ri") |>
    recipes::prep()
  prepped$steps[[1]]$measures <- c("ri", "ri_extra")
  prepped$steps[[1]]$output_col <- "converted RI"
  expect_error(
    recipes::bake(prepped, new_data = data),
    "exactly one RI measure"
  )
})
