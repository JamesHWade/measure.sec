sec_baseline_support_curve <- function(
  values,
  location = seq_along(values)
) {
  measure::new_measure_tbl(location = location, value = values)
}

sec_baseline_support_data <- function(
  curve,
  column = "ri",
  location = seq_along(curve)
) {
  data <- tibble::tibble(sample_id = "sample")
  data[[column]] <- measure::new_measure_list(list(
    sec_baseline_support_curve(curve, location)
  ))
  data
}

sec_baseline_support_bake <- function(
  curve,
  left_frac = 0.2,
  right_frac = 0.2,
  method = "linear",
  column = "ri",
  location = seq_along(curve)
) {
  training <- sec_baseline_support_data(
    seq_along(curve),
    column = column,
    location = location
  )
  data <- sec_baseline_support_data(
    curve,
    column = column,
    location = location
  )
  recipe <- recipes::recipe(training) |>
    step_sec_baseline(
      measures = column,
      left_frac = left_frac,
      right_frac = right_frac,
      method = method
    ) |>
    recipes::prep()
  recipes::bake(recipe, new_data = data)
}

expect_sec_baseline_support_error <- function(
  object,
  missing_edges,
  column = "ri"
) {
  while (
    !inherits(object, "measure_sec_baseline_support") &&
      !is.null(object$parent)
  ) {
    object <- object$parent
  }
  expect_s3_class(object, "measure_sec_baseline_support")
  expect_identical(object$missing_edges, missing_edges)
  expect_identical(object$column, column)
  expect_identical(object$curve, 1L)
  object
}

test_that("delay output reports missing baseline support", {
  n <- 201L
  location <- seq(0, 2, by = 0.01)
  data <- tibble::tibble(sample_id = "sample")
  data$reference <- measure::new_measure_list(list(
    sec_baseline_support_curve(rep(0, n), location)
  ))
  data$target <- measure::new_measure_list(list(
    sec_baseline_support_curve(
      dnorm(location, mean = 1, sd = 0.2) + 0.01 * location,
      location
    )
  ))

  delayed_recipe <- recipes::recipe(data) |>
    step_sec_detector_delay(
      reference = "reference",
      targets = "target",
      delay_volumes = c(target = 0.15),
      method = "interpolate"
    ) |>
    recipes::prep()
  delayed <- recipes::bake(delayed_recipe, new_data = data)

  expect_length(delayed$target[[1]]$location, n)
  expect_equal(sum(is.finite(delayed$target[[1]]$value)), 186L)
  expect_equal(sum(is.na(delayed$target[[1]]$value)), 15L)

  baseline_recipe <- recipes::recipe(data) |>
    step_sec_baseline(
      measures = "target",
      left_frac = 0.05,
      right_frac = 0.05
    ) |>
    recipes::prep()
  error <- tryCatch(
    recipes::bake(baseline_recipe, new_data = delayed),
    error = identity
  )
  expect_sec_baseline_support_error(
    error,
    missing_edges = "right",
    column = "target"
  )
})

test_that("edge support errors identify missing sides", {
  n <- 20L
  finite <- seq_len(n)
  cases <- list(
    left = replace(finite, 1:4, NA_real_),
    right = replace(finite, 17:20, NA_real_),
    both = replace(replace(finite, 1:4, NA_real_), 17:20, NA_real_)
  )

  for (method in c("linear", "median", "spline")) {
    for (edge in names(cases)) {
      error <- tryCatch(
        sec_baseline_support_bake(
          cases[[edge]],
          left_frac = 0.2,
          right_frac = 0.2,
          method = method
        ),
        error = identity
      )
      missing_edges <- if (edge == "both") c("left", "right") else edge
      support_error <- expect_sec_baseline_support_error(error, missing_edges)
      expect_identical(support_error$method, method)
    }
  }
})

test_that("linear fits finite asymmetric support", {
  n <- 20L
  location <- seq(10, 29)
  values <- 5 + 0.2 * seq_len(n)
  values[c(1, 3, 17, 19)] <- c(NA_real_, NaN, Inf, -Inf)

  result <- sec_baseline_support_bake(
    values,
    left_frac = 0.25,
    right_frac = 0.2,
    location = location
  )
  corrected <- result$ri[[1]]$value

  expect_equal(
    corrected[is.finite(values)],
    rep(0, sum(is.finite(values))),
    tolerance = 1e-12
  )
  expect_identical(corrected[!is.finite(values)], values[!is.finite(values)])
  expect_identical(result$ri[[1]]$location, location)

  control <- sec_baseline_support_bake(
    5 + 0.2 * seq_len(n),
    left_frac = 0.25,
    right_frac = 0.2
  )
  expect_equal(control$ri[[1]]$value, rep(0, n), tolerance = 1e-12)
})

test_that("median uses the supported index", {
  location <- seq_len(20L)
  values <- 2 + 0.5 * location
  values[c(3, 4)] <- NA_real_

  result <- sec_baseline_support_bake(
    values,
    left_frac = 0.25,
    right_frac = 0.15,
    method = "median",
    location = location
  )
  corrected <- result$ri[[1]]$value

  expect_equal(
    corrected[is.finite(values)],
    rep(0, sum(is.finite(values))),
    tolerance = 1e-12
  )
  expect_identical(corrected[c(3, 4)], values[c(3, 4)])
})

test_that("spline falls back with too little support", {
  values <- 5 + 0.2 * seq_len(20L)
  values[c(1, 20)] <- c(NA_real_, NA_real_)

  linear <- sec_baseline_support_bake(
    values,
    left_frac = 0.1,
    right_frac = 0.1,
    method = "linear"
  )
  expect_warning(
    spline <- sec_baseline_support_bake(
      values,
      left_frac = 0.1,
      right_frac = 0.1,
      method = "spline"
    ),
    "Too few valid baseline points"
  )
  expect_equal(spline$ri[[1]]$value, linear$ri[[1]]$value)
})

test_that("spline fits supported edges", {
  location <- seq_len(40L)
  values <- 0.5 + 0.02 * (location - 20)^2

  linear <- sec_baseline_support_bake(
    values,
    left_frac = 0.2,
    right_frac = 0.2,
    method = "linear",
    location = location
  )
  expect_no_warning({
    spline <- sec_baseline_support_bake(
      values,
      left_frac = 0.2,
      right_frac = 0.2,
      method = "spline",
      location = location
    )
  })

  expect_true(
    max(abs(spline$ri[[1]]$value - linear$ri[[1]]$value)) > 0.1
  )
  expect_length(spline$ri[[1]]$value, length(location))
})

test_that("missing states survive correction", {
  location <- seq(100, 119)
  values <- 2 + 0.1 * seq_along(location)
  missing <- c(NA_real_, NaN, Inf, -Inf)
  values[c(8, 9, 12, 13)] <- missing

  result <- sec_baseline_support_bake(
    values,
    location = location,
    left_frac = 0.2,
    right_frac = 0.2
  )
  corrected <- result$ri[[1]]

  expect_identical(corrected$location, location)
  expect_identical(corrected$value[c(8, 9, 12, 13)], missing)
  expect_equal(
    corrected$value[is.finite(values)],
    rep(0, sum(is.finite(values))),
    tolerance = 1e-12
  )
})
