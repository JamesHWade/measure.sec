test_that("MW slices interpolate selected support boundaries", {
  signal <- list(
    location = c(1, 2, 4),
    value = c(1, 3, 1)
  )
  mw_curve <- list(
    location = c(1, 2, 4),
    value = c(100, 200, 400)
  )

  slices <- .sec_mw_slices(
    signal,
    mw_curve,
    mw_scale = "absolute",
    integration_range = c(1.5, 3)
  )

  expect_equal(slices$location, c(1.5, 2, 3))
  expect_equal(slices$signal, c(2, 3, 2))
  expect_equal(slices$mw, c(150, 200, 300))
  expect_equal(slices$weight, c(0.5, 2.25, 1))
  expect_identical(slices$segment, c(1L, 1L, 1L))
})

test_that("MW interpolation respects the declared scale", {
  signal <- list(location = c(1, 2), value = c(1, 1))
  log_curve <- list(
    location = c(1, 2),
    value = log10(c(100, 200))
  )

  log_slices <- .sec_mw_slices(
    signal,
    log_curve,
    mw_scale = "log10",
    integration_range = c(1.5, 2)
  )
  absolute_slices <- .sec_mw_slices(
    signal,
    list(location = c(1, 2), value = c(100, 200)),
    mw_scale = "absolute",
    integration_range = c(1.5, 2)
  )

  expect_equal(log_slices$mw[1], sqrt(100 * 200))
  expect_equal(absolute_slices$mw[1], 150)
})

test_that("legacy log-MW and linear calibration map clipped bounds", {
  signal <- list(location = c(1, 2), value = c(1, 1))

  axis_log_mw <- .sec_mw_slices(
    signal,
    integration_range = c(1.5, 2)
  )
  linear_calibration <- .sec_mw_slices(
    signal,
    calibration = c(2, -1),
    integration_range = c(1.5, 2)
  )

  expect_equal(axis_log_mw$mw, c(10^1.5, 100))
  expect_equal(linear_calibration$mw, c(100, 1000))
})

test_that("MW averages integrate exact window boundaries", {
  result <- .calc_mw_averages(
    location = c(1, 2, 4),
    value = c(1, 2, 1),
    calibration = NULL,
    integration_range = c(1, 2),
    output_cols = "mw",
    mw_measure = list(
      location = c(1, 2, 4),
      value = c(100000, 200000, 400000)
    ),
    mw_scale = "absolute"
  )

  expect_equal(unname(result[["mw"]]), 166666.6667, tolerance = 1e-4)
})

test_that("integration exclusions split support and merge overlaps", {
  slices <- .sec_mw_slices(
    signal = list(location = 1:10, value = rep(1, 10)),
    excluded_regions = tibble::tibble(
      start = c(2, 2.5, 6),
      end = c(3, 4, 7),
      purpose = c("integration", "both", "integration")
    )
  )

  expect_equal(slices$location, c(1, 2, 4, 5, 6, 7, 8, 9, 10))
  expect_identical(
    slices$segment,
    c(1L, 1L, 2L, 2L, 2L, 3L, 3L, 3L, 3L)
  )
  expect_equal(sum(slices$weight), 6)
  expect_equal(slices$weight[c(2, 3, 6)], c(0.5, 0.5, 0.5))
})

test_that("baseline-only exclusions do not change MW support", {
  slices <- .sec_mw_slices(
    signal = list(location = 1:10, value = rep(1, 10)),
    excluded_regions = tibble::tibble(
      start = 2,
      end = 8,
      purpose = "baseline"
    )
  )

  expect_equal(slices$location, 1:10)
  expect_identical(slices$segment, rep(1L, 10))
  expect_equal(sum(slices$weight), 9)
})

test_that("MW integration preserves descending irregular axes", {
  slices <- .sec_mw_slices(
    signal = list(location = c(5, 2, 1), value = c(1, 3, 2)),
    integration_range = c(1.5, 4)
  )

  expect_equal(slices$location, c(4, 2, 1.5))
  expect_equal(slices$signal, c(1.6666667, 3, 2.5), tolerance = 1e-6)
  expect_equal(slices$weight, c(1.6666667, 3.75, 0.625), tolerance = 1e-6)
  expect_identical(slices$segment, c(1L, 1L, 1L))
})

test_that("invalid MW in an excluded interval is ignored", {
  slices <- .sec_mw_slices(
    signal = list(location = 1:4, value = rep(1, 4)),
    mw_curve = list(
      location = 1:4,
      value = c(100, NA_real_, 300, 400)
    ),
    mw_scale = "absolute",
    excluded_regions = tibble::tibble(
      start = 1.5,
      end = 2.5,
      purpose = "integration"
    )
  )

  expect_equal(slices$location, c(1, 1.5, 2.5, 3, 4))
  expect_equal(slices$mw, c(100, 150, 250, 300, 400))
  expect_identical(slices$segment, c(1L, 1L, 2L, 2L, 2L))
})

test_that("positive MW support with a calibration gap errors", {
  expect_snapshot(error = TRUE, {
    .sec_mw_slices(
      signal = list(location = 1:4, value = rep(1, 4)),
      mw_curve = list(
        location = 1:4,
        value = c(100, NA_real_, 300, 400)
      ),
      mw_scale = "absolute"
    )
  })
})

test_that("off-grid support outside valid calibration endpoints errors", {
  expect_snapshot(error = TRUE, {
    .sec_mw_slices(
      signal = list(location = 0:2, value = rep(1, 3)),
      mw_curve = list(location = 0:2, value = c(NA_real_, 1, 2)),
      mw_scale = "absolute",
      integration_range = c(0.5, 2)
    )
  })

  expect_snapshot(error = TRUE, {
    .sec_mw_slices(
      signal = list(location = 0:2, value = rep(1, 3)),
      mw_curve = list(location = 0:2, value = c(1, 2, NA_real_)),
      mw_scale = "absolute",
      integration_range = c(0, 1.5)
    )
  })
})

test_that("zero-width support does not require an MW value", {
  slices <- .sec_mw_slices(
    signal = list(location = 1, value = 1),
    mw_curve = list(location = 1, value = NA_real_),
    mw_scale = "absolute"
  )

  expect_equal(slices$weight, 0)
  expect_true(is.na(slices$mw))
})

test_that("row integration metadata uses overrides and validates metadata", {
  new_data <- tibble::tibble(
    .integration_window = list(tibble::tibble(start = 1, end = 2)),
    .excluded_regions = list(tibble::tibble(
      start = c(1.25, 1.5),
      end = c(1.4, 1.75),
      purpose = c("baseline", "integration")
    ))
  )

  metadata <- .sec_row_integration_metadata(new_data, 1)
  expect_equal(metadata$integration_range, c(1, 2))
  expect_equal(metadata$excluded_regions$start, 1.5)
  expect_equal(metadata$excluded_regions$end, 1.75)

  absent <- .sec_row_integration_metadata(tibble::tibble(sample_id = "s"), 1)
  null_window <- .sec_row_integration_metadata(
    tibble::tibble(.integration_window = list(NULL)),
    1
  )
  expect_null(absent$integration_range)
  expect_null(null_window$integration_range)

  override <- .sec_row_integration_metadata(
    new_data,
    1,
    integration_range = c(2, 3)
  )
  expect_equal(override$integration_range, c(2, 3))

  malformed <- new_data
  malformed$.integration_window[[1]] <- tibble::tibble(start = 2, end = 2)
  expect_snapshot(error = TRUE, {
    .sec_row_integration_metadata(malformed, 1)
  })

  malformed <- new_data
  malformed$.excluded_regions[[1]]$purpose[2] <- "unknown"
  expect_snapshot(error = TRUE, {
    .sec_row_integration_metadata(malformed, 1)
  })
})

test_that("paired MW axes validate before support selection", {
  expect_snapshot(error = TRUE, {
    .sec_mw_slices(
      signal = list(location = 1:3, value = rep(1, 3)),
      mw_curve = list(location = c(1, 2, 4), value = c(100, 200, 400)),
      mw_scale = "absolute",
      integration_range = c(10, 11)
    )
  })
})

test_that("MW steps read row metadata after recipe serialization and reorder", {
  locations <- 1:4
  make_curve <- function(value) {
    measure::new_measure_tbl(location = locations, value = value)
  }
  training <- tibble::tibble(sample_id = c("A", "B"))
  training$ri <- measure::new_measure_list(list(
    make_curve(rep(1, 4)),
    make_curve(rep(1, 4))
  ))
  training$mw_curve <- measure::new_measure_list(list(
    make_curve(c(100, 200, 300, 400)),
    make_curve(c(100, 200, 300, 400))
  ))

  exclusions <- tibble::tibble(
    start = c(2, 1),
    end = c(3, 3),
    sample_id = c("A", "B")
  )
  recipe <- recipes::recipe(~., data = training) |>
    step_sec_exclude_regions(
      measures = "ri",
      regions = exclusions,
      purpose = "integration"
    ) |>
    step_sec_mw_averages(
      measures = "ri",
      mw_col = "mw_curve",
      mw_scale = "absolute",
      output_cols = "mw"
    ) |>
    step_sec_mw_fractions(
      measures = "ri",
      mw_col = "mw_curve",
      mw_scale = "absolute",
      cutoffs = 250
    )
  prepped <- recipes::prep(recipe)
  path <- tempfile()
  on.exit(unlink(path), add = TRUE)
  saveRDS(prepped, path)
  prepped <- readRDS(path)

  result <- recipes::bake(prepped, new_data = training[2:1, ])

  expect_equal(as.character(result$sample_id), c("B", "A"))
  expect_equal(result$mw_mw, c(350, 250))
  expect_equal(result$frac_below_250, c(0, 0.5))
})

test_that("integration-window producer supplies MW defaults and override", {
  locations <- c(1, 2, 4)
  data <- tibble::tibble(sample_id = "sample")
  data$ri <- measure::new_measure_list(list(
    measure::new_measure_tbl(location = locations, value = c(1, 2, 1))
  ))
  data$mw_curve <- measure::new_measure_list(list(
    measure::new_measure_tbl(
      location = locations,
      value = c(100000, 200000, 400000)
    )
  ))

  base_recipe <- recipes::recipe(~., data = data) |>
    step_sec_integration_window(
      measures = "ri",
      start = 1,
      end = 2,
      auto_detect = FALSE
    )
  row_window <- base_recipe |>
    step_sec_mw_averages(
      measures = "ri",
      mw_col = "mw_curve",
      mw_scale = "absolute",
      output_cols = "mw"
    ) |>
    recipes::prep() |>
    recipes::bake(new_data = NULL)
  override <- base_recipe |>
    step_sec_mw_averages(
      measures = "ri",
      mw_col = "mw_curve",
      mw_scale = "absolute",
      integration_range = c(2, 4),
      output_cols = "mw"
    ) |>
    recipes::prep() |>
    recipes::bake(new_data = NULL)

  expect_equal(row_window$mw_mw, 166666.6667, tolerance = 1e-4)
  expect_equal(override$mw_mw, 266666.6667, tolerance = 1e-4)
})

test_that("empty MW support yields NA moments and fractions", {
  average <- .calc_mw_averages(
    location = 1:3,
    value = c(1, 2, 1),
    calibration = NULL,
    integration_range = c(5, 6),
    output_cols = "mw"
  )
  fractions <- .calc_mw_fractions(
    location = 1:3,
    value = c(1, 2, 1),
    cutoffs = 100,
    calibration = NULL,
    integration_range = c(5, 6)
  )

  expect_true(is.na(average[["mw"]]))
  expect_true(all(is.na(fractions)))
})
