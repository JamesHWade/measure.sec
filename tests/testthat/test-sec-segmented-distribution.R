# Focused tests for segmented SEC molecular-weight distributions.

mw_curve_area <- function(curve) {
  if (length(curve$location) < 2) {
    return(0)
  }

  sum(
    diff(curve$location) *
      (head(curve$value, -1) + tail(curve$value, -1)) /
      2
  )
}

mw_curve_at <- function(curve, location) {
  stats::approx(curve$location, curve$value, xout = location)$y
}

mw_empty_exclusions <- function() {
  tibble::tibble(
    start = numeric(),
    end = numeric(),
    purpose = character()
  )
}

mw_distribution_data <- function(location, signals, sample_id = "sample") {
  if (!is.list(signals)) {
    signals <- rep(list(signals), length(sample_id))
  }

  data <- tibble::tibble(sample_id = sample_id)
  data$ri <- measure::new_measure_list(lapply(signals, function(signal) {
    measure::new_measure_tbl(location, signal)
  }))
  data
}

test_that("the default unexcluded distribution grid remains unchanged", {
  skip_if_not_installed("measure")

  location <- seq(3, 6, by = 0.25)
  training <- mw_distribution_data(location, rep(2, length(location)))
  prepped <- recipes::recipe(training) |>
    step_sec_mw_distribution(
      measures = "ri",
      type = "both",
      n_points = 10,
      normalize = FALSE
    ) |>
    recipes::prep()
  metadata_training <- training
  metadata_training$.excluded_regions <- list(mw_empty_exclusions())
  metadata_prepped <- recipes::recipe(metadata_training) |>
    step_sec_mw_distribution(
      measures = "ri",
      type = "both",
      n_points = 10,
      normalize = FALSE
    ) |>
    recipes::prep()

  fresh <- mw_distribution_data(location, rep(2, length(location)))
  without_metadata <- recipes::bake(prepped, new_data = fresh)
  with_empty_exclusions <- fresh
  with_empty_exclusions$.excluded_regions <- list(mw_empty_exclusions())
  with_metadata <- recipes::bake(
    metadata_prepped,
    new_data = with_empty_exclusions
  )

  expect_identical(without_metadata$ri[[1]], with_metadata$ri[[1]])
  expect_identical(
    without_metadata$ri_cumulative[[1]],
    with_metadata$ri_cumulative[[1]]
  )
  expect_length(without_metadata$ri[[1]]$location, 10)
  expect_equal(
    diff(without_metadata$ri[[1]]$location),
    rep(1 / 3, 9),
    tolerance = 1e-10
  )
})

test_that("segmented distributions conserve mass and plateau across exclusions", {
  skip_if_not_installed("measure")

  location <- 0:10
  data <- mw_distribution_data(location, rep(2, length(location)))
  rec <- recipes::recipe(data) |>
    step_sec_exclude_regions(
      regions = tibble::tibble(start = 3.25, end = 5.75),
      purpose = "integration"
    ) |>
    step_sec_mw_distribution(
      measures = "ri",
      type = "both",
      n_points = 10,
      normalize = FALSE
    )

  result <- recipes::bake(recipes::prep(rec), new_data = NULL)
  differential <- result$ri[[1]]
  cumulative <- result$ri_cumulative[[1]]

  expect_equal(mw_curve_area(differential), 15, tolerance = 1e-10)
  expect_equal(
    mw_curve_at(cumulative, c(3.25, 4.5, 5.75)),
    rep(6.5 / 15, 3),
    tolerance = 1e-10
  )
  expect_equal(
    differential$value[
      differential$location > 3.25 &
        differential$location < 5.75
    ],
    rep(
      0,
      sum(
        differential$location > 3.25 &
          differential$location < 5.75
      )
    )
  )
  expect_equal(
    differential$value[differential$location == 3.25],
    c(2, 0),
    tolerance = 1e-10
  )
  expect_equal(
    differential$value[differential$location == 5.75],
    c(0, 2),
    tolerance = 1e-10
  )
})

test_that("segmented distributions preserve normalization and cropped mass", {
  skip_if_not_installed("measure")

  location <- 0:10
  data <- mw_distribution_data(location, rep(2, length(location)))
  build_recipe <- function(normalize) {
    recipes::recipe(data) |>
      step_sec_exclude_regions(
        regions = tibble::tibble(start = 3.25, end = 5.75),
        purpose = "integration"
      ) |>
      step_sec_mw_distribution(
        measures = "ri",
        type = "both",
        n_points = 10,
        mw_range = c(1e2, 1e8),
        normalize = normalize
      )
  }

  raw <- recipes::bake(recipes::prep(build_recipe(FALSE)), new_data = NULL)
  normalized <- recipes::bake(
    recipes::prep(build_recipe(TRUE)),
    new_data = NULL
  )

  expect_equal(mw_curve_area(raw$ri[[1]]), 7, tolerance = 1e-10)
  expect_equal(mw_curve_area(normalized$ri[[1]]), 7 / 15, tolerance = 1e-10)
  expect_equal(
    mw_curve_at(raw$ri_cumulative[[1]], c(2, 8)),
    c(4 / 15, 11 / 15),
    tolerance = 1e-10
  )
})

test_that("overlapping off-grid exclusions are merged without bridging", {
  skip_if_not_installed("measure")

  location <- 0:10
  data <- mw_distribution_data(location, rep(1, length(location)))
  rec <- recipes::recipe(data) |>
    step_sec_exclude_regions(
      regions = tibble::tibble(
        start = c(3.25, 4.5),
        end = c(5.5, 6.75)
      ),
      purpose = "integration"
    ) |>
    step_sec_mw_distribution(
      measures = "ri",
      type = "both",
      n_points = 10,
      normalize = FALSE
    )

  result <- recipes::bake(recipes::prep(rec), new_data = NULL)
  differential <- result$ri[[1]]
  cumulative <- result$ri_cumulative[[1]]

  expect_equal(mw_curve_area(differential), 6.5, tolerance = 1e-10)
  expect_equal(
    mw_curve_at(cumulative, c(3.25, 5, 6.75)),
    rep(3.25 / 6.5, 3),
    tolerance = 1e-10
  )
})

test_that("coarse grids retain narrow gaps and zero-mass segments", {
  skip_if_not_installed("measure")

  location <- 0:10
  narrow_data <- mw_distribution_data(location, rep(1, length(location)))
  narrow_recipe <- recipes::recipe(narrow_data) |>
    step_sec_exclude_regions(
      regions = tibble::tibble(start = 4.9, end = 5.1),
      purpose = "integration"
    ) |>
    step_sec_mw_distribution(
      measures = "ri",
      type = "both",
      n_points = 10,
      normalize = FALSE
    )
  narrow <- recipes::bake(
    recipes::prep(narrow_recipe),
    new_data = NULL
  )
  narrow_diff <- narrow$ri[[1]]
  narrow_cdf <- narrow$ri_cumulative[[1]]

  expect_equal(mw_curve_area(narrow_diff), 9.8, tolerance = 1e-10)
  expect_equal(
    mw_curve_at(narrow_cdf, c(4.9, 5, 5.1)),
    rep(0.5, 3),
    tolerance = 1e-10
  )
  expect_equal(
    narrow_diff$value[narrow_diff$location == 4.9],
    c(1, 0),
    tolerance = 1e-10
  )
  expect_equal(
    narrow_diff$value[narrow_diff$location == 5.1],
    c(0, 1),
    tolerance = 1e-10
  )

  zero_segment_data <- mw_distribution_data(
    location,
    c(rep(1, 5), rep(0, 6))
  )
  zero_segment_recipe <- recipes::recipe(zero_segment_data) |>
    step_sec_exclude_regions(
      regions = tibble::tibble(start = 4, end = 6),
      purpose = "integration"
    ) |>
    step_sec_mw_distribution(
      measures = "ri",
      type = "both",
      n_points = 10,
      normalize = FALSE
    )
  zero_segment <- recipes::bake(
    recipes::prep(zero_segment_recipe),
    new_data = NULL
  )

  expect_equal(mw_curve_area(zero_segment$ri[[1]]), 4, tolerance = 1e-10)
  expect_equal(
    mw_curve_at(zero_segment$ri_cumulative[[1]], c(4, 6, 8, 10)),
    rep(1, 4),
    tolerance = 1e-10
  )
  expect_true(all(is.finite(zero_segment$ri[[1]]$value)))
  expect_equal(
    zero_segment$ri[[1]]$value[zero_segment$ri[[1]]$location >= 6],
    rep(0, sum(zero_segment$ri[[1]]$location >= 6))
  )
})

test_that("disconnected components retain their own mass after smoothing", {
  skip_if_not_installed("measure")

  location <- c(0, 3.9, 4, 6, 6.1, 10)
  signal <- c(1, 1, 100, 1, 1, 1)
  data <- mw_distribution_data(location, signal)
  rec <- recipes::recipe(data) |>
    step_sec_exclude_regions(
      regions = tibble::tibble(start = 4, end = 6),
      purpose = "integration"
    ) |>
    step_sec_mw_distribution(
      measures = "ri",
      calibration = c(1, 5),
      type = "both",
      n_points = NULL,
      normalize = TRUE
    )

  result <- recipes::bake(recipes::prep(rec), new_data = NULL)
  differential <- result$ri[[1]]
  cumulative <- result$ri_cumulative[[1]]
  left_area <- mw_curve_area(differential[differential$location <= 9, ])
  right_area <- mw_curve_area(differential[differential$location >= 11, ])
  left_mass <- mw_curve_at(cumulative, 9)

  expect_equal(left_area, left_mass, tolerance = 1e-10)
  expect_equal(right_area, 1 - left_mass, tolerance = 1e-10)
})

test_that("excluded heads and tails stay zero outside retained support", {
  skip_if_not_installed("measure")

  location <- 0:10
  training <- mw_distribution_data(location, rep(1, length(location)))
  training$.excluded_regions <- list(mw_empty_exclusions())
  prepped <- recipes::recipe(training) |>
    step_sec_mw_distribution(
      measures = "ri",
      type = "both",
      n_points = NULL,
      mw_range = c(1, 1e15),
      normalize = TRUE
    ) |>
    recipes::prep()

  fresh <- mw_distribution_data(
    location,
    list(rep(1, length(location)), rep(1, length(location))),
    sample_id = c("tail", "head")
  )
  fresh$.excluded_regions <- list(
    tibble::tibble(start = 4, end = 10, purpose = "integration"),
    tibble::tibble(start = 0, end = 6, purpose = "integration")
  )
  result <- recipes::bake(prepped, new_data = fresh)
  tail_density <- result$ri[[1]]
  head_density <- result$ri[[2]]

  expect_equal(mw_curve_area(tail_density), 1, tolerance = 1e-10)
  expect_equal(mw_curve_area(head_density), 1, tolerance = 1e-10)
  expect_equal(
    mw_curve_area(tail_density[tail_density$location >= 4, ]),
    0,
    tolerance = 1e-10
  )
  expect_equal(
    mw_curve_area(head_density[head_density$location <= 6, ]),
    0,
    tolerance = 1e-10
  )

  full_support <- recipes::recipe(training) |>
    step_sec_mw_distribution(
      measures = "ri",
      type = "both",
      n_points = NULL,
      mw_range = c(1, 1e15),
      normalize = TRUE
    ) |>
    recipes::prep() |>
    recipes::bake(new_data = NULL)
  expect_equal(mw_curve_area(full_support$ri[[1]]), 1, tolerance = 1e-10)
  expect_equal(
    mw_curve_area(full_support$ri[[1]][
      full_support$ri[[1]]$location >= 10,
    ]),
    0,
    tolerance = 1e-10
  )
})

test_that("segmented distributions support descending irregular MW axes", {
  skip_if_not_installed("measure")

  location <- c(10, 7, 3, 0)
  mw <- c(100, 1000, 1e4, 1e6)
  data <- mw_distribution_data(location, rep(1, length(location)))
  data$mw <- measure::new_measure_list(list(
    measure::new_measure_tbl(location, mw)
  ))
  rec <- recipes::recipe(data) |>
    step_sec_exclude_regions(
      regions = tibble::tibble(start = 3.5, end = 7.5),
      purpose = "integration"
    ) |>
    step_sec_mw_distribution(
      measures = "ri",
      mw_col = "mw",
      mw_scale = "absolute",
      type = "both",
      n_points = 10,
      normalize = FALSE
    )

  result <- recipes::bake(recipes::prep(rec), new_data = NULL)
  cdf <- result$ri_cumulative[[1]]
  gap_left <- log10(stats::approx(location, mw, xout = 7.5)$y)
  gap_right <- log10(stats::approx(location, mw, xout = 3.5)$y)

  expect_equal(mw_curve_area(result$ri[[1]]), 6, tolerance = 1e-10)
  expect_equal(
    mw_curve_at(cdf, c(gap_left, mean(c(gap_left, gap_right)), gap_right)),
    rep(mw_curve_at(cdf, gap_left), 3),
    tolerance = 1e-10
  )
  expect_true(all(diff(cdf$location) > 0))
})

test_that("baking uses each fresh row's metadata after serialization", {
  skip_if_not_installed("measure")

  location <- 0:10
  training <- mw_distribution_data(location, rep(1, length(location)))
  training$.integration_window <- list(
    tibble::tibble(start = 0, end = 10)
  )
  training$.excluded_regions <- list(mw_empty_exclusions())
  prepped <- recipes::recipe(training) |>
    step_sec_mw_distribution(
      measures = "ri",
      type = "cumulative",
      n_points = 11
    ) |>
    recipes::prep()
  serialized <- withr::local_tempfile(fileext = ".rds")
  saveRDS(prepped, serialized)
  prepped <- readRDS(serialized)

  fresh <- mw_distribution_data(
    location,
    list(rep(1, length(location)), rep(1, length(location))),
    sample_id = c("later", "earlier")
  )
  fresh$.integration_window <- list(
    tibble::tibble(start = 0, end = 8),
    tibble::tibble(start = 2, end = 10)
  )
  fresh$.excluded_regions <- list(
    tibble::tibble(start = 2, end = 3, purpose = "integration"),
    tibble::tibble(start = 6, end = 7, purpose = "integration")
  )

  result <- recipes::bake(prepped, new_data = fresh)

  expect_equal(result$sample_id, c("later", "earlier"))
  expect_equal(range(result$ri[[1]]$location), c(0, 8))
  expect_equal(range(result$ri[[2]]$location), c(2, 10))
  expect_equal(
    mw_curve_at(result$ri[[1]], c(2, 3)),
    rep(2 / 7, 2),
    tolerance = 1e-10
  )
  expect_equal(
    mw_curve_at(result$ri[[2]], c(6, 7)),
    rep(4 / 7, 2),
    tolerance = 1e-10
  )
  expect_gt(
    mw_curve_at(result$ri[[2]], 2.5),
    mw_curve_at(result$ri[[2]], 2)
  )

  override_recipe <- recipes::recipe(training) |>
    step_sec_mw_distribution(
      measures = "ri",
      type = "cumulative",
      n_points = 11,
      integration_range = c(1, 5)
    ) |>
    recipes::prep()
  fresh$.excluded_regions <- list(
    mw_empty_exclusions(),
    mw_empty_exclusions()
  )
  overridden <- recipes::bake(override_recipe, new_data = fresh)
  expect_equal(range(overridden$ri[[1]]$location), c(1, 5))
  expect_equal(range(overridden$ri[[2]]$location), c(1, 5))
})

test_that("baseline-only exclusions are ignored and malformed metadata fails", {
  skip_if_not_installed("measure")

  location <- 0:10
  data <- mw_distribution_data(location, rep(1, length(location)))
  baseline_recipe <- recipes::recipe(data) |>
    step_sec_exclude_regions(
      regions = tibble::tibble(start = 2, end = 3),
      purpose = "baseline"
    ) |>
    step_sec_mw_distribution(
      measures = "ri",
      type = "cumulative",
      n_points = 11
    )
  baseline <- recipes::bake(
    recipes::prep(baseline_recipe),
    new_data = NULL
  )
  expect_gt(
    mw_curve_at(baseline$ri[[1]], 3),
    mw_curve_at(baseline$ri[[1]], 2)
  )

  metadata_training <- data
  metadata_training$.excluded_regions <- list(mw_empty_exclusions())
  prepped <- recipes::recipe(metadata_training) |>
    step_sec_mw_distribution(
      measures = "ri",
      type = "cumulative",
      n_points = 11
    ) |>
    recipes::prep()
  malformed <- data
  malformed$.excluded_regions <- list(tibble::tibble(
    start = "not numeric",
    end = 3,
    purpose = "integration"
  ))
  expect_snapshot(error = TRUE, recipes::bake(prepped, new_data = malformed))
})

test_that("missing calibration is tolerated only inside excluded support", {
  skip_if_not_installed("measure")

  location <- 1:7
  signal_data <- mw_distribution_data(location, rep(1, length(location)))
  signal_data$mw <- measure::new_measure_list(list(
    measure::new_measure_tbl(
      location,
      c(100, 80, 60, NA_real_, 20, 10, 5)
    )
  ))
  rec <- recipes::recipe(signal_data) |>
    step_sec_exclude_regions(
      regions = tibble::tibble(start = 3, end = 5),
      purpose = "integration"
    ) |>
    step_sec_mw_distribution(
      measures = "ri",
      mw_col = "mw",
      mw_scale = "absolute",
      type = "both",
      n_points = 10,
      normalize = FALSE
    )

  result <- recipes::bake(recipes::prep(rec), new_data = NULL)
  expect_equal(mw_curve_area(result$ri[[1]]), 4, tolerance = 1e-10)

  expect_snapshot(
    error = TRUE,
    recipes::recipe(signal_data) |>
      step_sec_mw_distribution(
        measures = "ri",
        mw_col = "mw",
        mw_scale = "absolute",
        type = "both",
        n_points = 10,
        normalize = FALSE
      ) |>
      recipes::prep() |>
      recipes::bake(new_data = NULL)
  )
})
