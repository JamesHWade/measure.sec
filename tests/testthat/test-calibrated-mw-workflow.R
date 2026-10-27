test_that("calibrated SEC averages use the per-slice molecular weights", {
  location <- c(1, 1.5, 2, 2.5, 3)
  signal <- c(1, 2, 3, 2, 1)
  log_mw <- 6 - 0.4 * location + 0.03 * location^2
  standards <- data.frame(location = location, log_mw = log_mw)
  data <- tibble::tibble(sample_id = "sample")
  data$ri <- measure::new_measure_list(list(
    measure::new_measure_tbl(location, signal)
  ))

  rec <- recipes::recipe(data) |>
    step_sec_conventional_cal(
      measures = "ri",
      standards = standards,
      fit_type = "quadratic",
      extrapolation = "none"
    ) |>
    step_sec_mw_averages(measures = "ri")
  result <- recipes::bake(recipes::prep(rec), new_data = NULL)

  # Trapezoidal sample contributions on the elution axis.
  weights <- signal * c(0.25, 0.5, 0.5, 0.5, 0.25)
  mw <- 10^log_mw
  expect_equal(result$mw_mw, sum(weights * mw) / sum(weights), tolerance = 1e-8)
  expect_equal(result$mw_mn, sum(weights) / sum(weights / mw), tolerance = 1e-8)
  expect_equal(
    result$mw_mz,
    sum(weights * mw^2) / sum(weights * mw),
    tolerance = 1e-8
  )
})

test_that("SEC averages can pair an absolute MW curve with detector signal", {
  data <- tibble::tibble(sample_id = "sample")
  data$ri <- measure::new_measure_list(list(
    measure::new_measure_tbl(c(1, 2, 4), c(1, 2, 1))
  ))
  data$mw_mals <- measure::new_measure_list(list(
    measure::new_measure_tbl(c(1, 2, 4), c(1e5, 2e5, 4e5))
  ))
  rec <- recipes::recipe(data) |>
    step_sec_mw_averages(
      measures = "ri",
      mw_col = "mw_mals",
      mw_scale = "absolute"
    )
  result <- recipes::bake(recipes::prep(rec), new_data = NULL)
  weights <- c(0.5, 3, 1)
  expect_equal(result$mw_mw, sum(weights * c(1e5, 2e5, 4e5)) / sum(weights))
})

test_that("SEC averages reject a misaligned calibrated curve", {
  data <- tibble::tibble(sample_id = "sample")
  data$ri <- measure::new_measure_list(list(
    measure::new_measure_tbl(c(1, 2, 3), c(1, 2, 1))
  ))
  data$mw <- measure::new_measure_list(list(
    measure::new_measure_tbl(c(1, 2.1, 3), c(5, 4, 3))
  ))
  rec <- recipes::recipe(data) |>
    step_sec_mw_averages(measures = "ri", mw_col = "mw", mw_scale = "log10")
  expect_error(recipes::prep(rec), "align|location")
})

test_that("SEC averages do not silently discard uncalibrated signal", {
  data <- tibble::tibble(sample_id = "sample")
  data$ri <- measure::new_measure_list(list(
    measure::new_measure_tbl(c(1, 2, 3), c(1, 2, 1))
  ))
  data$mw <- measure::new_measure_list(list(
    measure::new_measure_tbl(c(1, 2, 3), c(5, NA_real_, 3))
  ))
  rec <- recipes::recipe(data) |>
    step_sec_mw_averages(measures = "ri", mw_col = "mw", mw_scale = "log10")
  expect_error(recipes::prep(rec), "calibrat|molecular weight")
})

test_that("SEC fractions use the calibrated MW and elution-axis weights", {
  location <- c(1, 1.5, 2, 2.5, 3)
  signal <- c(1, 2, 3, 2, 1)
  log_mw <- 6 - 0.4 * location + 0.03 * location^2
  data <- tibble::tibble(sample_id = "sample")
  data$ri <- measure::new_measure_list(list(measure::new_measure_tbl(
    location,
    signal
  )))
  rec <- recipes::recipe(data) |>
    step_sec_conventional_cal(
      measures = "ri",
      standards = data.frame(location, log_mw),
      fit_type = "quadratic",
      extrapolation = "none"
    ) |>
    step_sec_mw_fractions(measures = "ri", cutoffs = 2e5)
  result <- recipes::bake(recipes::prep(rec), new_data = NULL)
  weights <- signal * c(0.25, 0.5, 0.5, 0.5, 0.25)
  expected <- sum(weights[10^log_mw < 2e5]) / sum(weights)
  expect_equal(result[[paste0("frac_below_", 2e5)]], expected)
})

test_that("nonlinear calibration preserves mass in SEC distributions", {
  location <- seq(1, 3, length.out = 201)
  log_mw <- 6 - 0.4 * location + 0.03 * location^2
  data <- tibble::tibble(sample_id = "sample")
  data$ri <- measure::new_measure_list(list(
    measure::new_measure_tbl(location, rep(1, length(location)))
  ))
  standards <- data.frame(
    location = location[c(1, 51, 101, 151, 201)],
    log_mw = log_mw[c(1, 51, 101, 151, 201)]
  )
  rec <- recipes::recipe(data) |>
    step_sec_conventional_cal(
      measures = "ri",
      standards = standards,
      fit_type = "quadratic",
      extrapolation = "none"
    ) |>
    step_sec_mw_distribution(measures = "ri", type = "both", n_points = 201)
  result <- recipes::bake(recipes::prep(rec), new_data = NULL)
  density <- result$ri_differential[[1]]
  cum <- result$ri_cumulative[[1]]
  area <- sum(
    diff(density$location) *
      (head(density$value, -1) + tail(density$value, -1)) /
      2
  )
  expect_equal(area, 1, tolerance = 1e-8)
  expect_equal(cum$value[1], 0, tolerance = 1e-8)
  expect_equal(tail(cum$value, 1), 1, tolerance = 1e-8)
  expect_equal(
    stats::approx(cum$location, cum$value, xout = log_mw[101])$y,
    0.5,
    tolerance = 1e-3
  )
})
