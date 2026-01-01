# ==============================================================================
# Tests for detector processing steps
# ==============================================================================

test_that("step_sec_detector_delay corrects inter-detector delays", {
  skip_if_not_installed("measure")

  # Create test data with offset signals
  set.seed(42)
  time <- seq(5, 15, by = 0.1)

  # RI signal centered at t=10
  ri_signal <- dnorm(time, mean = 10, sd = 0.5)

  # UV signal offset by 0.5 min (sees sample earlier)
  uv_signal <- dnorm(time, mean = 9.5, sd = 0.5)

  test_data <- tibble::tibble(sample_id = "test")

  test_data$ri <- measure::new_measure_list(
    list(measure::new_measure_tbl(location = time, value = ri_signal))
  )
  test_data$uv <- measure::new_measure_list(
    list(measure::new_measure_tbl(location = time, value = uv_signal))
  )

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_detector_delay(
      reference = "ri",
      delay_volumes = c(uv = -0.5),
      flow_rate = 1.0
    )

  prepped <- recipes::prep(rec)
  result <- recipes::bake(prepped, new_data = NULL)

  expect_s3_class(result, "tbl_df")
  expect_true("ri" %in% names(result))
  expect_true("uv" %in% names(result))
})

test_that("step_sec_ri applies dn/dc normalization", {
  skip_if_not_installed("measure")

  # Create simple test data
  time <- seq(5, 15, by = 0.1)
  ri_signal <- dnorm(time, mean = 10, sd = 0.5)

  test_data <- tibble::tibble(sample_id = "test")

  test_data$ri <- measure::new_measure_list(
    list(measure::new_measure_tbl(location = time, value = ri_signal))
  )

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_ri(measures = "ri", dn_dc = 0.185)

  prepped <- recipes::prep(rec)
  result <- recipes::bake(prepped, new_data = NULL)

  expect_s3_class(result, "tbl_df")
  expect_true("ri" %in% names(result))

  # Signal should be divided by dn/dc
  original_max <- max(ri_signal)
  processed_max <- max(result$ri[[1]]$value)
  expect_equal(processed_max, original_max / 0.185, tolerance = 0.01)
})

test_that("step_sec_uv applies extinction coefficient", {
  skip_if_not_installed("measure")

  time <- seq(5, 15, by = 0.1)
  uv_signal <- dnorm(time, mean = 10, sd = 0.5)

  test_data <- tibble::tibble(sample_id = "test")

  test_data$uv <- measure::new_measure_list(
    list(measure::new_measure_tbl(location = time, value = uv_signal))
  )

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_uv(measures = "uv", extinction_coef = 1.5, path_length = 1.0)

  prepped <- recipes::prep(rec)
  result <- recipes::bake(prepped, new_data = NULL)

  expect_s3_class(result, "tbl_df")

  # Signal should be divided by extinction coef * path length
  original_max <- max(uv_signal)
  processed_max <- max(result$uv[[1]]$value)
  expect_equal(processed_max, original_max / 1.5, tolerance = 0.01)
})

test_that("step_sec_concentration converts to concentration", {
  skip_if_not_installed("measure")

  time <- seq(5, 15, by = 0.1)
  ri_signal <- dnorm(time, mean = 10, sd = 0.5)

  test_data <- tibble::tibble(sample_id = "test")

  test_data$ri <- measure::new_measure_list(
    list(measure::new_measure_tbl(location = time, value = ri_signal))
  )

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_concentration(
      measures = "ri",
      detector = "ri",
      injection_mass = 0.2,
      flow_rate = 1.0
    )

  prepped <- recipes::prep(rec)
  result <- recipes::bake(prepped, new_data = NULL)

  expect_s3_class(result, "tbl_df")
  expect_true("ri" %in% names(result))
})

test_that("step_sec_mals requires dn_dc parameter", {
  skip_if_not_installed("measure")

  time <- seq(5, 15, by = 0.1)
  test_data <- tibble::tibble(sample_id = "test")

  test_data$mals <- measure::new_measure_list(
    list(measure::new_measure_tbl(
      location = time,
      value = dnorm(time, 10, 0.5)
    ))
  )
  test_data$ri <- measure::new_measure_list(
    list(measure::new_measure_tbl(
      location = time,
      value = dnorm(time, 10, 0.5)
    ))
  )

  expect_error(
    recipes::recipe(~., data = test_data) |>
      step_sec_mals(mals_col = "mals"),
    "dn_dc"
  )
})

test_that("step_sec_viscometer processes DP signal", {
  skip_if_not_installed("measure")

  time <- seq(5, 15, by = 0.1)
  dp_signal <- dnorm(time, mean = 10, sd = 0.5) * 100 # DP in arbitrary units

  test_data <- tibble::tibble(sample_id = "test")

  test_data$dp <- measure::new_measure_list(
    list(measure::new_measure_tbl(location = time, value = dp_signal))
  )

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_viscometer(dp_col = "dp", viscometer_constant = 0.01)

  prepped <- recipes::prep(rec)
  result <- recipes::bake(prepped, new_data = NULL)

  expect_s3_class(result, "tbl_df")
  expect_true("specific_visc" %in% names(result))
})

test_that("step_sec_intrinsic_visc calculates from viscosity and concentration", {
  skip_if_not_installed("measure")

  time <- seq(5, 15, by = 0.1)
  specific_visc <- dnorm(time, mean = 10, sd = 0.5) * 0.1
  concentration <- dnorm(time, mean = 10, sd = 0.5)

  test_data <- tibble::tibble(sample_id = "test")

  # Create measure objects
  test_data$visc <- measure::new_measure_list(
    list(measure::new_measure_tbl(location = time, value = specific_visc))
  )
  test_data$conc <- measure::new_measure_list(
    list(measure::new_measure_tbl(location = time, value = concentration))
  )

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_intrinsic_visc(
      specific_visc_col = "visc",
      concentration_col = "conc"
    )

  prepped <- recipes::prep(rec)
  result <- recipes::bake(prepped, new_data = NULL)

  expect_s3_class(result, "tbl_df")
  expect_true("intrinsic_visc" %in% names(result))
})
