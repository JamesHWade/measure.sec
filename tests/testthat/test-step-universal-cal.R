# ==============================================================================
# Tests for step_sec_universal_cal
# ==============================================================================

test_that("step_sec_universal_cal requires Mark-Houwink parameters", {
  skip_if_not_installed("measure")

  log_mw <- seq(3, 6, by = 0.05)
  conc <- dnorm(log_mw, mean = 4.5, sd = 0.5)

  test_data <- tibble::tibble(sample_id = "test")
  test_data$ri <- measure::new_measure_list(
    list(measure::new_measure_tbl(location = log_mw, value = conc))
  )

  # Missing K_sample and a_sample

  expect_error(
    recipes::recipe(~., data = test_data) |>
      step_sec_universal_cal(
        calibration = data.frame(x = 1:10, log_eta_m = 10:1)
      ),
    "K_sample.*a_sample"
  )
})

test_that("step_sec_universal_cal validates K_sample", {
  skip_if_not_installed("measure")

  log_mw <- seq(3, 6, by = 0.05)
  conc <- dnorm(log_mw, mean = 4.5, sd = 0.5)

  test_data <- tibble::tibble(sample_id = "test")
  test_data$ri <- measure::new_measure_list(
    list(measure::new_measure_tbl(location = log_mw, value = conc))
  )

  # Negative K_sample
  expect_error(
    recipes::recipe(~., data = test_data) |>
      step_sec_universal_cal(
        calibration = data.frame(x = 1:10, log_eta_m = 10:1),
        K_sample = -0.0001,
        a_sample = 0.7
      ),
    "K_sample.*positive"
  )

  # Zero K_sample
  expect_error(
    recipes::recipe(~., data = test_data) |>
      step_sec_universal_cal(
        calibration = data.frame(x = 1:10, log_eta_m = 10:1),
        K_sample = 0,
        a_sample = 0.7
      ),
    "K_sample.*positive"
  )
})

test_that("step_sec_universal_cal validates a_sample range", {
  skip_if_not_installed("measure")

  log_mw <- seq(3, 6, by = 0.05)
  conc <- dnorm(log_mw, mean = 4.5, sd = 0.5)

  test_data <- tibble::tibble(sample_id = "test")
  test_data$ri <- measure::new_measure_list(
    list(measure::new_measure_tbl(location = log_mw, value = conc))
  )

  # a_sample > 1
  expect_error(
    recipes::recipe(~., data = test_data) |>
      step_sec_universal_cal(
        calibration = data.frame(x = 1:10, log_eta_m = 10:1),
        K_sample = 0.0001,
        a_sample = 1.5
      ),
    "a_sample.*between 0 and 1"
  )

  # a_sample <= 0
  expect_error(
    recipes::recipe(~., data = test_data) |>
      step_sec_universal_cal(
        calibration = data.frame(x = 1:10, log_eta_m = 10:1),
        K_sample = 0.0001,
        a_sample = 0
      ),
    "a_sample.*between 0 and 1"
  )
})

test_that("step_sec_universal_cal requires calibration", {
  skip_if_not_installed("measure")

  log_mw <- seq(3, 6, by = 0.05)
  conc <- dnorm(log_mw, mean = 4.5, sd = 0.5)

  test_data <- tibble::tibble(sample_id = "test")
  test_data$ri <- measure::new_measure_list(
    list(measure::new_measure_tbl(location = log_mw, value = conc))
  )

  # No calibration provided
  expect_error(
    recipes::recipe(~., data = test_data) |>
      step_sec_universal_cal(
        K_sample = 0.0001,
        a_sample = 0.7
      ),
    "calibration"
  )
})

test_that("step_sec_universal_cal print method works", {
  skip_if_not_installed("measure")

  log_mw <- seq(3, 6, by = 0.05)
  conc <- dnorm(log_mw, mean = 4.5, sd = 0.5)

  test_data <- tibble::tibble(sample_id = "test")
  test_data$ri <- measure::new_measure_list(
    list(measure::new_measure_tbl(location = log_mw, value = conc))
  )

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_universal_cal(
      calibration = data.frame(x = 1:10, log_eta_m = 10:1),
      K_sample = 0.000128,
      a_sample = 0.690
    )

  expect_output(print(rec), "universal")
})

test_that("step_sec_universal_cal tidy method returns expected structure", {
  skip_if_not_installed("measure")

  log_mw <- seq(3, 6, by = 0.05)
  conc <- dnorm(log_mw, mean = 4.5, sd = 0.5)

  test_data <- tibble::tibble(sample_id = "test")
  test_data$ri <- measure::new_measure_list(
    list(measure::new_measure_tbl(location = log_mw, value = conc))
  )

  rec <- recipes::recipe(~., data = test_data) |>
    step_sec_universal_cal(
      calibration = data.frame(x = 1:10, log_eta_m = 10:1),
      K_sample = 0.000128,
      a_sample = 0.690
    )

  tidy_result <- recipes::tidy(rec, number = 1)

  expect_s3_class(tidy_result, "tbl_df")
  expect_true("K_sample" %in% names(tidy_result))
  expect_true("a_sample" %in% names(tidy_result))
  expect_true("K_standard" %in% names(tidy_result))
  expect_true("a_standard" %in% names(tidy_result))
})
