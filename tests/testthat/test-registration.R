# Test that measure.sec registers correctly with measure

test_that("measure.sec registers as a pack", {
  skip_if_not_installed("measure")

  # The pack should be registered on package load
  packs <- measure::measure_packs()
  expect_true("measure.sec" %in% packs$name)

  # Check technique is correct
  sec_pack <- packs[packs$name == "measure.sec", ]
  expect_equal(sec_pack$technique, "SEC/GPC")
})

test_that("SEC steps are registered", {
  skip_if_not_installed("measure")

  steps <- measure::measure_steps(techniques = "SEC/GPC")
  expect_gt(nrow(steps), 0)

  # Check specific steps are registered
  step_names <- steps$step_name
  expect_true("step_sec_baseline" %in% step_names)
  expect_true("step_sec_mw_averages" %in% step_names)
  expect_true("step_sec_mw_fractions" %in% step_names)
  expect_true("step_sec_mw_distribution" %in% step_names)
})

test_that("SEC steps have correct categories",
 {
  skip_if_not_installed("measure")

  steps <- measure::measure_steps(techniques = "SEC/GPC")

  baseline_step <- steps[steps$step_name == "step_sec_baseline", ]
  expect_equal(baseline_step$category, "baseline")

  mw_avg_step <- steps[steps$step_name == "step_sec_mw_averages", ]
  expect_equal(mw_avg_step$category, "calculation")
})
