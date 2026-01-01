# ==============================================================================
# Tests for QC functions
# ==============================================================================

test_that("measure_sec_resolution calculates USP resolution", {
  # Two peaks with known separation
  Rs <- measure_sec_resolution(
    retention_1 = 8.0,
    retention_2 = 10.0,
    width_1 = 0.5,
    width_2 = 0.5,
    method = "usp"
  )

  # Rs = 2 * (10 - 8) / (0.5 + 0.5) = 2 * 2 / 1 = 4
 expect_equal(Rs, 4.0)

  # Test with EP method
  Rs_ep <- measure_sec_resolution(
    retention_1 = 8.0,
    retention_2 = 10.0,
    width_1 = 0.5,
    width_2 = 0.5,
    method = "ep"
  )

  # Rs = 1.18 * 2 / 1 = 2.36
  expect_equal(Rs_ep, 2.36, tolerance = 0.01)
})

test_that("measure_sec_resolution handles reversed peak order", {
  Rs1 <- measure_sec_resolution(8.0, 10.0, 0.5, 0.5)
  Rs2 <- measure_sec_resolution(10.0, 8.0, 0.5, 0.5)

  expect_equal(Rs1, Rs2)
})

test_that("measure_sec_resolution rejects invalid inputs", {
  expect_error(measure_sec_resolution(8.0, 10.0, -0.5, 0.5), "positive")
  expect_error(measure_sec_resolution(8.0, 10.0, 0.5, 0), "positive")
})

test_that("measure_sec_plate_count calculates N correctly", {
  # Known example: retention = 10 min, width at half height = 0.5 min
  # N = 5.54 * (10/0.5)^2 = 5.54 * 400 = 2216
  N <- measure_sec_plate_count(
    retention = 10.0,
    width = 0.5,
    width_type = "half_height"
  )

  expect_equal(N, 2216, tolerance = 1)

  # Test baseline width method
  # N = 16 * (10/0.5)^2 = 16 * 400 = 6400
  N_baseline <- measure_sec_plate_count(
    retention = 10.0,
    width = 0.5,
    width_type = "baseline"
  )

  expect_equal(N_baseline, 6400)
})

test_that("measure_sec_plate_count calculates effective plates with dead time", {
  # With dead time = 2 min, effective retention = 8 min
  # N_eff = 5.54 * (8/0.5)^2 = 5.54 * 256 = 1418.24
  N_eff <- measure_sec_plate_count(
    retention = 10.0,
    width = 0.5,
    dead_time = 2.0
  )

  expect_equal(N_eff, 1418.24, tolerance = 1)
})

test_that("measure_sec_plate_count rejects invalid inputs", {
  expect_error(measure_sec_plate_count(-10, 0.5), "positive")
  expect_error(measure_sec_plate_count(10, 0), "positive")
  expect_error(measure_sec_plate_count(10, 0.5, dead_time = 15), "less than")
})

test_that("measure_sec_asymmetry calculates USP tailing factor", {
  # Leading = 0.2, Tailing = 0.3
  # USP Tf = (0.2 + 0.3) / (2 * 0.2) = 0.5 / 0.4 = 1.25
  Tf <- measure_sec_asymmetry(
    leading = 0.2,
    tailing = 0.3,
    method = "usp"
  )

  expect_equal(Tf, 1.25)
})

test_that("measure_sec_asymmetry calculates EP asymmetry factor", {
  # EP As = tailing / leading = 0.3 / 0.2 = 1.5
  As <- measure_sec_asymmetry(
    leading = 0.2,
    tailing = 0.3,
    method = "ep"
  )

  expect_equal(As, 1.5)
})

test_that("measure_sec_recovery calculates percent correctly", {
  recovery <- measure_sec_recovery(
    detected_mass = 0.195,
    injected_mass = 0.200
  )

  expect_equal(recovery, 97.5)
})

test_that("measure_sec_recovery rejects invalid inputs", {
  expect_error(measure_sec_recovery(-0.1, 0.2), "non-negative")
  expect_error(measure_sec_recovery(0.2, 0), "positive")
})

test_that("measure_sec_suitability performs comprehensive testing", {
  peaks <- data.frame(
    name = c("dimer", "monomer"),
    retention = c(8.5, 10.0),
    width = c(0.3, 0.35),
    area = c(5, 95)
  )

  sst <- measure_sec_suitability(
    peaks = peaks,
    reference_peaks = c("dimer", "monomer")
  )

  expect_s3_class(sst, "sec_suitability")
  expect_true("results" %in% names(sst))
  expect_true("passed" %in% names(sst))
  expect_true("resolution" %in% names(sst$results))
  expect_true("plate_count" %in% names(sst$results))
})

test_that("measure_sec_suitability print method works", {
  peaks <- data.frame(
    name = c("dimer", "monomer"),
    retention = c(8.5, 10.0),
    width = c(0.3, 0.35),
    area = c(5, 95)
  )

  sst <- measure_sec_suitability(peaks = peaks)

  expect_output(print(sst), "System Suitability")
  expect_output(print(sst), "plate")
})

test_that("measure_sec_suitability requires valid peaks data", {
  expect_error(
    measure_sec_suitability(peaks = "not a data frame"),
    "data frame"
  )

  expect_error(
    measure_sec_suitability(peaks = data.frame(x = 1)),
    "Missing"
  )
})
