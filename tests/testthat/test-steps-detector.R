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
		list(
			measure::new_measure_tbl(
				location = time,
				value = dnorm(time, 10, 0.5)
			)
		)
	)
	test_data$ri <- measure::new_measure_list(
		list(
			measure::new_measure_tbl(
				location = time,
				value = dnorm(time, 10, 0.5)
			)
		)
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

test_that("step_sec_dad applies extinction coefficients and ratios", {
	skip_if_not_installed("measure")

	time <- seq(1, 10, by = 0.5)
	uv_254 <- dnorm(time, mean = 5, sd = 1)
	uv_280 <- dnorm(time, mean = 5.5, sd = 1)

	test_data <- tibble::tibble(sample_id = "test")
	test_data$uv_254 <- measure::new_measure_list(
		list(measure::new_measure_tbl(location = time, value = uv_254))
	)
	test_data$uv_280 <- measure::new_measure_list(
		list(measure::new_measure_tbl(location = time, value = uv_280))
	)

	rec <- recipes::recipe(~., data = test_data) |>
		step_sec_dad(
			measures = c("uv_254", "uv_280"),
			wavelengths = c(254, 280),
			extinction_coefs = c(`254` = 2, `280` = 1),
			reference_wavelength = 280,
			output_prefix = "uv"
		)

	prepped <- recipes::prep(rec)
	result <- recipes::bake(prepped, new_data = NULL)

	expect_true("uv_254" %in% names(result))
	expect_true("uv_280" %in% names(result))
	expect_true("uv_254_to_280" %in% names(result))

	processed_max <- max(result$uv_254[[1]]$value)
	expect_equal(processed_max, max(uv_254) / 2, tolerance = 0.01)
})

test_that("step_sec_lals computes MW and validates angle", {
	skip_if_not_installed("measure")

	time <- seq(1, 10, by = 0.5)
	lals_signal <- dnorm(time, mean = 5, sd = 1)
	ri_signal <- dnorm(time, mean = 5, sd = 1)

	test_data <- tibble::tibble(sample_id = "test")
	test_data$lals <- measure::new_measure_list(
		list(measure::new_measure_tbl(location = time, value = lals_signal))
	)
	test_data$ri <- measure::new_measure_list(
		list(measure::new_measure_tbl(location = time, value = ri_signal))
	)

	rec <- recipes::recipe(~., data = test_data) |>
		step_sec_lals(
			measures = "lals",
			concentration_col = "ri",
			dn_dc = 0.185,
			calibration_constant = 1.0
		)

	prepped <- recipes::prep(rec)
	result <- recipes::bake(prepped, new_data = NULL)

	expect_true("mw_lals" %in% names(result))
	expect_true(!all(is.na(result$mw_lals[[1]]$value)))

	expect_error(
		recipes::recipe(~., data = test_data) |>
			step_sec_lals(
				measures = "lals",
				concentration_col = "ri",
				dn_dc = 0.185,
				angle = 25
			),
		"angle"
	)
})

test_that("step_sec_rals computes MW", {
	skip_if_not_installed("measure")

	time <- seq(1, 10, by = 0.5)
	rals_signal <- dnorm(time, mean = 5, sd = 1)
	ri_signal <- dnorm(time, mean = 5, sd = 1)

	test_data <- tibble::tibble(sample_id = "test")
	test_data$rals <- measure::new_measure_list(
		list(measure::new_measure_tbl(location = time, value = rals_signal))
	)
	test_data$ri <- measure::new_measure_list(
		list(measure::new_measure_tbl(location = time, value = ri_signal))
	)

	rec <- recipes::recipe(~., data = test_data) |>
		step_sec_rals(
			measures = "rals",
			concentration_col = "ri",
			dn_dc = 0.185,
			calibration_constant = 1.0
		)

	prepped <- recipes::prep(rec)
	result <- recipes::bake(prepped, new_data = NULL)

	expect_true("mw_rals" %in% names(result))
	expect_true(!all(is.na(result$mw_rals[[1]]$value)))
})

test_that("step_sec_dls estimates diffusion and Rh", {
	skip_if_not_installed("measure")

	tau <- seq(1e-6, 1e-3, length.out = 50)
	laser_wavelength <- 633
	angle <- 90
	solvent_ri <- 1.333
	diffusion <- 1e-10

	q_val <- (4 * pi * solvent_ri / (laser_wavelength * 1e-9)) * sin(pi / 4)
	gamma <- diffusion * q_val^2
	beta <- 0.8
	g2 <- 1 + beta * exp(-2 * gamma * tau)

	test_data <- tibble::tibble(sample_id = "test")
	test_data$dls <- measure::new_measure_list(
		list(measure::new_measure_tbl(location = tau, value = g2))
	)
	test_data$rg <- measure::new_measure_list(
		list(measure::new_measure_tbl(location = tau, value = rep(10, length(tau))))
	)

	rec <- recipes::recipe(~., data = test_data) |>
		step_sec_dls(
			measures = "dls",
			temperature = 25,
			viscosity = 1.0,
			laser_wavelength = laser_wavelength,
			angle = angle,
			solvent_ri = solvent_ri
		)

	prepped <- recipes::prep(rec)
	result <- recipes::bake(prepped, new_data = NULL)

	expect_true("rh" %in% names(result))
	expect_true("diffusion_coef" %in% names(result))
	expect_true("rg_rh" %in% names(result))

	rh_vals <- result$rh[[1]]$value
	expect_true(all(is.finite(rh_vals)))

	k_b <- 1.380649e-23
	expected_rh_m <- k_b * (25 + 273.15) / (6 * pi * 1e-3 * diffusion)
	expect_equal(mean(rh_vals), expected_rh_m * 1e9, tolerance = 0.2)
})

# ==============================================================================
# Parameter validation tests
# ==============================================================================

test_that("step_sec_lals requires dn_dc or optical_constant", {
	skip_if_not_installed("measure")

	time <- seq(1, 10, by = 0.5)
	test_data <- tibble::tibble(sample_id = "test")
	test_data$lals <- measure::new_measure_list(
		list(measure::new_measure_tbl(location = time, value = dnorm(time, 5, 1)))
	)
	test_data$ri <- measure::new_measure_list(
		list(measure::new_measure_tbl(location = time, value = dnorm(time, 5, 1)))
	)

	expect_error(
		recipes::recipe(~., data = test_data) |>
			step_sec_lals(measures = "lals", concentration_col = "ri"),
		"dn_dc|optical_constant"
	)
})

test_that("step_sec_rals requires dn_dc or optical_constant", {
	skip_if_not_installed("measure")

	time <- seq(1, 10, by = 0.5)
	test_data <- tibble::tibble(sample_id = "test")
	test_data$rals <- measure::new_measure_list(
		list(measure::new_measure_tbl(location = time, value = dnorm(time, 5, 1)))
	)
	test_data$ri <- measure::new_measure_list(
		list(measure::new_measure_tbl(location = time, value = dnorm(time, 5, 1)))
	)

	expect_error(
		recipes::recipe(~., data = test_data) |>
			step_sec_rals(measures = "rals", concentration_col = "ri"),
		"dn_dc|optical_constant"
	)
})

test_that("step_sec_dls validates temperature above absolute zero", {
	skip_if_not_installed("measure")

	time <- seq(1, 10, by = 0.5)
	test_data <- tibble::tibble(sample_id = "test")
	test_data$dls <- measure::new_measure_list(
		list(measure::new_measure_tbl(location = time, value = dnorm(time, 5, 1)))
	)

	expect_error(
		recipes::recipe(~., data = test_data) |>
			step_sec_dls(measures = "dls", temperature = -300),
		"absolute zero"
	)
})

test_that("step_sec_dls validates positive viscosity", {
	skip_if_not_installed("measure")

	time <- seq(1, 10, by = 0.5)
	test_data <- tibble::tibble(sample_id = "test")
	test_data$dls <- measure::new_measure_list(
		list(measure::new_measure_tbl(location = time, value = dnorm(time, 5, 1)))
	)

	expect_error(
		recipes::recipe(~., data = test_data) |>
			step_sec_dls(measures = "dls", temperature = 25, viscosity = -1),
		"positive"
	)
})

test_that("step_sec_dls validates angle range", {
	skip_if_not_installed("measure")

	time <- seq(1, 10, by = 0.5)
	test_data <- tibble::tibble(sample_id = "test")
	test_data$dls <- measure::new_measure_list(
		list(measure::new_measure_tbl(location = time, value = dnorm(time, 5, 1)))
	)

	expect_error(
		recipes::recipe(~., data = test_data) |>
			step_sec_dls(measures = "dls", temperature = 25, angle = 200),
		"angle"
	)
})

test_that("step_sec_dls validates output_cols is character", {
	skip_if_not_installed("measure")

	time <- seq(1, 10, by = 0.5)
	test_data <- tibble::tibble(sample_id = "test")
	test_data$dls <- measure::new_measure_list(
		list(measure::new_measure_tbl(location = time, value = dnorm(time, 5, 1)))
	)

	expect_error(
		recipes::recipe(~., data = test_data) |>
			step_sec_dls(measures = "dls", temperature = 25, output_cols = c(1, 2)),
		"output_cols"
	)
})

test_that("step_sec_dls validates output_cols length", {
	skip_if_not_installed("measure")

	time <- seq(1, 10, by = 0.5)
	test_data <- tibble::tibble(sample_id = "test")
	test_data$dls <- measure::new_measure_list(
		list(measure::new_measure_tbl(location = time, value = dnorm(time, 5, 1)))
	)

	expect_error(
		recipes::recipe(~., data = test_data) |>
			step_sec_dls(measures = "dls", temperature = 25, output_cols = "rh"),
		"output_cols"
	)
})

test_that("step_sec_dad validates positive path_length", {
	skip_if_not_installed("measure")

	time <- seq(1, 10, by = 0.5)
	test_data <- tibble::tibble(sample_id = "test")
	test_data$uv_254 <- measure::new_measure_list(
		list(measure::new_measure_tbl(location = time, value = dnorm(time, 5, 1)))
	)

	expect_error(
		recipes::recipe(~., data = test_data) |>
			step_sec_dad(measures = "uv_254", wavelengths = 254, path_length = -1),
		"path_length"
	)
})

# ==============================================================================
# tidy() method tests
# ==============================================================================

test_that("tidy.step_sec_dad returns expected structure", {
	skip_if_not_installed("measure")

	time <- seq(1, 10, by = 0.5)
	test_data <- tibble::tibble(sample_id = "test")
	test_data$uv_254 <- measure::new_measure_list(
		list(measure::new_measure_tbl(location = time, value = dnorm(time, 5, 1)))
	)
	test_data$uv_280 <- measure::new_measure_list(
		list(measure::new_measure_tbl(location = time, value = dnorm(time, 5, 1)))
	)

	rec <- recipes::recipe(~., data = test_data) |>
		step_sec_dad(
			measures = c("uv_254", "uv_280"),
			wavelengths = c(254, 280),
			reference_wavelength = 280
		)

	prepped <- recipes::prep(rec)
	tidy_result <- recipes::tidy(prepped, number = 1)

	expect_s3_class(tidy_result, "tbl_df")
	expect_true("measures" %in% names(tidy_result))
	expect_true("wavelengths" %in% names(tidy_result))
	expect_true("reference_wavelength" %in% names(tidy_result))
	expect_true("id" %in% names(tidy_result))
})

test_that("tidy.step_sec_lals returns expected structure", {
	skip_if_not_installed("measure")

	time <- seq(1, 10, by = 0.5)
	test_data <- tibble::tibble(sample_id = "test")
	test_data$lals <- measure::new_measure_list(
		list(measure::new_measure_tbl(location = time, value = dnorm(time, 5, 1)))
	)
	test_data$ri <- measure::new_measure_list(
		list(measure::new_measure_tbl(location = time, value = dnorm(time, 5, 1)))
	)

	rec <- recipes::recipe(~., data = test_data) |>
		step_sec_lals(measures = "lals", concentration_col = "ri", dn_dc = 0.185)

	prepped <- suppressWarnings(recipes::prep(rec))
	tidy_result <- recipes::tidy(prepped, number = 1)

	expect_s3_class(tidy_result, "tbl_df")
	expect_true("lals_col" %in% names(tidy_result))
	expect_true("concentration_col" %in% names(tidy_result))
	expect_true("angle" %in% names(tidy_result))
	expect_true("dn_dc" %in% names(tidy_result))
	expect_true("id" %in% names(tidy_result))
})

test_that("tidy.step_sec_rals returns expected structure", {
	skip_if_not_installed("measure")

	time <- seq(1, 10, by = 0.5)
	test_data <- tibble::tibble(sample_id = "test")
	test_data$rals <- measure::new_measure_list(
		list(measure::new_measure_tbl(location = time, value = dnorm(time, 5, 1)))
	)
	test_data$ri <- measure::new_measure_list(
		list(measure::new_measure_tbl(location = time, value = dnorm(time, 5, 1)))
	)

	rec <- recipes::recipe(~., data = test_data) |>
		step_sec_rals(measures = "rals", concentration_col = "ri", dn_dc = 0.185)

	prepped <- suppressWarnings(recipes::prep(rec))
	tidy_result <- recipes::tidy(prepped, number = 1)

	expect_s3_class(tidy_result, "tbl_df")
	expect_true("rals_col" %in% names(tidy_result))
	expect_true("concentration_col" %in% names(tidy_result))
	expect_true("angle" %in% names(tidy_result))
	expect_true("dn_dc" %in% names(tidy_result))
	expect_true("id" %in% names(tidy_result))
})

test_that("tidy.step_sec_dls returns expected structure", {
	skip_if_not_installed("measure")

	time <- seq(1, 10, by = 0.5)
	test_data <- tibble::tibble(sample_id = "test")
	test_data$dls <- measure::new_measure_list(
		list(measure::new_measure_tbl(location = time, value = dnorm(time, 5, 1)))
	)

	rec <- recipes::recipe(~., data = test_data) |>
		step_sec_dls(measures = "dls", temperature = 25)

	prepped <- recipes::prep(rec)
	tidy_result <- recipes::tidy(prepped, number = 1)

	expect_s3_class(tidy_result, "tbl_df")
	expect_true("dls_col" %in% names(tidy_result))
	expect_true("temperature" %in% names(tidy_result))
	expect_true("angle" %in% names(tidy_result))
	expect_true("output_cols" %in% names(tidy_result))
	expect_true("id" %in% names(tidy_result))
})

# ==============================================================================
# Calibration constant warning tests
# ==============================================================================

test_that("step_sec_lals warns when calibration_constant is NULL", {
	skip_if_not_installed("measure")

	time <- seq(1, 10, by = 0.5)
	test_data <- tibble::tibble(sample_id = "test")
	test_data$lals <- measure::new_measure_list(
		list(measure::new_measure_tbl(location = time, value = dnorm(time, 5, 1)))
	)
	test_data$ri <- measure::new_measure_list(
		list(measure::new_measure_tbl(location = time, value = dnorm(time, 5, 1)))
	)

	rec <- recipes::recipe(~., data = test_data) |>
		step_sec_lals(measures = "lals", concentration_col = "ri", dn_dc = 0.185)

	# Warning is raised during prep (which internally bakes training data)
	expect_warning(
		recipes::prep(rec),
		"calibration_constant"
	)
})

test_that("step_sec_rals warns when calibration_constant is NULL", {
	skip_if_not_installed("measure")

	time <- seq(1, 10, by = 0.5)
	test_data <- tibble::tibble(sample_id = "test")
	test_data$rals <- measure::new_measure_list(
		list(measure::new_measure_tbl(location = time, value = dnorm(time, 5, 1)))
	)
	test_data$ri <- measure::new_measure_list(
		list(measure::new_measure_tbl(location = time, value = dnorm(time, 5, 1)))
	)

	rec <- recipes::recipe(~., data = test_data) |>
		step_sec_rals(measures = "rals", concentration_col = "ri", dn_dc = 0.185)

	# Warning is raised during prep (which internally bakes training data)
	expect_warning(
		recipes::prep(rec),
		"calibration_constant"
	)
})
