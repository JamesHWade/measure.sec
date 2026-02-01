# ==============================================================================
# step_sec_lals
#
# Low-angle light scattering processing
# ==============================================================================

#' Low-Angle Light Scattering Processing for SEC
#'
#' `step_sec_lals()` creates a *specification* of a recipe step that processes
#' low-angle light scattering (LALS) signals for absolute molecular weight.
#'
#' @param recipe A recipe object.
#' @param measures Character vector of LALS measure columns. If `NULL`, the
#'   step searches for measure columns containing "lals".
#' @param concentration_col Name of the concentration measure column (from
#'   `step_sec_concentration()` or similar).
#' @param angle Detection angle in degrees (must be < 20 for LALS).
#' @param laser_wavelength Laser wavelength in nm.
#' @param dn_dc Refractive index increment (mL/g). Required unless
#'   `optical_constant` is provided.
#' @param solvent_ri Solvent refractive index. Default is 1.333 (water).
#' @param optical_constant Optional optical constant K; overrides dn/dc.
#' @param calibration_constant LALS instrument calibration constant. If `NULL`,
#'   results are in relative units.
#' @param output_mw Name for the molecular weight output column.
#' @param min_signal Minimum signal threshold (as fraction of max) below which
#'   MW is set to NA.
#' @param role Role for generated columns.
#' @param trained Logical indicating if the step has been trained.
#' @param skip Logical. Should the step be skipped when baking?
#' @param id Unique step identifier.
#'
#' @return An updated recipe with the new step added.
#'
#' @details
#' LALS provides absolute MW using a low-angle detector (e.g., ~7 degrees).
#' It assumes P(theta) ~ 1 for small particles and is most accurate when
#' angular dependence is minimal.
#'
#' **When to use LALS vs MALS:**
#' \itemize{
#'   \item **LALS**: Preferred for smaller molecules (Rg < ~10 nm) or when
#'     multi-angle data is not available. Single-angle measurement is faster
#'     and simpler but cannot determine Rg.
#'   \item **MALS**: Required for large molecules where angular dependence is
#'     significant. Provides both Mw and Rg from extrapolation to zero angle.
#' }
#'
#' @family sec-detectors
#' @export
#'
#' @examples
#' \dontrun{
#' library(recipes)
#' library(measure)
#'
#' rec <- recipe(~., data = sec_triple_detect) |>
#'   step_measure_input_long(ri_signal, location = vars(elution_time), col_name = "ri") |>
#'   step_measure_input_long(lals_signal, location = vars(elution_time), col_name = "lals") |>
#'   step_sec_concentration(detector = "ri", injection_mass = 0.2, flow_rate = 1.0) |>
#'   step_sec_lals(measures = "lals", concentration_col = "ri", dn_dc = 0.185) |>
#'   prep()
#' }
step_sec_lals <- function(
	recipe,
	measures = NULL,
	concentration_col = NULL,
	angle = 7,
	laser_wavelength = 670,
	dn_dc = NULL,
	solvent_ri = 1.333,
	optical_constant = NULL,
	calibration_constant = NULL,
	output_mw = "mw_lals",
	min_signal = 0.01,
	role = NA,
	trained = FALSE,
	skip = FALSE,
	id = recipes::rand_id("sec_lals")
) {
	if (!is.null(optical_constant) && !is.null(dn_dc)) {
		cli::cli_warn(
			"{.arg optical_constant} takes precedence over {.arg dn_dc}."
		)
	}

	if (!is.numeric(angle) || angle <= 0 || angle >= 20) {
		cli::cli_abort("{.arg angle} must be between 0 and 20 degrees for LALS.")
	}

	if (is.null(optical_constant) && is.null(dn_dc)) {
		cli::cli_abort(
			"Either {.arg dn_dc} or {.arg optical_constant} must be provided."
		)
	}

	recipes::add_step(
		recipe,
		step_sec_lals_new(
			measures = measures,
			concentration_col = concentration_col,
			angle = angle,
			laser_wavelength = laser_wavelength,
			dn_dc = dn_dc,
			solvent_ri = solvent_ri,
			optical_constant = optical_constant,
			calibration_constant = calibration_constant,
			output_mw = output_mw,
			min_signal = min_signal,
			role = role,
			trained = trained,
			skip = skip,
			id = id
		)
	)
}

step_sec_lals_new <- function(
	measures,
	concentration_col,
	angle,
	laser_wavelength,
	dn_dc,
	solvent_ri,
	optical_constant,
	calibration_constant,
	output_mw,
	min_signal,
	role,
	trained,
	skip,
	id
) {
	recipes::step(
		subclass = "sec_lals",
		measures = measures,
		concentration_col = concentration_col,
		angle = angle,
		laser_wavelength = laser_wavelength,
		dn_dc = dn_dc,
		solvent_ri = solvent_ri,
		optical_constant = optical_constant,
		calibration_constant = calibration_constant,
		output_mw = output_mw,
		min_signal = min_signal,
		role = role,
		trained = trained,
		skip = skip,
		id = id
	)
}

#' @export
prep.step_sec_lals <- function(x, training, info = NULL, ...) {
	check_for_measure(training)

	measure_cols <- find_measure_cols(training)

	if (is.null(x$measures)) {
		lals_cols <- measure_cols[grepl("lals", measure_cols, ignore.case = TRUE)]
		if (length(lals_cols) == 0) {
			cli::cli_abort(
				"No LALS column found. Specify {.arg measures} explicitly."
			)
		}
		measures <- lals_cols[1]
	} else {
		measures <- x$measures
	}

	if (length(measures) != 1) {
		cli::cli_abort("{.arg measures} must specify a single LALS column.")
	}

	if (!measures %in% measure_cols) {
		cli::cli_abort(
			"LALS column {.val {measures}} not found in measure columns."
		)
	}

	if (is.null(x$concentration_col)) {
		conc_cols <- measure_cols[
			grepl(
				"ri|conc",
				measure_cols,
				ignore.case = TRUE
			)
		]
		if (length(conc_cols) == 0) {
			cli::cli_abort(
				"No concentration column found. Specify {.arg concentration_col} explicitly."
			)
		}
		concentration_col <- conc_cols[1]
	} else {
		concentration_col <- x$concentration_col
		if (!concentration_col %in% measure_cols) {
			cli::cli_abort(
				"Concentration column {.val {concentration_col}} not found in measure columns."
			)
		}
	}

	step_sec_lals_new(
		measures = measures,
		concentration_col = concentration_col,
		angle = x$angle,
		laser_wavelength = x$laser_wavelength,
		dn_dc = x$dn_dc,
		solvent_ri = x$solvent_ri,
		optical_constant = x$optical_constant,
		calibration_constant = x$calibration_constant,
		output_mw = x$output_mw,
		min_signal = x$min_signal,
		role = x$role,
		trained = TRUE,
		skip = x$skip,
		id = x$id
	)
}

#' @export
bake.step_sec_lals <- function(object, new_data, ...) {
	lals_col <- object$measures
	concentration_col <- object$concentration_col
	dn_dc <- object$dn_dc
	solvent_ri <- object$solvent_ri
	laser_wavelength <- object$laser_wavelength
	calibration_constant <- object$calibration_constant
	min_signal <- object$min_signal
	output_mw <- object$output_mw

	if (is.null(calibration_constant)) {
		cli::cli_warn(
			c(
				"No {.arg calibration_constant} provided.",
				"i" = "Results will be in relative units, not absolute MW.",
				"i" = "Provide a calibration constant for absolute values."
			)
		)
		calibration_constant <- 1.0
	}

	if (is.null(object$optical_constant)) {
		K <- .optical_constant(dn_dc, solvent_ri, laser_wavelength)
	} else {
		K <- object$optical_constant
	}

	mw_list <- purrr::pmap(
		list(new_data[[lals_col]], new_data[[concentration_col]]),
		function(lals_m, conc_m) {
			lals_val <- lals_m$value
			conc_val <- conc_m$value
			location <- lals_m$location

			max_conc <- max(abs(conc_val), na.rm = TRUE)
			threshold <- min_signal * max_conc

			mw <- rep(NA_real_, length(lals_val))
			valid <- abs(conc_val) > threshold &
				lals_val > 0 &
				!is.na(lals_val) &
				!is.na(conc_val)

			if (any(valid)) {
				R_theta <- lals_val[valid] * calibration_constant
				mw[valid] <- R_theta / (K * conc_val[valid])

				mw[valid][mw[valid] < 100] <- NA_real_
				mw[valid][mw[valid] > 1e10] <- NA_real_
			}

			new_measure_tbl(location = location, value = mw)
		}
	)

	new_data[[output_mw]] <- new_measure_list(mw_list)

	tibble::as_tibble(new_data)
}

#' @export
print.step_sec_lals <- function(
	x,
	width = max(20, options()$width - 30),
	...
) {
	title <- "SEC LALS processing"
	if (x$trained) {
		cat(title, " on ", x$measures, " -> ", x$output_mw, sep = "")
	} else {
		cat(title)
	}
	cat("\n")
	invisible(x)
}

#' @rdname tidy.step_sec
#' @export
#' @keywords internal
tidy.step_sec_lals <- function(x, ...) {
	tibble::tibble(
		lals_col = x$measures %||% NA_character_,
		concentration_col = x$concentration_col %||% NA_character_,
		angle = x$angle,
		laser_wavelength = x$laser_wavelength,
		dn_dc = x$dn_dc %||% NA_real_,
		output_mw = x$output_mw,
		id = x$id
	)
}

#' @rdname required_pkgs.step_sec
#' @export
#' @keywords internal
required_pkgs.step_sec_lals <- function(x, ...) {
	c("measure.sec", "measure")
}
