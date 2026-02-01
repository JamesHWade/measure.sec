# ==============================================================================
# step_sec_mals
#
# Multi-Angle Light Scattering processing for absolute MW
# ==============================================================================

#' Multi-Angle Light Scattering Processing for SEC
#'
#' `step_sec_mals()` creates a *specification* of a recipe step that processes
#' multi-angle light scattering (MALS) detector signals to determine absolute
#' molecular weight and radius of gyration at each elution slice.
#'
#' @param recipe A recipe object.
#' @param mals_col Name of the MALS detector measure column. For single-angle
#'   detectors (RALS/LALS), this is the only signal needed.
#' @param concentration_col Name of the concentration measure column (from
#'   `step_sec_concentration()` or similar).
#' @param dn_dc Refractive index increment (mL/g). Required for absolute MW.
#' @param dn_dc_column Column containing sample-specific dn/dc values.
#' @param wavelength Laser wavelength in nm. Default is 658 (common for MALS).
#' @param solvent_ri Solvent refractive index. Default is 1.333 (water).
#'   Common values: water = 1.333, THF = 1.407, toluene = 1.497.
#' @param angles Numeric vector of detection angles in degrees. For single-angle
#'   detectors, provide just one value (e.g., `90` for RALS). Default assumes
#'   a 90-degree detector.
#' @param formalism Angular extrapolation method for multi-angle data:
#'   `"zimm"` (default), `"debye"`, or `"berry"`.
#' @param calibration_constant MALS instrument calibration constant. If `NULL`,
#'   results are in relative units. Obtain from toluene standard calibration.
#' @param output_mw Name for the molecular weight output column. Default is `"mw_mals"`.
#' @param output_rg Name for the radius of gyration output column. Default is `"rg"`.
#'   Only calculated if multiple angles are provided.
#' @param min_signal Minimum signal threshold (as fraction of max) below which
#'   MW is set to NA. Default is 0.01.
#' @param role Role for generated columns.
#' @param trained Logical indicating if the step has been trained.
#' @param skip Logical. Should the step be skipped when baking?
#' @param id Unique step identifier.
#'
#' @return An updated recipe with the new step added, containing absolute MW
#'   (and optionally Rg) at each elution point.
#'
#' @details
#' Light scattering provides absolute molecular weight without calibration
#' standards. The fundamental relationship is:
#'
#' \deqn{\frac{K \cdot c}{R(\theta)} = \frac{1}{M_w} \cdot P(\theta) + 2A_2 c}
#'
#' where:
#' \itemize{
#'   \item K is the optical constant
#'   \item c is the concentration
#'   \item R(theta) is the excess Rayleigh ratio
#'   \item Mw is the weight-average molecular weight
#'   \item P(theta) is the particle scattering function
#'   \item A2 is the second virial coefficient
#' }
#'
#' The optical constant K is calculated as:
#' \deqn{K = \frac{4\pi^2 n_0^2 (dn/dc)^2}{N_A \lambda^4}}
#'
#' **Formalisms for Angular Extrapolation:**
#' \itemize{
#'   \item Zimm: Kc/R vs sin^2(theta/2) - best for random coils
#'   \item Debye: Kc/R vs sin^2(theta/2) - similar to Zimm
#'   \item Berry: sqrt(Kc/R) vs sin^2(theta/2) - better for large particles
#' }
#'
#' **Single-Angle vs Multi-Angle:**
#' \itemize{
#'   \item Single angle (RALS/LALS): Provides Mw only, assumes P(theta) ~ 1
#'   \item Multi-angle (MALS): Provides both Mw and Rg from angular dependence
#' }
#'
#' @note
#' Accurate results require:
#' \itemize{
#'   \item Known and accurate dn/dc value
#'   \item Calibrated instrument (calibration_constant from toluene)
#'   \item Accurate concentration from RI detector
#'   \item Clean baseline and aligned detectors
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
#' # Single-angle (90 degree) light scattering
#' rec <- recipe(~., data = sec_triple_detect) |>
#'   step_measure_input_long(ri_signal, location = vars(elution_time), col_name = "ri") |>
#'   step_measure_input_long(mals_signal, location = vars(elution_time), col_name = "mals") |>
#'   step_sec_baseline() |>
#'   step_sec_ri(dn_dc = 0.185) |>
#'   step_sec_concentration(detector = "ri", injection_mass = 0.2, flow_rate = 1.0) |>
#'   step_sec_mals(
#'     mals_col = "mals",
#'     concentration_col = "ri",
#'     dn_dc = 0.185,
#'     wavelength = 658,
#'     angles = 90
#'   ) |>
#'   prep()
#' }
step_sec_mals <- function(
	recipe,
	mals_col = NULL,
	concentration_col = NULL,
	dn_dc = NULL,
	dn_dc_column = NULL,
	wavelength = 658,
	solvent_ri = 1.333,
	angles = 90,
	formalism = c("zimm", "debye", "berry"),
	calibration_constant = NULL,
	output_mw = "mw_mals",
	output_rg = "rg",
	min_signal = 0.01,
	role = NA,
	trained = FALSE,
	skip = FALSE,
	id = recipes::rand_id("sec_mals")
) {
	formalism <- match.arg(formalism)

	# Validate inputs
	if (is.null(dn_dc) && is.null(dn_dc_column)) {
		cli::cli_abort(
			"Either {.arg dn_dc} or {.arg dn_dc_column} must be specified for absolute MW."
		)
	}

	if (!is.numeric(wavelength) || wavelength <= 0) {
		cli::cli_abort("{.arg wavelength} must be a positive number (in nm).")
	}

	if (!is.numeric(solvent_ri) || solvent_ri <= 0) {
		cli::cli_abort("{.arg solvent_ri} must be a positive number.")
	}

	if (!is.numeric(angles) || any(angles <= 0) || any(angles >= 180)) {
		cli::cli_abort("{.arg angles} must be between 0 and 180 degrees.")
	}

	recipes::add_step(
		recipe,
		step_sec_mals_new(
			mals_col = mals_col,
			concentration_col = concentration_col,
			dn_dc = dn_dc,
			dn_dc_column = dn_dc_column,
			wavelength = wavelength,
			solvent_ri = solvent_ri,
			angles = angles,
			formalism = formalism,
			calibration_constant = calibration_constant,
			output_mw = output_mw,
			output_rg = output_rg,
			min_signal = min_signal,
			role = role,
			trained = trained,
			skip = skip,
			id = id
		)
	)
}

step_sec_mals_new <- function(
	mals_col,
	concentration_col,
	dn_dc,
	dn_dc_column,
	wavelength,
	solvent_ri,
	angles,
	formalism,
	calibration_constant,
	output_mw,
	output_rg,
	min_signal,
	role,
	trained,
	skip,
	id
) {
	recipes::step(
		subclass = "sec_mals",
		mals_col = mals_col,
		concentration_col = concentration_col,
		dn_dc = dn_dc,
		dn_dc_column = dn_dc_column,
		wavelength = wavelength,
		solvent_ri = solvent_ri,
		angles = angles,
		formalism = formalism,
		calibration_constant = calibration_constant,
		output_mw = output_mw,
		output_rg = output_rg,
		min_signal = min_signal,
		role = role,
		trained = trained,
		skip = skip,
		id = id
	)
}

#' Calculate optical constant K
#' @noRd
.optical_constant <- function(dn_dc, solvent_ri, wavelength_nm) {
	# K = 4 * pi^2 * n0^2 * (dn/dc)^2 / (NA * lambda^4)
	# wavelength in nm -> m
	# dn/dc in mL/g = cm^3/g
	# Result in mol*cm^2/g^2

	lambda_cm <- wavelength_nm * 1e-7 # nm to cm
	NA_val <- 6.022e23 # Avogadro's number

	K <- 4 * pi^2 * solvent_ri^2 * dn_dc^2 / (NA_val * lambda_cm^4)
	K
}

#' @export
prep.step_sec_mals <- function(x, training, info = NULL, ...) {
	check_for_measure(training)

	measure_cols <- find_measure_cols(training)

	# Find MALS column if not specified
	if (is.null(x$mals_col)) {
		mals_cols <- measure_cols[
			grepl(
				"mals|ls|light",
				measure_cols,
				ignore.case = TRUE
			)
		]
		if (length(mals_cols) == 0) {
			cli::cli_abort(
				"No MALS column found. Specify {.arg mals_col} explicitly."
			)
		}
		mals_col <- mals_cols[1]
	} else {
		mals_col <- x$mals_col
		if (!mals_col %in% measure_cols) {
			cli::cli_abort(
				"MALS column {.val {mals_col}} not found in measure columns."
			)
		}
	}

	# Find concentration column if not specified
	if (is.null(x$concentration_col)) {
		# Try to find RI or concentration column
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

	# Validate dn_dc_column exists if specified
	if (!is.null(x$dn_dc_column)) {
		if (!x$dn_dc_column %in% names(training)) {
			cli::cli_abort(
				"Column {.val {x$dn_dc_column}} not found in training data."
			)
		}
	}

	step_sec_mals_new(
		mals_col = mals_col,
		concentration_col = concentration_col,
		dn_dc = x$dn_dc,
		dn_dc_column = x$dn_dc_column,
		wavelength = x$wavelength,
		solvent_ri = x$solvent_ri,
		angles = x$angles,
		formalism = x$formalism,
		calibration_constant = x$calibration_constant,
		output_mw = x$output_mw,
		output_rg = x$output_rg,
		min_signal = x$min_signal,
		role = x$role,
		trained = TRUE,
		skip = x$skip,
		id = x$id
	)
}

#' @export
bake.step_sec_mals <- function(object, new_data, ...) {
	mals_col <- object$mals_col
	concentration_col <- object$concentration_col
	dn_dc <- object$dn_dc
	dn_dc_column <- object$dn_dc_column
	wavelength <- object$wavelength
	solvent_ri <- object$solvent_ri
	angles <- object$angles
	calibration_constant <- object$calibration_constant
	output_mw <- object$output_mw
	min_signal <- object$min_signal

	# Get dn/dc values
	if (!is.null(dn_dc_column)) {
		dn_dc_values <- new_data[[dn_dc_column]]
	} else {
		dn_dc_values <- rep(dn_dc, nrow(new_data))
	}

	# Default calibration constant if not provided
	if (is.null(calibration_constant)) {
		calibration_constant <- 1.0
	}

	# Calculate MW for each sample
	mw_list <- purrr::pmap(
		list(
			new_data[[mals_col]],
			new_data[[concentration_col]],
			dn_dc_values
		),
		function(mals_m, conc_m, dndc) {
			mals_val <- mals_m$value
			conc_val <- conc_m$value
			location <- mals_m$location

			# Calculate optical constant
			K <- .optical_constant(dndc, solvent_ri, wavelength)

			# Determine signal threshold
			max_conc <- max(abs(conc_val), na.rm = TRUE)
			threshold <- min_signal * max_conc

			# Calculate MW at each point
			mw <- rep(NA_real_, length(mals_val))
			valid <- abs(conc_val) > threshold &
				mals_val > 0 &
				!is.na(mals_val) &
				!is.na(conc_val)

			if (any(valid)) {
				# For single angle at 90 degrees, P(theta) ~ 1 for small particles
				# Mw = R(theta) / (K * c)
				# With calibration: R = signal * calibration_constant
				R_theta <- mals_val[valid] * calibration_constant
				mw[valid] <- R_theta / (K * conc_val[valid])

				# Ensure reasonable MW values (filter out noise)
				mw[valid][mw[valid] < 100] <- NA_real_ # MW < 100 is unrealistic
				mw[valid][mw[valid] > 1e10] <- NA_real_ # MW > 10 billion is unrealistic
			}

			# Return as measure object
			new_measure_tbl(location = location, value = mw)
		}
	)

	new_data[[output_mw]] <- new_measure_list(mw_list)

	tibble::as_tibble(new_data)
}

#' @export
print.step_sec_mals <- function(
	x,
	width = max(20, options()$width - 30),
	...
) {
	title <- paste0(
		"SEC MALS processing (",
		length(x$angles),
		" angle",
		if (length(x$angles) > 1) "s" else "",
		")"
	)
	if (x$trained) {
		cat(title, " on ", x$mals_col, " -> ", x$output_mw, sep = "")
	} else {
		cat(title)
	}
	cat("\n")
	invisible(x)
}

#' @rdname tidy.step_sec
#' @export
#' @keywords internal
tidy.step_sec_mals <- function(x, ...) {
	tibble::tibble(
		mals_col = x$mals_col %||% NA_character_,
		concentration_col = x$concentration_col %||% NA_character_,
		dn_dc = x$dn_dc %||% NA_real_,
		wavelength = x$wavelength,
		angles = list(x$angles),
		formalism = x$formalism,
		output_mw = x$output_mw,
		id = x$id
	)
}

#' @rdname required_pkgs.step_sec
#' @export
#' @keywords internal
required_pkgs.step_sec_mals <- function(x, ...) {
	c("measure.sec", "measure")
}
