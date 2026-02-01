# ==============================================================================
# step_sec_uv
#
# UV detector processing with extinction coefficient handling
# ==============================================================================

#' UV Detector Processing for SEC
#'
#' `step_sec_uv()` creates a *specification* of a recipe step that processes
#' UV detector signals for SEC analysis, including application of extinction
#' coefficients for concentration determination.
#'
#' @param recipe A recipe object.
#' @param measures Character vector of UV detector column names to process.
#'   If `NULL`, will look for columns containing "uv" in the name.
#' @param extinction_coef Extinction coefficient in mL/(mg*cm) or L/(g*cm).
#'   Can be:
#'   - A single numeric value applied to all samples
#'   - `NULL` to skip normalization (signal remains in AU)
#' @param extinction_column Character name of a column containing sample-specific
#'   extinction coefficients. Overrides `extinction_coef` if provided.
#' @param wavelength UV detection wavelength in nm. For documentation only.
#' @param path_length Path length of the flow cell in cm. Default is 1.0.
#' @param output_col Name for the output column. Default is to modify in place.
#' @param role Role for generated columns.
#' @param trained Logical indicating if the step has been trained.
#' @param skip Logical. Should the step be skipped when baking?
#' @param id Unique step identifier.
#'
#' @return An updated recipe with the new step added.
#'
#' @details
#' The UV detector measures absorbance according to the Beer-Lambert law:
#'
#' \deqn{A = \varepsilon \times c \times l}
#'
#' where:
#' - A is absorbance (AU)
#' - epsilon is the molar extinction coefficient in mL/(mg*cm)
#' - c is the concentration (mg/mL)
#' - l is the path length (cm)
#'
#' This step can convert UV absorbance to concentration-proportional signals
#' by dividing by the extinction coefficient and path length.
#'
#' **Common UV applications in SEC:**
#' - Proteins at 280 nm (aromatic amino acids)
#' - Nucleic acids at 260 nm
#' - Conjugated polymers
#' - UV-active end groups or labels
#'
#' **UV vs RI for concentration:**
#' - UV is more sensitive for chromophore-containing analytes
#' - UV response depends on chemical composition (may vary with MW)
#' - RI is more universal but less sensitive
#' - For accurate MW, combine both detectors
#'
#' @family sec-detectors
#' @export
#'
#' @examples
#' \dontrun{
#' library(recipes)
#' library(measure)
#'
#' # Apply fixed extinction coefficient
#' rec <- recipe(~., data = sec_triple_detect) |>
#'   step_measure_input_long(uv_signal, location = vars(elution_time), col_name = "uv") |>
#'   step_sec_uv(extinction_coef = 1.0, wavelength = 280) |>
#'   prep()
#'
#' # Use sample-specific extinction coefficients
#' rec <- recipe(~., data = sec_triple_detect) |>
#'   step_measure_input_long(uv_signal, location = vars(elution_time), col_name = "uv") |>
#'   step_sec_uv(extinction_column = "ext_coef") |>
#'   prep()
#' }
step_sec_uv <- function(
	recipe,
	measures = NULL,
	extinction_coef = NULL,
	extinction_column = NULL,
	wavelength = 280,
	path_length = 1.0,
	output_col = NULL,
	role = NA,
	trained = FALSE,
	skip = FALSE,
	id = recipes::rand_id("sec_uv")
) {
	# Validate inputs
	if (!is.null(extinction_coef) && !is.null(extinction_column)) {
		cli::cli_warn(
			"{.arg extinction_column} takes precedence over {.arg extinction_coef}."
		)
	}

	if (!is.null(extinction_coef)) {
		if (
			!is.numeric(extinction_coef) ||
				length(extinction_coef) != 1 ||
				extinction_coef <= 0
		) {
			cli::cli_abort("{.arg extinction_coef} must be a positive numeric value.")
		}
	}

	if (!is.numeric(path_length) || path_length <= 0) {
		cli::cli_abort("{.arg path_length} must be a positive numeric value.")
	}

	recipes::add_step(
		recipe,
		step_sec_uv_new(
			measures = measures,
			extinction_coef = extinction_coef,
			extinction_column = extinction_column,
			wavelength = wavelength,
			path_length = path_length,
			output_col = output_col,
			role = role,
			trained = trained,
			skip = skip,
			id = id
		)
	)
}

step_sec_uv_new <- function(
	measures,
	extinction_coef,
	extinction_column,
	wavelength,
	path_length,
	output_col,
	role,
	trained,
	skip,
	id
) {
	recipes::step(
		subclass = "sec_uv",
		measures = measures,
		extinction_coef = extinction_coef,
		extinction_column = extinction_column,
		wavelength = wavelength,
		path_length = path_length,
		output_col = output_col,
		role = role,
		trained = trained,
		skip = skip,
		id = id
	)
}

#' @export
prep.step_sec_uv <- function(x, training, info = NULL, ...) {
	check_for_measure(training)

	# Find UV columns if not specified
	if (is.null(x$measures)) {
		measure_cols <- find_measure_cols(training)
		# Look for columns containing "uv" (case insensitive)
		uv_cols <- measure_cols[grepl("uv", measure_cols, ignore.case = TRUE)]
		if (length(uv_cols) == 0) {
			cli::cli_abort(
				"No UV detector columns found. Specify {.arg measures} explicitly."
			)
		}
		measures <- uv_cols
	} else {
		measures <- x$measures
	}

	# Validate extinction_column exists if specified
	if (!is.null(x$extinction_column)) {
		if (!x$extinction_column %in% names(training)) {
			cli::cli_abort(
				"Column {.val {x$extinction_column}} not found in training data."
			)
		}
	}

	step_sec_uv_new(
		measures = measures,
		extinction_coef = x$extinction_coef,
		extinction_column = x$extinction_column,
		wavelength = x$wavelength,
		path_length = x$path_length,
		output_col = x$output_col,
		role = x$role,
		trained = TRUE,
		skip = x$skip,
		id = x$id
	)
}

#' @export
bake.step_sec_uv <- function(object, new_data, ...) {
	measures <- object$measures
	extinction_coef <- object$extinction_coef
	extinction_column <- object$extinction_column
	path_length <- object$path_length

	for (col in measures) {
		# Get extinction coefficient values for each sample
		if (!is.null(extinction_column)) {
			ext_values <- new_data[[extinction_column]]
		} else if (!is.null(extinction_coef)) {
			ext_values <- rep(extinction_coef, nrow(new_data))
		} else {
			ext_values <- rep(1.0, nrow(new_data)) # No normalization
		}

		# Process each sample
		new_data[[col]] <- new_measure_list(
			purrr::map2(new_data[[col]], ext_values, function(m, ext) {
				# Convert absorbance to concentration-proportional signal
				# A = epsilon * c * l, so c ~ A / (epsilon * l)
				if (!is.na(ext) && ext > 0) {
					m$value <- m$value / (ext * path_length)
				}
				m
			})
		)
	}

	tibble::as_tibble(new_data)
}

#' @export
print.step_sec_uv <- function(
	x,
	width = max(20, options()$width - 30),
	...
) {
	title <- paste0("SEC UV detector processing (", x$wavelength, " nm)")
	if (!is.null(x$extinction_coef)) {
		title <- paste0(title, ", e = ", x$extinction_coef)
	} else if (!is.null(x$extinction_column)) {
		title <- paste0(title, ", e from ", x$extinction_column)
	}

	if (x$trained) {
		cols_str <- paste(x$measures, collapse = ", ")
		cat(title, " on ", cols_str, sep = "")
	} else {
		cat(title)
	}
	cat("\n")
	invisible(x)
}

#' @rdname tidy.step_sec
#' @export
#' @keywords internal
tidy.step_sec_uv <- function(x, ...) {
	tibble::tibble(
		measures = list(x$measures),
		extinction_coef = x$extinction_coef %||% NA_real_,
		extinction_column = x$extinction_column %||% NA_character_,
		wavelength = x$wavelength,
		path_length = x$path_length,
		id = x$id
	)
}

#' @rdname required_pkgs.step_sec
#' @export
#' @keywords internal
required_pkgs.step_sec_uv <- function(x, ...) {
	c("measure.sec", "measure")
}
