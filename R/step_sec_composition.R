# ==============================================================================
# step_sec_composition
#
# Calculate copolymer composition from UV/RI signals
# ==============================================================================

#' Calculate Copolymer Composition from Detector Signals
#'
#' `step_sec_composition()` creates a *specification* of a recipe step that
#' calculates the weight fraction of components in a copolymer or blend using
#' UV and RI detector signals with known response factors.
#'
#' @param recipe A recipe object.
#' @param uv_col Name of the UV detector measure column.
#' @param ri_col Name of the RI detector measure column.
#' @param component_a_uv UV response factor for component A (extinction coefficient
#'   in mL/(mg*cm) or relative units).
#' @param component_a_ri RI response factor for component A (dn/dc in mL/g or
#'   relative units).
#' @param component_b_uv UV response factor for component B.
#' @param component_b_ri RI response factor for component B.
#' @param output_col Name for the output composition column. Default is
#'   `"composition_a"`. Contains weight fraction of component A (0-1).
#' @param min_signal Minimum signal threshold (as fraction of max) below which
#'   composition is set to NA. Default is 0.01.
#' @param clip Logical. Clip composition values to \[0, 1\] range? Default is `TRUE`.
#' @param role Role for generated columns.
#' @param trained Logical indicating if the step has been trained.
#' @param skip Logical. Should the step be skipped when baking?
#' @param id Unique step identifier.
#'
#' @return An updated recipe with the new step added.
#'
#' @details
#' For a two-component system, the detector signals are:
#'
#' \deqn{UV = \varepsilon_A \cdot c_A + \varepsilon_B \cdot c_B}
#' \deqn{RI = (dn/dc)_A \cdot c_A + (dn/dc)_B \cdot c_B}
#'
#' where c_A and c_B are the concentrations of components A and B.
#'
#' The weight fraction of component A is calculated by solving this system:
#'
#' \deqn{w_A = \frac{R_{obs} - R_B}{R_A - R_B}}
#'
#' where R_obs is the observed UV/RI ratio, and R_A, R_B are the ratios for
#' pure components (e/dn/dc).
#'
#' **Common Applications:**
#' \itemize{
#'   \item Styrene-acrylate copolymers (styrene is UV-active)
#'   \item Block copolymers with different chromophore content
#'   \item PEGylated proteins (protein at 280nm, PEG is UV-transparent)
#'   \item Polymer blends with known compositions
#' }
#'
#' **Example Response Factors:**
#' \itemize{
#'   \item Polystyrene: UV (254nm) ~ 1.0, dn/dc ~ 0.185
#'   \item PMMA: UV (254nm) ~ 0.01, dn/dc ~ 0.084
#'   \item PEG: UV (280nm) ~ 0, dn/dc ~ 0.135
#'   \item Proteins: UV (280nm) ~ 1.0, dn/dc ~ 0.185
#' }
#'
#' @note
#' Response factors must be in consistent units. The absolute values don't
#' matter as long as the ratios are correct for pure components.
#'
#' @family sec-composition
#' @export
#'
#' @examples
#' \dontrun{
#' library(recipes)
#' library(measure)
#'
#' # Styrene-MMA copolymer composition
#' rec <- recipe(~., data = copolymer_data) |>
#'   step_measure_input_long(ri_signal, location = vars(elution_time), col_name = "ri") |>
#'   step_measure_input_long(uv_signal, location = vars(elution_time), col_name = "uv") |>
#'   step_sec_baseline() |>
#'   step_sec_composition(
#'     uv_col = "uv",
#'     ri_col = "ri",
#'     component_a_uv = 1.0,    # Styrene (UV-active)
#'     component_a_ri = 0.185,  # Styrene dn/dc
#'     component_b_uv = 0.01,   # MMA (weak UV)
#'     component_b_ri = 0.084,  # MMA dn/dc
#'     output_col = "styrene_fraction"
#'   ) |>
#'   prep()
#' }
step_sec_composition <- function(
	recipe,
	uv_col = NULL,
	ri_col = NULL,
	component_a_uv,
	component_a_ri,
	component_b_uv,
	component_b_ri,
	output_col = "composition_a",
	min_signal = 0.01,
	clip = TRUE,
	role = NA,
	trained = FALSE,
	skip = FALSE,
	id = recipes::rand_id("sec_composition")
) {
	# Validate response factors
	if (
		missing(component_a_uv) ||
			missing(component_a_ri) ||
			missing(component_b_uv) ||
			missing(component_b_ri)
	) {
		cli::cli_abort(
			"All response factors must be specified: {.arg component_a_uv}, {.arg component_a_ri}, {.arg component_b_uv}, {.arg component_b_ri}."
		)
	}

	# Calculate response ratios for pure components
	ratio_a <- component_a_uv / component_a_ri
	ratio_b <- component_b_uv / component_b_ri

	if (abs(ratio_a - ratio_b) < 1e-10) {
		cli::cli_abort(
			c(
				"Components have identical UV/RI ratios.",
				"i" = "Composition cannot be determined when response ratios are equal."
			)
		)
	}

	recipes::add_step(
		recipe,
		step_sec_composition_new(
			uv_col = uv_col,
			ri_col = ri_col,
			component_a_uv = component_a_uv,
			component_a_ri = component_a_ri,
			component_b_uv = component_b_uv,
			component_b_ri = component_b_ri,
			ratio_a = ratio_a,
			ratio_b = ratio_b,
			output_col = output_col,
			min_signal = min_signal,
			clip = clip,
			role = role,
			trained = trained,
			skip = skip,
			id = id
		)
	)
}

step_sec_composition_new <- function(
	uv_col,
	ri_col,
	component_a_uv,
	component_a_ri,
	component_b_uv,
	component_b_ri,
	ratio_a,
	ratio_b,
	output_col,
	min_signal,
	clip,
	role,
	trained,
	skip,
	id
) {
	recipes::step(
		subclass = "sec_composition",
		uv_col = uv_col,
		ri_col = ri_col,
		component_a_uv = component_a_uv,
		component_a_ri = component_a_ri,
		component_b_uv = component_b_uv,
		component_b_ri = component_b_ri,
		ratio_a = ratio_a,
		ratio_b = ratio_b,
		output_col = output_col,
		min_signal = min_signal,
		clip = clip,
		role = role,
		trained = trained,
		skip = skip,
		id = id
	)
}

#' @export
prep.step_sec_composition <- function(x, training, info = NULL, ...) {
	check_for_measure(training)

	measure_cols <- find_measure_cols(training)

	# Find UV column if not specified
	if (is.null(x$uv_col)) {
		uv_cols <- measure_cols[grepl("uv", measure_cols, ignore.case = TRUE)]
		if (length(uv_cols) == 0) {
			cli::cli_abort("No UV column found. Specify {.arg uv_col} explicitly.")
		}
		uv_col <- uv_cols[1]
	} else {
		uv_col <- x$uv_col
		if (!uv_col %in% measure_cols) {
			cli::cli_abort("UV column {.val {uv_col}} not found in measure columns.")
		}
	}

	# Find RI column if not specified
	if (is.null(x$ri_col)) {
		ri_cols <- measure_cols[grepl("ri", measure_cols, ignore.case = TRUE)]
		if (length(ri_cols) == 0) {
			cli::cli_abort("No RI column found. Specify {.arg ri_col} explicitly.")
		}
		ri_col <- ri_cols[1]
	} else {
		ri_col <- x$ri_col
		if (!ri_col %in% measure_cols) {
			cli::cli_abort("RI column {.val {ri_col}} not found in measure columns.")
		}
	}

	step_sec_composition_new(
		uv_col = uv_col,
		ri_col = ri_col,
		component_a_uv = x$component_a_uv,
		component_a_ri = x$component_a_ri,
		component_b_uv = x$component_b_uv,
		component_b_ri = x$component_b_ri,
		ratio_a = x$ratio_a,
		ratio_b = x$ratio_b,
		output_col = x$output_col,
		min_signal = x$min_signal,
		clip = x$clip,
		role = x$role,
		trained = TRUE,
		skip = x$skip,
		id = x$id
	)
}

#' @export
bake.step_sec_composition <- function(object, new_data, ...) {
	uv_col <- object$uv_col
	ri_col <- object$ri_col
	ratio_a <- object$ratio_a
	ratio_b <- object$ratio_b
	output_col <- object$output_col
	min_signal <- object$min_signal
	clip <- object$clip

	# Calculate composition for each sample
	comp_list <- purrr::map2(
		new_data[[uv_col]],
		new_data[[ri_col]],
		function(uv_m, ri_m) {
			uv_val <- uv_m$value
			ri_val <- ri_m$value
			location <- uv_m$location

			# Determine signal threshold
			max_ri <- max(abs(ri_val), na.rm = TRUE)
			threshold <- min_signal * max_ri

			# Calculate observed ratio and composition
			composition <- rep(NA_real_, length(uv_val))
			valid <- abs(ri_val) > threshold & !is.na(uv_val) & !is.na(ri_val)

			if (any(valid)) {
				observed_ratio <- uv_val[valid] / ri_val[valid]

				# Calculate weight fraction of component A
				# w_A = (R_obs - R_B) / (R_A - R_B)
				composition[valid] <- (observed_ratio - ratio_b) / (ratio_a - ratio_b)

				# Clip to \[0, 1\] if requested
				if (clip) {
					composition[valid] <- pmax(0, pmin(1, composition[valid]))
				}
			}

			# Return as measure object
			new_measure_tbl(location = location, value = composition)
		}
	)

	new_data[[output_col]] <- new_measure_list(comp_list)

	tibble::as_tibble(new_data)
}

#' @export
print.step_sec_composition <- function(
	x,
	width = max(20, options()$width - 30),
	...
) {
	title <- "SEC composition analysis"
	if (x$trained) {
		cat(
			title,
			" (",
			x$uv_col,
			"/",
			x$ri_col,
			" -> ",
			x$output_col,
			")",
			sep = ""
		)
	} else {
		cat(title)
	}
	cat("\n")
	invisible(x)
}

#' @rdname tidy.step_sec
#' @export
#' @keywords internal
tidy.step_sec_composition <- function(x, ...) {
	tibble::tibble(
		uv_col = x$uv_col %||% NA_character_,
		ri_col = x$ri_col %||% NA_character_,
		output_col = x$output_col,
		ratio_a = x$ratio_a,
		ratio_b = x$ratio_b,
		id = x$id
	)
}

#' @rdname required_pkgs.step_sec
#' @export
#' @keywords internal
required_pkgs.step_sec_composition <- function(x, ...) {
	c("measure.sec", "measure")
}
