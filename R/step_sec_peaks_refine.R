# ==============================================================================
# step_sec_peaks_refine
#
# Tighten peak boundaries from step_sec_peaks_detect using height-fraction
# or other methods. Addresses the common issue where finderskeepers boundaries
# extend to full baseline return, including noisy tails that bias MW results.
# ==============================================================================

#' Refine SEC/GPC Peak Boundaries
#'
#' `step_sec_peaks_refine()` creates a *specification* of a recipe step that
#' tightens peak boundaries set by [step_sec_peaks_detect()]. The default
#' `"height_fraction"` method sets boundaries where the corrected signal drops
#' below a fraction of each peak's apex height, matching the approach used by
#' industry GPC software (e.g., Waters Empower's 0.5% default).
#'
#' @param recipe A recipe object. The step will be added to the sequence of
#'   operations for this recipe.
#' @param method Boundary refinement method. Currently only
#'   `"height_fraction"` is supported:
#'   - `"height_fraction"` (default): Set boundaries where the corrected signal
#'     drops below `cutoff * apex_height`.
#' @param cutoff Numeric threshold for the chosen method. For
#'   `"height_fraction"`, this is the fraction of apex height below which
#'   the signal is considered outside the peak. Default is `0.005` (0.5%).
#'   Common values:
#'   - `0.005`: 0.5% of peak height (Empower default, most common)
#'   - `0.01`: 1% (more aggressive trimming)
#'   - `0.001`: 0.1% (looser, closer to baseline-return)
#' @param peaks_col Name of the peaks column to refine. Default is `".peaks"`.
#' @param measures_col Name of the measure column containing the
#'   baseline-corrected signal. Default is `".measures"`. If not found,
#'   the first available measure column is used automatically.
#' @param role Not used by this step since no new variables are created.
#' @param trained A logical to indicate if the quantities for preprocessing
#'   have been estimated.
#' @param skip A logical. Should the step be skipped when the recipe is baked?
#' @param id A character string that is unique to this step to identify it.
#'
#' @return An updated version of `recipe` with the new step added to the
#'   sequence of any existing operations. The `.peaks` column will be updated
#'   with refined `left_base` and `right_base` values. Two new columns are
#'   added to each `peaks_tbl`: `original_left_base` and `original_right_base`
#'   preserving the boundaries from [step_sec_peaks_detect()].
#'
#' @details
#' Wide peak boundaries from baseline-return detection cause systematic errors
#' in molecular weight calculations:
#'
#' - **Mn bias**: Low-signal noise at high elution volumes contributes small
#'   M_i values, pulling Mn downward.
#' - **PDI inflation**: With Mw stable but Mn depressed, dispersity (Mw/Mn)
#'   increases artificially.
#' - **Noise integration**: Baseline noise fluctuating around zero gets included
#'   when `corrected > 0`.
#'
#' The `"height_fraction"` method addresses this by walking inward from each
#' boundary until the signal exceeds a fraction of the peak's apex height.
#'
#' This step should run after [step_sec_baseline()] and
#' [step_sec_peaks_detect()]. The signal in `measures_col` should be
#' baseline-corrected for best results.
#'
#' **Edge cases handled:**
#' - Peak with apex height <= 0: boundaries unchanged
#' - Fewer than 5 data points within boundaries: boundaries unchanged
#' - No points above threshold: boundaries unchanged
#' - Refined boundaries are guaranteed to never be wider than originals
#'
#' # Tidying
#'
#' When you [`tidy()`][recipes::tidy.recipe()] this step, a tibble with columns
#' `method`, `cutoff`, `peaks_col`, `measures_col`, and `id` is returned.
#'
#' @seealso [step_sec_peaks_detect()] for peak detection,
#'   [step_sec_baseline()] for baseline correction.
#'
#' @family sec-chromatography
#' @export
#'
#' @examples
#' \dontrun{
#' library(recipes)
#' library(measure)
#'
#' # Refine peaks with default 0.5% height fraction
#' rec <- recipe(~., data = sec_data) |>
#'   step_measure_input_long(ri_signal, location = vars(elution_volume)) |>
#'   step_sec_baseline() |>
#'   step_sec_peaks_detect() |>
#'   step_sec_peaks_refine() |>
#'   prep()
#'
#' # More aggressive trimming at 1%
#' rec <- recipe(~., data = sec_data) |>
#'   step_measure_input_long(ri_signal, location = vars(elution_volume)) |>
#'   step_sec_baseline() |>
#'   step_sec_peaks_detect() |>
#'   step_sec_peaks_refine(cutoff = 0.01) |>
#'   prep()
#'
#' # Looser boundaries at 0.1%
#' rec <- recipe(~., data = sec_data) |>
#'   step_measure_input_long(ri_signal, location = vars(elution_volume)) |>
#'   step_sec_baseline() |>
#'   step_sec_peaks_detect() |>
#'   step_sec_peaks_refine(cutoff = 0.001) |>
#'   prep()
#' }
step_sec_peaks_refine <- function(
	recipe,
	method = "height_fraction",
	cutoff = 0.005,
	peaks_col = ".peaks",
	measures_col = ".measures",
	role = NA,
	trained = FALSE,
	skip = FALSE,
	id = recipes::rand_id("sec_peaks_refine")
) {
	method <- rlang::arg_match(method, "height_fraction")

	if (
		!is.numeric(cutoff) || length(cutoff) != 1 || cutoff <= 0 || cutoff >= 1
	) {
		cli::cli_abort(
			"{.arg cutoff} must be a number between 0 and 1 (exclusive)."
		)
	}

	recipes::add_step(
		recipe,
		step_sec_peaks_refine_new(
			method = method,
			cutoff = cutoff,
			peaks_col = peaks_col,
			measures_col = measures_col,
			role = role,
			trained = trained,
			skip = skip,
			id = id
		)
	)
}

step_sec_peaks_refine_new <- function(
	method,
	cutoff,
	peaks_col,
	measures_col,
	role,
	trained,
	skip,
	id
) {
	recipes::step(
		subclass = "sec_peaks_refine",
		method = method,
		cutoff = cutoff,
		peaks_col = peaks_col,
		measures_col = measures_col,
		role = role,
		trained = trained,
		skip = skip,
		id = id
	)
}

#' @export
prep.step_sec_peaks_refine <- function(x, training, info = NULL, ...) {
	# Verify peaks column exists
	if (!x$peaks_col %in% names(training)) {
		cli::cli_abort(
			c(
				"Column {.val {x$peaks_col}} not found.",
				"i" = "Run {.fn step_sec_peaks_detect} first to detect peaks."
			)
		)
	}

	# Resolve measures column
	measures_col <- x$measures_col
	if (!measures_col %in% names(training)) {
		measure_cols <- find_measure_cols(training)
		if (length(measure_cols) > 0) {
			measures_col <- measure_cols[1]
			cli::cli_inform(
				c(
					"i" = "Column {.val {x$measures_col}} not found.",
					"i" = "Using {.val {measures_col}} instead."
				)
			)
		} else {
			cli::cli_abort("No measure columns found in data.")
		}
	}

	step_sec_peaks_refine_new(
		method = x$method,
		cutoff = x$cutoff,
		peaks_col = x$peaks_col,
		measures_col = measures_col,
		role = x$role,
		trained = TRUE,
		skip = x$skip,
		id = x$id
	)
}

#' @export
bake.step_sec_peaks_refine <- function(object, new_data, ...) {
	peaks_col <- object$peaks_col
	measures_col <- object$measures_col
	cutoff <- object$cutoff
	n_rows <- nrow(new_data)

	refined_peaks <- vector("list", n_rows)

	for (i in seq_len(n_rows)) {
		peaks <- new_data[[peaks_col]][[i]]
		signal <- new_data[[measures_col]][[i]]

		if (nrow(peaks) == 0) {
			refined_peaks[[i]] <- peaks
			next
		}

		refined_peaks[[i]] <- .refine_peaks_height_fraction(
			peaks,
			signal$location,
			signal$value,
			cutoff
		)
	}

	new_data[[peaks_col]] <- measure:::new_peaks_list(refined_peaks)
	tibble::as_tibble(new_data)
}

#' @export
print.step_sec_peaks_refine <- function(
	x,
	width = max(20, options()$width - 30),
	...
) {
	cat(
		glue::glue("SEC peak boundary refinement ({x$method}, cutoff={x$cutoff})")
	)
	cat("\n")
	invisible(x)
}

#' @rdname tidy.step_sec
#' @export
#' @keywords internal
tidy.step_sec_peaks_refine <- function(x, ...) {
	tibble::tibble(
		method = x$method,
		cutoff = x$cutoff,
		peaks_col = x$peaks_col,
		measures_col = x$measures_col %||% NA_character_,
		id = x$id
	)
}

#' @rdname required_pkgs.step_sec
#' @export
#' @keywords internal
required_pkgs.step_sec_peaks_refine <- function(x, ...) {
	c("measure.sec", "measure")
}

# ==============================================================================
# Height fraction refinement
# ==============================================================================

#' Refine peak boundaries using height fraction method
#'
#' For each peak, walks inward from the current boundaries until the signal
#' exceeds `cutoff * apex_height`. Preserves original boundaries as
#' `original_left_base` and `original_right_base`.
#'
#' @param peaks A `peaks_tbl` from peak detection.
#' @param location Numeric vector of x-axis values (elution volume).
#' @param value Numeric vector of y-axis values (baseline-corrected signal).
#' @param cutoff Fraction of apex height for boundary threshold.
#'
#' @return Modified `peaks_tbl` with refined boundaries and original boundary
#'   columns.
#' @noRd
.refine_peaks_height_fraction <- function(peaks, location, value, cutoff) {
	n_peaks <- nrow(peaks)

	# Store original boundaries
	orig_left <- peaks$left_base
	orig_right <- peaks$right_base

	for (j in seq_len(n_peaks)) {
		# Find indices within the current peak boundaries
		in_region <- which(location >= orig_left[j] & location <= orig_right[j])

		# Skip if fewer than 5 data points
		if (length(in_region) < 5) next

		region_values <- value[in_region]
		region_locations <- location[in_region]

		# Find apex height within this region
		apex_height <- max(region_values, na.rm = TRUE)

		# Skip if apex height is non-positive
		if (is.na(apex_height) || apex_height <= 0) next

		threshold <- cutoff * apex_height

		# Find points above threshold
		above <- region_values >= threshold
		if (!any(above)) next

		# Walk inward from left: first point >= threshold
		above_idx <- which(above)
		new_left <- region_locations[above_idx[1]]

		# Walk inward from right: last point >= threshold
		new_right <- region_locations[above_idx[length(above_idx)]]

		# Enforce that refined boundaries never widen beyond originals
		peaks$left_base[j] <- max(new_left, orig_left[j])
		peaks$right_base[j] <- min(new_right, orig_right[j])
	}

	# Add original boundary columns
	peaks$original_left_base <- orig_left
	peaks$original_right_base <- orig_right

	peaks
}
