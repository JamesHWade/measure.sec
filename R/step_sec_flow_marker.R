# ==============================================================================
# step_sec_flow_marker
#
# Flow marker detection and correction for SEC/GPC data
# ==============================================================================

#' Flow Marker Correction for SEC/GPC
#'
#' `step_sec_flow_marker()` creates a *specification* of a recipe step that
#' detects a flow marker peak (typically toluene or other small molecule) and
#' applies a linear correction to align retention times/volumes across runs.
#'
#' @param recipe A recipe object.
#' @param measures Character vector of measure column names to correct.
#'   If `NULL`, uses all measure columns.
#' @param marker_range Numeric vector of length 2 specifying the expected
#'   elution range `c(min, max)` where the flow marker peak is expected.
#'   Required unless `auto_detect = TRUE` with `expected_volume` specified.
#' @param target_volume The target elution volume to align the flow marker to.
#'   If `NULL` (default), uses the first value of `marker_range`.
#' @param auto_detect Logical. If `TRUE`, automatically detect the flow marker
#'   peak within `marker_range`. If `FALSE`, uses the peak maximum within range.
#'   Default is `TRUE`.
#' @param min_peak_height Minimum peak height (in signal units) for a valid

#'   flow marker detection. Peaks below this threshold are ignored.
#'   Default is `NULL` (no minimum).
#' @param store_correction Logical. If `TRUE`, stores the correction factor
#'   as a column named `flow_marker_correction`. Default is `TRUE`.
#' @param role Role for generated columns.
#' @param trained Logical indicating if the step has been trained.
#' @param skip Logical. Should the step be skipped when baking?
#' @param id Unique step identifier.
#'
#' @return An updated recipe with the new step added.
#'
#' @details
#' Flow marker correction compensates for small variations in flow rate between
#' runs. A flow marker is a small molecule (often toluene) that elutes near the
#' total permeation volume and provides a reference point for alignment.
#'
#' **Correction Algorithm:**
#'
#' 1. Find the flow marker peak maximum within `marker_range`
#' 2. Calculate the shift: `correction = observed_volume - target_volume`
#' 3. Apply linear correction: `corrected_volume = original_volume - correction`
#'
#' **Auto-Detection:**
#'
#' When `auto_detect = TRUE`, the step uses second-derivative analysis to find
#' the sharpest peak in the specified range, which is typically the flow marker.
#' This is more robust than simply finding the maximum signal.
#'
#' **Prerequisites:**
#' - Should be applied before calibration and MW calculations
#' - Best applied after baseline correction
#'
#' @family sec-preprocessing
#' @export
#'
#' @examples
#' \dontrun{
#' library(recipes)
#' library(measure)
#'
#' # Apply flow marker correction with known range
#' rec <- recipe(~., data = sec_data) |>
#'   step_measure_input_long(ri, location = vars(elution_volume)) |>
#'   step_sec_baseline() |>
#'   step_sec_flow_marker(
#'     marker_range = c(18, 20),
#'     target_volume = 18.5
#'   ) |>
#'   step_sec_conventional_cal(standards = ps_standards) |>
#'   prep()
#'
#' # Auto-detect flow marker with default target (first value of range)
#' rec <- recipe(~., data = sec_data) |>
#'   step_measure_input_long(ri, location = vars(elution_volume)) |>
#'   step_sec_flow_marker(marker_range = c(18, 20)) |>
#'   prep()
#' }
step_sec_flow_marker <- function(
	recipe,
	measures = NULL,
	marker_range = NULL,
	target_volume = NULL,
	auto_detect = TRUE,
	min_peak_height = NULL,
	store_correction = TRUE,
	role = NA,
	trained = FALSE,
	skip = FALSE,
	id = recipes::rand_id("sec_flow_marker")
) {
	# Validate marker_range

	if (is.null(marker_range)) {
		cli::cli_abort(
			c(
				"{.arg marker_range} is required.",
				"i" = "Specify the expected elution range for the flow marker, e.g., {.code c(18, 20)}."
			)
		)
	}

	if (!is.numeric(marker_range) || length(marker_range) != 2) {
		cli::cli_abort(
			"{.arg marker_range} must be a numeric vector of length 2."
		)
	}

	if (marker_range[1] >= marker_range[2]) {
		cli::cli_abort(
			"{.arg marker_range} must be in order: {.code c(min, max)}."
		)
	}

	# Set default target volume
	if (is.null(target_volume)) {
		target_volume <- marker_range[1]
	}

	if (!is.numeric(target_volume) || length(target_volume) != 1) {
		cli::cli_abort(
			"{.arg target_volume} must be a single numeric value."
		)
	}

	if (!is.null(min_peak_height)) {
		if (!is.numeric(min_peak_height) || length(min_peak_height) != 1) {
			cli::cli_abort(
				"{.arg min_peak_height} must be a single numeric value."
			)
		}
	}

	recipes::add_step(
		recipe,
		step_sec_flow_marker_new(
			measures = measures,
			marker_range = marker_range,
			target_volume = target_volume,
			auto_detect = auto_detect,
			min_peak_height = min_peak_height,
			store_correction = store_correction,
			role = role,
			trained = trained,
			skip = skip,
			id = id
		)
	)
}

step_sec_flow_marker_new <- function(
	measures,
	marker_range,
	target_volume,
	auto_detect,
	min_peak_height,
	store_correction,
	role,
	trained,
	skip,
	id
) {
	recipes::step(
		subclass = "sec_flow_marker",
		measures = measures,
		marker_range = marker_range,
		target_volume = target_volume,
		auto_detect = auto_detect,
		min_peak_height = min_peak_height,
		store_correction = store_correction,
		role = role,
		trained = trained,
		skip = skip,
		id = id
	)
}

#' @export
prep.step_sec_flow_marker <- function(x, training, info = NULL, ...) {
	check_for_measure(training)

	if (is.null(x$measures)) {
		measures <- find_measure_cols(training)
	} else {
		measures <- x$measures
	}

	if (length(measures) == 0) {
		cli::cli_abort("No measure columns found in data.")
	}

	step_sec_flow_marker_new(
		measures = measures,
		marker_range = x$marker_range,
		target_volume = x$target_volume,
		auto_detect = x$auto_detect,
		min_peak_height = x$min_peak_height,
		store_correction = x$store_correction,
		role = x$role,
		trained = TRUE,
		skip = x$skip,
		id = x$id
	)
}

#' Detect flow marker peak in a chromatogram
#'
#' @param location Numeric vector of elution positions
#' @param value Numeric vector of signal values
#' @param marker_range Numeric vector c(min, max) for search range
#' @param auto_detect Logical, use second derivative detection
#' @param min_peak_height Minimum peak height threshold
#'
#' @return List with flow_marker_volume, correction, and confidence
#' @noRd
.detect_flow_marker <- function(
	location,
	value,
	marker_range,
	auto_detect = TRUE,
	min_peak_height = NULL
) {
	# Filter to marker range
	in_range <- location >= marker_range[1] & location <= marker_range[2]
	range_location <- location[in_range]
	range_value <- value[in_range]

	if (length(range_location) < 3) {
		return(
			list(
				flow_marker_volume = NA_real_,
				peak_height = NA_real_,
				confidence = 0
			)
		)
	}

	# Apply minimum height filter if specified
	if (!is.null(min_peak_height)) {
		valid <- abs(range_value) >= min_peak_height
		if (sum(valid) < 3) {
			return(
				list(
					flow_marker_volume = NA_real_,
					peak_height = NA_real_,
					confidence = 0
				)
			)
		}
		# Filter to only valid points that meet threshold
		range_location <- range_location[valid]
		range_value <- range_value[valid]
	}

	if (auto_detect && length(range_value) >= 10) {
		# Use second derivative to find sharpest peak
		# This is more robust for finding the flow marker among other features
		d2_signal <- diff(diff(range_value))

		if (length(d2_signal) > 0 && !all(is.na(d2_signal))) {
			# Find the position with largest negative second derivative (peak)
			peak_idx <- which.min(d2_signal) + 1

			# Ensure valid index
			peak_idx <- max(1, min(peak_idx, length(range_location)))

			flow_marker_volume <- range_location[peak_idx]
			peak_height <- range_value[peak_idx]

			# Calculate confidence based on sharpness
			d2_range <- max(abs(d2_signal), na.rm = TRUE)
			if (d2_range > 0) {
				confidence <- abs(d2_signal[max(1, peak_idx - 1)]) / d2_range
			} else {
				confidence <- 0.5
			}
		} else {
			# Fallback to simple maximum
			peak_idx <- which.max(abs(range_value))
			flow_marker_volume <- range_location[peak_idx]
			peak_height <- range_value[peak_idx]
			confidence <- 0.5
		}
	} else {
		# Simple maximum detection
		peak_idx <- which.max(abs(range_value))
		flow_marker_volume <- range_location[peak_idx]
		peak_height <- range_value[peak_idx]
		confidence <- 0.5
	}

	list(
		flow_marker_volume = flow_marker_volume,
		peak_height = peak_height,
		confidence = confidence
	)
}

#' @export
bake.step_sec_flow_marker <- function(object, new_data, ...) {
	measures <- object$measures
	marker_range <- object$marker_range
	target_volume <- object$target_volume
	auto_detect <- object$auto_detect
	min_peak_height <- object$min_peak_height
	store_correction <- object$store_correction

	# Use the first measure column for flow marker detection
	detection_col <- measures[1]

	# Calculate correction for each sample and apply to all measure columns
	corrections <- numeric(nrow(new_data))

	for (i in seq_len(nrow(new_data))) {
		m <- new_data[[detection_col]][[i]]

		detection <- .detect_flow_marker(
			location = m$location,
			value = m$value,
			marker_range = marker_range,
			auto_detect = auto_detect,
			min_peak_height = min_peak_height
		)

		if (!is.na(detection$flow_marker_volume)) {
			corrections[i] <- detection$flow_marker_volume - target_volume
		} else {
			corrections[i] <- 0
			cli::cli_warn(
				"Row {i}: Flow marker not detected in range [{marker_range[1]}, {marker_range[2]}]. No correction applied."
			)
		}
	}

	# Apply correction to all measure columns
	for (col in measures) {
		new_data[[col]] <- purrr::map2(
			new_data[[col]],
			corrections,
			function(m, corr) {
				new_measure_tbl(
					location = m$location - corr,
					value = m$value
				)
			}
		)
	}

	# Store correction if requested
	if (store_correction) {
		new_data$flow_marker_correction <- corrections
	}

	tibble::as_tibble(new_data)
}

#' @export
print.step_sec_flow_marker <- function(
	x,
	width = max(20, options()$width - 30),
	...
) {
	title <- glue::glue(
		"SEC flow marker correction (range: [{x$marker_range[1]}, {x$marker_range[2]}], target: {x$target_volume})"
	)

	if (x$trained) {
		cat(title, "\n", sep = "")
	} else {
		cat(title, " [untrained]\n", sep = "")
	}

	invisible(x)
}

#' @rdname tidy.step_sec
#' @export
#' @keywords internal
tidy.step_sec_flow_marker <- function(x, ...) {
	tibble::tibble(
		marker_range_min = x$marker_range[1],
		marker_range_max = x$marker_range[2],
		target_volume = x$target_volume,
		auto_detect = x$auto_detect,
		id = x$id
	)
}

#' @rdname required_pkgs.step_sec
#' @export
#' @keywords internal
required_pkgs.step_sec_flow_marker <- function(x, ...) {
	c("measure.sec", "measure")
}
