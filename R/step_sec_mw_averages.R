# ==============================================================================
# step_sec_mw_averages
#
# Calculate molecular weight averages (Mn, Mw, Mz, Mp, dispersity)
# ==============================================================================

#' Calculate Molecular Weight Averages for SEC/GPC
#'
#' `step_sec_mw_averages()` creates a *specification* of a recipe step that
#' calculates molecular weight averages from size exclusion chromatography data.
#'
#' @param recipe A recipe object.
#' @param measures An optional character vector of measure column names.
#' @param calibration Calibration method for converting x-axis to log(MW).
#'   Can be:
#'   - `NULL` (default): Assumes x-axis is already log10(MW)
#'   - A numeric vector of length 2: Linear calibration `c(slope, intercept)`
#'     where `log10(MW) = slope * x + intercept`
#'   - `"auto"`: Estimate from data range (assumes typical polymer range)
#' @param integration_range Optional numeric vector `c(min, max)` specifying
#'   the x-axis range for integration. If `NULL`, uses full range.
#' @param output_cols Character vector of metrics to calculate. Default
#'   includes all: `c("mn", "mw", "mz", "mp", "dispersity")`.
#' @param include_uncertainty Logical. If `TRUE`, calculates and outputs
#'   uncertainty estimates for MW averages. Requires `calibration_error` to
#'   be specified. Default is `FALSE`.
#' @param calibration_error Calibration error (RMSE) in log10(MW) units for
#'   uncertainty propagation. Required when `include_uncertainty = TRUE`.
#'   Can be obtained from `tidy()` output of `step_sec_conventional_cal()`.
#' @param prefix Prefix for output column names. Default is `"mw_"`.
#' @param role Role for generated columns. Default is `"predictor"`.
#' @param trained Logical indicating if the step has been trained.
#' @param skip Logical. Should the step be skipped when baking?
#' @param id Unique step identifier.
#'
#' @return An updated recipe with the new step added.
#'
#' @details
#' This step calculates standard molecular weight averages from SEC/GPC data:
#'
#' | Metric | Formula | Description |
#' |--------|---------|-------------|
#' | Mn | sum(w) / sum(w/M) | Number-average molecular weight |
#' | Mw | sum(wM) / sum(w) | Weight-average molecular weight |
#' | Mz | sum(wM^2) / sum(wM) | Z-average molecular weight |
#' | Mp | M at peak maximum | Peak molecular weight |
#' | D | Mw/Mn | Dispersity (polydispersity index) |
#'
#' The detector signal is assumed to be proportional to weight concentration.
#' For RI detection, this is typically valid. For UV detection, response factors
#' may need to be applied first.
#'
#' **Uncertainty Propagation:**
#'
#' When `include_uncertainty = TRUE`, the step calculates uncertainty estimates
#' based on calibration error propagation. The uncertainties account for:
#' - Calibration curve fit error (RMSE in log10 MW units)
#' - MW distribution width effects on different averages
#'
#' The propagation follows:
#' - Mn uncertainty is enhanced for wide distributions (most sensitive to low MW)
#' - Mw uncertainty equals the relative calibration error
#' - Mz uncertainty is enhanced for high MW sensitivity
#' - Dispersity uncertainty from error propagation of Mw/Mn
#'
#' **Prerequisites:**
#' - Data should be baseline corrected
#' - X-axis should represent retention time/volume or log(MW)
#' - Integration limits should exclude solvent peaks
#'
#' @family sec-chromatography
#' @export
#'
#' @examples
#' \dontrun{
#' library(recipes)
#' library(measure)
#'
#' # Assuming x-axis is already calibrated to log10(MW)
#' rec <- recipe(~., data = sec_triple_detect) |>
#'   step_measure_input_wide(starts_with("signal_")) |>
#'   step_sec_baseline() |>
#'   step_sec_mw_averages() |>
#'   prep()
#'
#' # With uncertainty propagation (calibration_error from tidy() of calibration step)
#' rec_with_unc <- recipe(~., data = sec_triple_detect) |>
#'   step_measure_input_wide(starts_with("signal_")) |>
#'   step_sec_baseline() |>
#'   step_sec_mw_averages(
#'     include_uncertainty = TRUE,
#'     calibration_error = 0.02  # RMSE in log10(MW)
#'   ) |>
#'   prep()
#' }
step_sec_mw_averages <- function(
	recipe,
	measures = NULL,
	calibration = NULL,
	integration_range = NULL,
	output_cols = c("mn", "mw", "mz", "mp", "dispersity"),
	include_uncertainty = FALSE,
	calibration_error = NULL,
	prefix = "mw_",
	role = "predictor",
	trained = FALSE,
	skip = FALSE,
	id = recipes::rand_id("sec_mw_averages")
) {
	valid_cols <- c("mn", "mw", "mz", "mp", "dispersity")
	if (!all(output_cols %in% valid_cols)) {
		invalid <- setdiff(output_cols, valid_cols)
		cli::cli_abort(
			"Invalid output columns: {.val {invalid}}. Must be one of: {.val {valid_cols}}"
		)
	}

	if (!is.null(integration_range)) {
		if (!is.numeric(integration_range) || length(integration_range) != 2) {
			cli::cli_abort(
				"{.arg integration_range} must be a numeric vector of length 2."
			)
		}
	}

	# Validate uncertainty parameters
	if (include_uncertainty && is.null(calibration_error)) {
		cli::cli_abort(
			c(
				"{.arg calibration_error} is required when {.arg include_uncertainty} is TRUE.",
				"i" = "Get calibration error (RMSE) from {.fn tidy} output of {.fn step_sec_conventional_cal}."
			)
		)
	}

	if (!is.null(calibration_error)) {
		if (!is.numeric(calibration_error) || length(calibration_error) != 1) {
			cli::cli_abort(
				"{.arg calibration_error} must be a single numeric value (RMSE in log10 MW units)."
			)
		}
		if (calibration_error <= 0) {
			cli::cli_abort(
				"{.arg calibration_error} must be positive."
			)
		}
	}

	recipes::add_step(
		recipe,
		step_sec_mw_averages_new(
			measures = measures,
			calibration = calibration,
			integration_range = integration_range,
			output_cols = output_cols,
			include_uncertainty = include_uncertainty,
			calibration_error = calibration_error,
			prefix = prefix,
			role = role,
			trained = trained,
			skip = skip,
			id = id
		)
	)
}

step_sec_mw_averages_new <- function(
	measures,
	calibration,
	integration_range,
	output_cols,
	include_uncertainty,
	calibration_error,
	prefix,
	role,
	trained,
	skip,
	id
) {
	recipes::step(
		subclass = "sec_mw_averages",
		measures = measures,
		calibration = calibration,
		integration_range = integration_range,
		output_cols = output_cols,
		include_uncertainty = include_uncertainty,
		calibration_error = calibration_error,
		prefix = prefix,
		role = role,
		trained = trained,
		skip = skip,
		id = id
	)
}

#' @export
prep.step_sec_mw_averages <- function(x, training, info = NULL, ...) {
	check_for_measure(training)

	if (is.null(x$measures)) {
		measure_cols <- find_measure_cols(training)
	} else {
		measure_cols <- x$measures
	}

	step_sec_mw_averages_new(
		measures = measure_cols,
		calibration = x$calibration,
		integration_range = x$integration_range,
		output_cols = x$output_cols,
		include_uncertainty = x$include_uncertainty,
		calibration_error = x$calibration_error,
		prefix = x$prefix,
		role = x$role,
		trained = TRUE,
		skip = x$skip,
		id = x$id
	)
}

#' Calculate MW uncertainties from calibration error
#'
#' Propagates calibration curve uncertainty to MW average uncertainties.
#' Based on error propagation theory and SEC-specific considerations.
#'
#' @param mw MW values from chromatogram
#' @param w Weight values (signal intensity)
#' @param mn Pre-calculated Mn value
#' @param mw_avg Pre-calculated Mw value
#' @param mz Pre-calculated Mz value
#' @param calibration_error RMSE in log10(MW) units
#'
#' @return Named list with uncertainties for mn, mw, mz, dispersity
#' @noRd
.calc_mw_uncertainties <- function(mw, w, mn, mw_avg, mz, calibration_error) {
	# Convert calibration RMSE from log10 to relative uncertainty

	# sigma(MW)/MW = ln(10) * sigma(log10 MW) ~= 2.303 * sigma(log10 MW)
	relative_cal_error <- log(10) * calibration_error

	# Calculate the MW distribution width (in log units)
	log_mw_std <- stats::sd(log10(mw))

	# If distribution width is NA or zero, use a default
	if (is.na(log_mw_std) || log_mw_std == 0) {
		log_mw_std <- 0.3 # Typical for moderate dispersity
	}

	# Mn uncertainty: Enhanced by distribution width
	# Wider distributions have more uncertainty in Mn
	mn_rel_uncertainty <- relative_cal_error * sqrt(1 + 2 * log_mw_std^2)
	mn_uncertainty <- mn * mn_rel_uncertainty

	# Mw uncertainty: Standard calibration uncertainty
	mw_rel_uncertainty <- relative_cal_error
	mw_uncertainty <- mw_avg * mw_rel_uncertainty

	# Mz uncertainty: Enhanced for high MW sensitivity
	mz_rel_uncertainty <- relative_cal_error * sqrt(1 + 4 * log_mw_std^2)
	mz_uncertainty <- mz * mz_rel_uncertainty

	# PDI/Dispersity uncertainty from error propagation
	dispersity <- mw_avg / mn
	dispersity_rel_uncertainty <- sqrt(
		mw_rel_uncertainty^2 + mn_rel_uncertainty^2
	)
	dispersity_uncertainty <- dispersity * dispersity_rel_uncertainty

	list(
		mn_uncertainty = signif(mn_uncertainty, 2),
		mw_uncertainty = signif(mw_uncertainty, 2),
		mz_uncertainty = signif(mz_uncertainty, 2),
		dispersity_uncertainty = signif(dispersity_uncertainty, 2)
	)
}

#' Calculate MW averages from a single chromatogram
#' @noRd
.calc_mw_averages <- function(
	location,
	value,
	calibration,
	integration_range,
	output_cols,
	include_uncertainty = FALSE,
	calibration_error = NULL
) {
	# Build full list of output column names (including uncertainties if requested)
	all_output_cols <- output_cols
	if (include_uncertainty) {
		all_output_cols <- c(
			output_cols,
			"mn_uncertainty",
			"mw_uncertainty",
			"mz_uncertainty",
			"dispersity_uncertainty"
		)
	}

	# Apply integration range
	if (!is.null(integration_range)) {
		idx <- location >= integration_range[1] & location <= integration_range[2]
		location <- location[idx]
		value <- value[idx]
	}

	if (length(location) < 2) {
		result <- stats::setNames(
			rep(NA_real_, length(all_output_cols)),
			all_output_cols
		)
		return(result)
	}

	# Convert x-axis to log10(MW) if calibration provided
	if (is.null(calibration)) {
		log_mw <- location
	} else if (is.numeric(calibration) && length(calibration) == 2) {
		# Linear calibration: log10(MW) = slope * x + intercept
		log_mw <- calibration[1] * location + calibration[2]
	} else if (identical(calibration, "auto")) {
		# Auto-calibration: assume typical polymer range (1e2 to 1e7)
		log_mw <- seq(7, 2, length.out = length(location))
	} else {
		log_mw <- location
	}

	# Convert to MW
	mw <- 10^log_mw

	# Weight is proportional to signal (assuming RI detector)
	# Ensure non-negative
	w <- pmax(value, 0)

	# Remove zero weights
	valid <- w > 0
	if (sum(valid) < 2) {
		result <- stats::setNames(
			rep(NA_real_, length(all_output_cols)),
			all_output_cols
		)
		return(result)
	}

	mw <- mw[valid]
	w <- w[valid]

	result <- numeric(length(output_cols))
	names(result) <- output_cols

	# Calculate moments
	sum_w <- sum(w)

	if ("mn" %in% output_cols) {
		# Mn = sum(w) / sum(w/M)
		result["mn"] <- sum_w / sum(w / mw)
	}

	if ("mw" %in% output_cols) {
		# Mw = sum(wM) / sum(w)
		result["mw"] <- sum(w * mw) / sum_w
	}

	if ("mz" %in% output_cols) {
		# Mz = sum(wM^2) / sum(wM)
		result["mz"] <- sum(w * mw^2) / sum(w * mw)
	}

	if ("mp" %in% output_cols) {
		# Mp = MW at peak maximum
		peak_idx <- which.max(w)
		result["mp"] <- mw[peak_idx]
	}

	if ("dispersity" %in% output_cols) {
		# D = Mw / Mn
		if ("mw" %in% names(result) && "mn" %in% names(result)) {
			result["dispersity"] <- result["mw"] / result["mn"]
		} else {
			mw_val <- sum(w * mw) / sum_w
			mn_val <- sum_w / sum(w / mw)
			result["dispersity"] <- mw_val / mn_val
		}
	}

	# Add uncertainty calculations if requested
	if (include_uncertainty && !is.null(calibration_error)) {
		# Calculate MW values if not already done (for uncertainty calculation)
		mn_val <- if ("mn" %in% names(result)) result["mn"] else sum_w / sum(w / mw)
		mw_val <- if ("mw" %in% names(result)) result["mw"] else sum(w * mw) / sum_w
		mz_val <- if ("mz" %in% names(result)) {
			result["mz"]
		} else {
			sum(w * mw^2) / sum(w * mw)
		}

		uncertainties <- .calc_mw_uncertainties(
			mw = mw,
			w = w,
			mn = mn_val,
			mw_avg = mw_val,
			mz = mz_val,
			calibration_error = calibration_error
		)

		# Add uncertainties to result
		result["mn_uncertainty"] <- uncertainties$mn_uncertainty
		result["mw_uncertainty"] <- uncertainties$mw_uncertainty
		result["mz_uncertainty"] <- uncertainties$mz_uncertainty
		result["dispersity_uncertainty"] <- uncertainties$dispersity_uncertainty
	}

	result
}

#' @export
bake.step_sec_mw_averages <- function(object, new_data, ...) {
	calibration <- object$calibration
	integration_range <- object$integration_range
	output_cols <- object$output_cols
	include_uncertainty <- object$include_uncertainty
	calibration_error <- object$calibration_error
	prefix <- object$prefix

	# Calculate MW averages for each sample
	all_results <- purrr::map(new_data[[object$measures[1]]], function(m) {
		.calc_mw_averages(
			m$location,
			m$value,
			calibration,
			integration_range,
			output_cols,
			include_uncertainty = include_uncertainty,
			calibration_error = calibration_error
		)
	})

	# Convert to data frame
	result_df <- do.call(rbind, all_results)
	result_df <- tibble::as_tibble(result_df)

	# Add prefix to column names
	names(result_df) <- paste0(prefix, names(result_df))

	# Bind to original data
	new_data <- dplyr::bind_cols(new_data, result_df)

	tibble::as_tibble(new_data)
}

#' @export
print.step_sec_mw_averages <- function(
	x,
	width = max(20, options()$width - 30),
	...
) {
	cols <- paste(x$output_cols, collapse = ", ")
	title <- paste0("SEC MW averages (", cols, ")")
	if (x$include_uncertainty) {
		title <- paste0(title, " + uncertainties")
	}
	if (x$trained) {
		cat(title, " on <internal measurements>", sep = "")
	} else {
		cat(title)
	}
	cat("\n")
	invisible(x)
}

#' @rdname tidy.step_sec
#' @export
#' @keywords internal
tidy.step_sec_mw_averages <- function(x, ...) {
	tibble::tibble(
		output_cols = list(x$output_cols),
		include_uncertainty = x$include_uncertainty,
		calibration_error = x$calibration_error %||% NA_real_,
		prefix = x$prefix,
		id = x$id
	)
}

#' @rdname required_pkgs.step_sec
#' @export
#' @keywords internal
required_pkgs.step_sec_mw_averages <- function(x, ...) {
	c("measure.sec", "measure")
}
