# ==============================================================================
# step_sec_aggregates
#
# Protein SEC aggregate and fragment quantitation
# ==============================================================================

#' Quantify Protein Aggregates and Fragments in SEC
#'
#' `step_sec_aggregates()` creates a *specification* of a recipe step that
#' quantifies high molecular weight species (HMWS/aggregates), monomers, and
#' low molecular weight species (LMWS/fragments) from protein SEC chromatograms.
#'
#' @param recipe A recipe object.
#' @param measures Character vector of measure columns to analyze. If `NULL`,
#'   analyzes all measure columns.
#' @param monomer_start Start of the monomer peak region (in location units,
#'   typically minutes). If `NULL`, automatically determined.
#' @param monomer_end End of the monomer peak region. If `NULL`, automatically
#'   determined.
#' @param method Method for peak boundary detection when `monomer_start` or
#'   `monomer_end` is `NULL`:
#'   \itemize{
#'     \item `"tallest"` (default): Uses the tallest peak as monomer
#'     \item `"manual"`: Requires explicit boundaries
#'   }
#' @param hmws_threshold Minimum fraction of monomer height to consider as
#'   HMWS signal. Default is 0.001 (0.1%). Below this, signal is considered
#'   baseline.
#' @param include_main_peak Logical. Include the main peak boundaries in output?
#'   Default is TRUE.
#' @param output_prefix Prefix for output columns. Default is `"purity_"`.
#'   Creates columns: `{prefix}hmws`, `{prefix}monomer`, `{prefix}lmws`.
#' @param role Role for generated columns.
#' @param trained Logical indicating if the step has been trained.
#' @param skip Logical. Should the step be skipped when baking?
#' @param id Unique step identifier.
#'
#' @return An updated recipe with new columns containing aggregate percentages:
#' \describe{
#'   \item{purity_hmws}{Percent high molecular weight species (aggregates)}
#'   \item{purity_monomer}{Percent monomer (main peak)}
#'   \item{purity_lmws}{Percent low molecular weight species (fragments)}
#'   \item{purity_main_start}{Start of main peak region (if include_main_peak)}
#'   \item{purity_main_end}{End of main peak region (if include_main_peak)}
#' }
#'
#' @details
#' Aggregate analysis is critical for biopharmaceutical characterization:
#'
#' \itemize{
#'   \item HMWS (High Molecular Weight Species): Elute before the monomer peak.
#'     Includes dimers, trimers, and higher-order aggregates.
#'   \item Monomer: The main therapeutic protein peak.
#'   \item LMWS (Low Molecular Weight Species): Elute after the monomer peak.
#'     Includes fragments, clips, and degradation products.
#' }
#'
#' **Calculation Method:**
#' \deqn{\% HMWS = \frac{A_{HMWS}}{A_{total}} \times 100}
#' \deqn{\% Monomer = \frac{A_{monomer}}{A_{total}} \times 100}
#' \deqn{\% LMWS = \frac{A_{LMWS}}{A_{total}} \times 100}
#'
#' where A represents integrated peak areas.
#'
#' **Regulatory Importance:**
#' \itemize{
#'   \item ICH Q6B requires aggregate content specification
#'   \item USP <129> provides guidance on aggregate testing
#'   \item Typical acceptance: HMWS < 5%, Monomer > 95%
#' }
#'
#' @note
#' For accurate results:
#' \itemize{
#'   \item Baseline correct the chromatogram first
#'   \item Ensure proper column resolution (especially for dimer separation)
#'   \item Use UV detection at 280 nm (or 214 nm for higher sensitivity)
#' }
#'
#' @family sec-protein
#' @export
#'
#' @examples
#' \dontrun{
#' library(recipes)
#' library(measure)
#'
#' # Analyze mAb aggregate content
#' rec <- recipe(~., data = mab_data) |>
#'   step_measure_input_long(uv280, location = vars(elution_time), col_name = "uv") |>
#'   step_sec_baseline() |>
#'   step_sec_aggregates(
#'     monomer_start = 8.5,
#'     monomer_end = 10.5
#'   ) |>
#'   prep()
#'
#' # Automatic peak detection
#' rec <- recipe(~., data = mab_data) |>
#'   step_measure_input_long(uv280, location = vars(elution_time), col_name = "uv") |>
#'   step_sec_baseline() |>
#'   step_sec_aggregates(method = "tallest") |>
#'   prep()
#' }
step_sec_aggregates <- function(
	recipe,
	measures = NULL,
	monomer_start = NULL,
	monomer_end = NULL,
	method = c("tallest", "manual"),
	hmws_threshold = 0.001,
	include_main_peak = TRUE,
	output_prefix = "purity_",
	role = NA,
	trained = FALSE,
	skip = FALSE,
	id = recipes::rand_id("sec_aggregates")
) {
	method <- match.arg(method)

	# Validate manual method requires boundaries
	if (method == "manual" && (is.null(monomer_start) || is.null(monomer_end))) {
		cli::cli_abort(
			"Method {.val manual} requires both {.arg monomer_start} and {.arg monomer_end}."
		)
	}

	if (!is.numeric(hmws_threshold) || hmws_threshold < 0 || hmws_threshold > 1) {
		cli::cli_abort("{.arg hmws_threshold} must be between 0 and 1.")
	}

	recipes::add_step(
		recipe,
		step_sec_aggregates_new(
			measures = measures,
			monomer_start = monomer_start,
			monomer_end = monomer_end,
			method = method,
			hmws_threshold = hmws_threshold,
			include_main_peak = include_main_peak,
			output_prefix = output_prefix,
			role = role,
			trained = trained,
			skip = skip,
			id = id
		)
	)
}

step_sec_aggregates_new <- function(
	measures,
	monomer_start,
	monomer_end,
	method,
	hmws_threshold,
	include_main_peak,
	output_prefix,
	role,
	trained,
	skip,
	id
) {
	recipes::step(
		subclass = "sec_aggregates",
		measures = measures,
		monomer_start = monomer_start,
		monomer_end = monomer_end,
		method = method,
		hmws_threshold = hmws_threshold,
		include_main_peak = include_main_peak,
		output_prefix = output_prefix,
		role = role,
		trained = trained,
		skip = skip,
		id = id
	)
}

#' @export
prep.step_sec_aggregates <- function(x, training, info = NULL, ...) {
	check_for_measure(training)

	# Find measure columns if not specified
	if (is.null(x$measures)) {
		measures <- find_measure_cols(training)
	} else {
		measures <- x$measures
	}

	step_sec_aggregates_new(
		measures = measures,
		monomer_start = x$monomer_start,
		monomer_end = x$monomer_end,
		method = x$method,
		hmws_threshold = x$hmws_threshold,
		include_main_peak = x$include_main_peak,
		output_prefix = x$output_prefix,
		role = x$role,
		trained = TRUE,
		skip = x$skip,
		id = x$id
	)
}

#' Find the tallest peak boundaries
#' @noRd
.find_main_peak <- function(location, value, threshold_frac = 0.05) {
	# Find the maximum and define peak region
	max_idx <- which.max(value)
	max_val <- value[max_idx]
	threshold <- threshold_frac * max_val

	n <- length(value)

	# Find start of peak (going backwards from max)
	start_idx <- max_idx
	for (i in seq(max_idx, 1, -1)) {
		if (value[i] < threshold) {
			start_idx <- i
			break
		}
		if (i == 1) start_idx <- 1
	}

	# Find end of peak (going forward from max)
	end_idx <- max_idx
	for (i in seq(max_idx, n)) {
		if (value[i] < threshold) {
			end_idx <- i
			break
		}
		if (i == n) end_idx <- n
	}

	list(
		start = location[start_idx],
		end = location[end_idx],
		start_idx = start_idx,
		end_idx = end_idx
	)
}

#' @export
bake.step_sec_aggregates <- function(object, new_data, ...) {
	measures <- object$measures
	monomer_start <- object$monomer_start
	monomer_end <- object$monomer_end
	method <- object$method
	hmws_threshold <- object$hmws_threshold
	include_main_peak <- object$include_main_peak
	output_prefix <- object$output_prefix

	# Initialize output columns
	n_rows <- nrow(new_data)
	hmws_pct <- numeric(n_rows)
	monomer_pct <- numeric(n_rows)
	lmws_pct <- numeric(n_rows)
	main_start <- numeric(n_rows)
	main_end <- numeric(n_rows)

	# Use first measure column for analysis
	measure_col <- measures[1]

	for (i in seq_len(n_rows)) {
		m <- new_data[[measure_col]][[i]]
		location <- m$location
		value <- m$value

		# Handle NA values
		value[is.na(value)] <- 0

		# Ensure non-negative
		value <- pmax(value, 0)

		# Determine peak boundaries
		if (
			method == "tallest" && (is.null(monomer_start) || is.null(monomer_end))
		) {
			peak_info <- .find_main_peak(location, value, threshold_frac = 0.05)
			m_start <- peak_info$start
			m_end <- peak_info$end
		} else {
			m_start <- monomer_start
			m_end <- monomer_end
		}

		# Calculate areas using trapezoidal integration
		dt <- diff(location)

		# Total area (simple sum approximation)
		total_area <- sum(value[-1] * dt + value[-length(value)] * dt) / 2

		if (total_area > 0) {
			# HMWS: area before monomer_start
			hmws_idx <- location < m_start
			if (any(hmws_idx)) {
				hmws_vals <- value[hmws_idx]
				hmws_locs <- location[hmws_idx]
				if (length(hmws_vals) > 1) {
					hmws_dt <- diff(hmws_locs)
					hmws_area <- sum(
						hmws_vals[-1] * hmws_dt + hmws_vals[-length(hmws_vals)] * hmws_dt
					) /
						2
				} else {
					hmws_area <- 0
				}
			} else {
				hmws_area <- 0
			}

			# Monomer: area between monomer_start and monomer_end
			mono_idx <- location >= m_start & location <= m_end
			if (any(mono_idx)) {
				mono_vals <- value[mono_idx]
				mono_locs <- location[mono_idx]
				if (length(mono_vals) > 1) {
					mono_dt <- diff(mono_locs)
					mono_area <- sum(
						mono_vals[-1] * mono_dt + mono_vals[-length(mono_vals)] * mono_dt
					) /
						2
				} else {
					mono_area <- 0
				}
			} else {
				mono_area <- 0
			}

			# LMWS: area after monomer_end
			lmws_idx <- location > m_end
			if (any(lmws_idx)) {
				lmws_vals <- value[lmws_idx]
				lmws_locs <- location[lmws_idx]
				if (length(lmws_vals) > 1) {
					lmws_dt <- diff(lmws_locs)
					lmws_area <- sum(
						lmws_vals[-1] * lmws_dt + lmws_vals[-length(lmws_vals)] * lmws_dt
					) /
						2
				} else {
					lmws_area <- 0
				}
			} else {
				lmws_area <- 0
			}

			# Calculate percentages
			hmws_pct[i] <- 100 * hmws_area / total_area
			monomer_pct[i] <- 100 * mono_area / total_area
			lmws_pct[i] <- 100 * lmws_area / total_area
		}

		main_start[i] <- m_start
		main_end[i] <- m_end
	}

	# Add output columns
	new_data[[paste0(output_prefix, "hmws")]] <- hmws_pct
	new_data[[paste0(output_prefix, "monomer")]] <- monomer_pct
	new_data[[paste0(output_prefix, "lmws")]] <- lmws_pct

	if (include_main_peak) {
		new_data[[paste0(output_prefix, "main_start")]] <- main_start
		new_data[[paste0(output_prefix, "main_end")]] <- main_end
	}

	tibble::as_tibble(new_data)
}

#' @export
print.step_sec_aggregates <- function(
	x,
	width = max(20, options()$width - 30),
	...
) {
	title <- "SEC aggregate quantitation"
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
tidy.step_sec_aggregates <- function(x, ...) {
	tibble::tibble(
		measures = list(x$measures),
		monomer_start = x$monomer_start %||% NA_real_,
		monomer_end = x$monomer_end %||% NA_real_,
		method = x$method,
		id = x$id
	)
}

#' @rdname required_pkgs.step_sec
#' @export
#' @keywords internal
required_pkgs.step_sec_aggregates <- function(x, ...) {
	c("measure.sec", "measure")
}
