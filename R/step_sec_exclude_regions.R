# ==============================================================================
# step_sec_exclude_regions
#
# Mark regions for exclusion from baseline fitting and/or MW integration
# ==============================================================================

#' SEC/GPC Region Exclusion
#'
#' `step_sec_exclude_regions()` creates a *specification* of a recipe step
#' that marks regions for exclusion from analysis. Excluded regions can be
#' solvent peaks, artifacts, flow markers, or other features that should not
#' be included in baseline fitting or molecular weight integration.
#'
#' @param recipe A recipe object. The step will be added to the sequence of
#'   operations for this recipe.
#' @param measures An optional character vector of measure column names to
#'   process. If `NULL` (the default), all measure columns will be processed.
#' @param regions A data frame or tibble specifying regions to exclude. Must
#'   contain columns `start` and `end` (x-axis values defining each region).
#'   Optional columns: `reason` (character describing why excluded),
#'   `sample_id` (for sample-specific exclusions).
#' @param purpose Character. What the exclusions apply to:
#'   - `"baseline"`: Exclude from baseline fitting only
#'   - `"integration"`: Exclude from MW integration only
#'   - `"both"` (default): Exclude from both baseline fitting and integration
#' @param role Not used by this step since no new variables are created.
#' @param trained A logical to indicate if the quantities for preprocessing
#'   have been estimated.
#' @param skip A logical. Should the step be skipped when the recipe is baked?
#' @param id A character string that is unique to this step to identify it.
#'
#' @return An updated version of `recipe` with the new step added to the
#'   sequence of any existing operations. An `.excluded_regions` column
#'   will be added containing exclusion information for each sample.
#'
#' @details
#' Excluded regions are stored in a list column `.excluded_regions` where
#' each element is a tibble with columns:
#' - `start`: Start of excluded region (x-axis value)
#' - `end`: End of excluded region (x-axis value)
#' - `purpose`: What this exclusion applies to
#' - `reason`: Optional description of why this region is excluded
#'
#' **Sample-specific exclusions:**
#' If the `regions` data frame includes a `sample_id` column, exclusions
#' can be applied only to specific samples. Rows without `sample_id` (or with
#' `sample_id = NA`) are applied globally to all samples.
#'
#' **Common uses:**
#' - Solvent peaks that elute at the end of the chromatogram
#' - System peaks or artifacts at specific retention times
#' - Flow marker peaks used for retention time correction
#' - Air bubbles or other transient artifacts
#'
#' # Tidying
#'
#' When you [`tidy()`][recipes::tidy.recipe()] this step, a tibble with columns
#' `terms`, `n_regions`, `purpose`, and `id` is returned.
#'
#' @seealso [step_sec_baseline()] which can use exclusion information,
#'   [step_sec_mw_averages()] which respects integration exclusions.
#'
#' @family sec-chromatography
#' @export
#'
#' @examples
#' \dontrun{
#' library(recipes)
#' library(measure)
#'
#' # Exclude a single region (solvent peak at end)
#' exclusions <- tibble::tibble(
#'   start = 18.5,
#'   end = 20.0,
#'   reason = "Solvent peak"
#' )
#'
#' rec <- recipe(~., data = sec_data) |>
#'   step_measure_input_long(ri_signal, location = vars(elution_volume)) |>
#'   step_sec_exclude_regions(regions = exclusions) |>
#'   prep()
#'
#' # Exclude multiple regions
#' exclusions <- tibble::tibble(
#'   start = c(8.0, 18.5),
#'   end = c(9.0, 20.0),
#'   reason = c("Void volume", "Solvent peak")
#' )
#'
#' rec <- recipe(~., data = sec_data) |>
#'   step_measure_input_long(ri_signal, location = vars(elution_volume)) |>
#'   step_sec_exclude_regions(regions = exclusions, purpose = "integration") |>
#'   prep()
#'
#' # Sample-specific exclusions
#' exclusions <- tibble::tibble(
#'   start = c(18.5, 15.0),
#'   end = c(20.0, 16.0),
#'   sample_id = c(NA, "sample_2"),  # NA = global, specific sample_id = per-sample
#'   reason = c("Solvent peak", "Artifact in sample_2")
#' )
#'
#' rec <- recipe(~., data = sec_data) |>
#'   step_measure_input_long(ri_signal, location = vars(elution_volume)) |>
#'   step_sec_exclude_regions(regions = exclusions) |>
#'   prep()
#' }
step_sec_exclude_regions <- function(
	recipe,
	measures = NULL,
	regions = NULL,
	purpose = c("both", "baseline", "integration"),
	role = NA,
	trained = FALSE,
	skip = FALSE,
	id = recipes::rand_id("sec_exclude_regions")
) {
	purpose <- rlang::arg_match(purpose)

	# Validate regions

	if (!is.null(regions)) {
		if (!is.data.frame(regions)) {
			cli::cli_abort("{.arg regions} must be a data frame or tibble.")
		}
		if (!all(c("start", "end") %in% names(regions))) {
			cli::cli_abort(
				"{.arg regions} must contain columns {.field start} and {.field end}."
			)
		}
		if (!is.numeric(regions$start) || !is.numeric(regions$end)) {
			cli::cli_abort(
				"Columns {.field start} and {.field end} must be numeric."
			)
		}
		# Check start < end for each region
		invalid_regions <- regions$start >= regions$end
		if (any(invalid_regions)) {
			n_invalid <- sum(invalid_regions)
			cli::cli_abort(
				c(
					"{n_invalid} region{?s} ha{?s/ve} {.field start} >= {.field end}.",
					"i" = "Each region's {.field start} must be less than {.field end}."
				)
			)
		}
	}

	recipes::add_step(
		recipe,
		step_sec_exclude_regions_new(
			measures = measures,
			regions = regions,
			purpose = purpose,
			role = role,
			trained = trained,
			skip = skip,
			id = id
		)
	)
}

step_sec_exclude_regions_new <- function(
	measures,
	regions,
	purpose,
	role,
	trained,
	skip,
	id
) {
	recipes::step(
		subclass = "sec_exclude_regions",
		measures = measures,
		regions = regions,
		purpose = purpose,
		role = role,
		trained = trained,
		skip = skip,
		id = id
	)
}

#' @export
prep.step_sec_exclude_regions <- function(x, training, info = NULL, ...) {
	check_for_measure(training)

	# Resolve which columns to process
	if (is.null(x$measures)) {
		measure_cols <- find_measure_cols(training)
	} else {
		measure_cols <- x$measures
	}

	step_sec_exclude_regions_new(
		measures = measure_cols,
		regions = x$regions,
		purpose = x$purpose,
		role = x$role,
		trained = TRUE,
		skip = x$skip,
		id = x$id
	)
}

#' @export
bake.step_sec_exclude_regions <- function(object, new_data, ...) {
	regions <- object$regions
	purpose <- object$purpose

	n_rows <- nrow(new_data)

	# Check if data has sample_id column for per-sample exclusions
	has_sample_id <- "sample_id" %in% names(new_data)

	# Build exclusion list for each row
	exclusion_list <- lapply(seq_len(n_rows), function(i) {
		if (is.null(regions) || nrow(regions) == 0) {
			# Return empty tibble with correct structure
			return(
				tibble::tibble(
					start = numeric(0),
					end = numeric(0),
					purpose = character(0),
					reason = character(0)
				)
			)
		}

		# Determine which exclusions apply to this row
		if (has_sample_id && "sample_id" %in% names(regions)) {
			sample_id_val <- new_data$sample_id[i]
			# Include global exclusions (NA sample_id) and sample-specific ones
			applicable <- is.na(regions$sample_id) |
				regions$sample_id == sample_id_val
			row_regions <- regions[applicable, , drop = FALSE]
		} else {
			# All regions apply
			row_regions <- regions
		}

		if (nrow(row_regions) == 0) {
			return(
				tibble::tibble(
					start = numeric(0),
					end = numeric(0),
					purpose = character(0),
					reason = character(0)
				)
			)
		}

		# Build output tibble
		tibble::tibble(
			start = row_regions$start,
			end = row_regions$end,
			purpose = rep(purpose, nrow(row_regions)),
			reason = if ("reason" %in% names(row_regions)) {
				row_regions$reason
			} else {
				rep(NA_character_, nrow(row_regions))
			}
		)
	})

	new_data[[".excluded_regions"]] <- exclusion_list

	tibble::as_tibble(new_data)
}

#' @export
print.step_sec_exclude_regions <- function(
	x,
	width = max(20, options()$width - 30),
	...
) {
	n_regions <- if (is.null(x$regions)) 0L else nrow(x$regions)

	if (n_regions == 0) {
		cat("SEC exclude regions [no regions defined]")
	} else {
		cat(
			glue::glue(
				"SEC exclude {n_regions} region{if (n_regions > 1) 's' else ''} ",
				"(purpose: {x$purpose})"
			)
		)
	}
	cat("\n")

	invisible(x)
}

#' @rdname tidy.step_sec
#' @export
#' @keywords internal
tidy.step_sec_exclude_regions <- function(x, ...) {
	if (recipes::is_trained(x)) {
		terms_val <- x$measures
	} else {
		terms_val <- "<all measure columns>"
	}

	n_regions <- if (is.null(x$regions)) 0L else nrow(x$regions)

	tibble::tibble(
		terms = terms_val,
		n_regions = n_regions,
		purpose = x$purpose,
		id = x$id
	)
}

#' @rdname required_pkgs.step_sec
#' @export
#' @keywords internal
required_pkgs.step_sec_exclude_regions <- function(x, ...) {
	c("measure.sec", "measure")
}
