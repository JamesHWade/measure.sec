# ==============================================================================
# SEC Data Export Functions
#
# Slice tables and summary tables for SEC analysis results
# ==============================================================================

#' Extract Slice-by-Slice SEC Data
#'
#' Extracts point-by-point (slice) data from SEC analysis results, creating
#' a long-format table suitable for export or further analysis.
#'
#' @param data A data frame containing SEC results with measure columns.
#' @param measures Character vector of measure column names to include.
#'   If `NULL`, includes all measure columns found.
#' @param sample_id Column name containing sample identifiers. If `NULL`,
#'   uses row numbers.
#' @param include_location Logical. Include the location (time/volume) column?
#'   Default is `TRUE`.
#' @param pivot Logical. Pivot measures to wide format (one column per measure)?
#'   Default is `FALSE` (long format).
#'
#' @return A tibble with slice-by-slice data:
#'   \describe{
#'     \item{sample_id}{Sample identifier}
#'     \item{slice}{Slice index (1, 2, 3, ...)}
#'     \item{location}{Elution time or volume}
#'     \item{measure}{Measure column name (if pivot = FALSE)}
#'     \item{value}{Signal value (if pivot = FALSE)}
#'     \item{<measure_names>}{Individual measure columns (if pivot = TRUE)}
#'   }
#'
#' @details
#' This function extracts the raw slice data from processed SEC chromatograms,
#' making it easy to:
#' \itemize{
#'   \item Export to CSV/Excel for external analysis
#'   \item Create custom plots
#'   \item Perform slice-level calculations
#'   \item Compare samples point-by-point
#' }
#'
#' **Typical Slice Data Columns:**
#' \itemize{
#'   \item Retention time/volume (location)
#'   \item Concentration (from RI or UV)
#'   \item Molecular weight (from calibration or MALS)
#'   \item Intrinsic viscosity (from viscometer)
#'   \item Radius of gyration (from MALS angles)
#'   \item Composition (from UV/RI ratio)
#' }
#'
#' @family sec-export
#' @export
#'
#' @examples
#' \dontrun{
#' library(recipes)
#' library(measure)
#'
#' # Process SEC data
#' prepped <- recipe(~., data = sec_triple_detect) |>
#'   step_measure_input_long(ri, location = vars(time), col_name = "ri") |>
#'   step_sec_baseline() |>
#'   prep() |>
#'   bake(new_data = NULL)
#'
#' # Extract slice table (long format)
#' slices <- measure_sec_slice_table(prepped, measures = "ri")
#'
#' # Extract slice table (wide format)
#' slices_wide <- measure_sec_slice_table(
#'   prepped,
#'   measures = c("ri", "mw"),
#'   pivot = TRUE
#' )
#'
#' # Export to CSV
#' write.csv(slices, "sec_slices.csv", row.names = FALSE)
#' }
measure_sec_slice_table <- function(
	data,
	measures = NULL,
	sample_id = NULL,
	include_location = TRUE,
	pivot = FALSE
) {
	if (!is.data.frame(data)) {
		cli::cli_abort("{.arg data} must be a data frame.")
	}

	# Find measure columns
	all_measure_cols <- find_measure_cols(data)

	if (length(all_measure_cols) == 0) {
		cli::cli_abort("No measure columns found in {.arg data}.")
	}

	if (is.null(measures)) {
		measures <- all_measure_cols
	} else {
		missing <- setdiff(measures, all_measure_cols)
		if (length(missing) > 0) {
			cli::cli_abort("Measure columns not found: {.val {missing}}.")
		}
	}

	# Determine sample IDs
	if (!is.null(sample_id) && sample_id %in% names(data)) {
		sample_ids <- data[[sample_id]]
	} else {
		sample_ids <- seq_len(nrow(data))
	}

	# Extract slice data for each sample and measure
	slice_list <- list()
	null_warnings <- character()

	for (i in seq_len(nrow(data))) {
		for (measure_col in measures) {
			m <- data[[measure_col]][[i]]

			if (is.null(m)) {
				null_warnings <- c(
					null_warnings,
					sprintf("row %d, measure '%s'", i, measure_col)
				)
				next
			}

			n_slices <- length(m$value)

			slice_df <- tibble::tibble(
				sample_id = rep(sample_ids[i], n_slices),
				slice = seq_len(n_slices),
				measure = rep(measure_col, n_slices),
				value = m$value
			)

			if (include_location) {
				slice_df$location <- m$location
				# Reorder columns
				slice_df <- slice_df[,
					c(
						"sample_id",
						"slice",
						"location",
						"measure",
						"value"
					)
				]
			}

			slice_list <- c(slice_list, list(slice_df))
		}
	}

	# Warn about NULL measures if any were found
	if (length(null_warnings) > 0) {
		cli::cli_warn(
			c(
				"Skipped {length(null_warnings)} NULL measure value{?s}:",
				"i" = "Affected: {.val {null_warnings}}"
			)
		)
	}

	result <- dplyr::bind_rows(slice_list)

	# Pivot to wide format if requested
	if (pivot && nrow(result) > 0) {
		if (include_location) {
			result <- tidyr::pivot_wider(
				result,
				id_cols = c("sample_id", "slice", "location"),
				names_from = "measure",
				values_from = "value"
			)
		} else {
			result <- tidyr::pivot_wider(
				result,
				id_cols = c("sample_id", "slice"),
				names_from = "measure",
				values_from = "value"
			)
		}
	}

	result
}

#' Generate SEC Summary Table
#'
#' Creates a summary table of SEC analysis results with key metrics for
#' each sample.
#'
#' @param data A data frame containing SEC results.
#' @param mw_col Column name containing molecular weight averages (list column
#'   with Mn, Mw, Mz, dispersity).
#' @param include_mw Logical. Include molecular weight averages? Default is `TRUE`.
#' @param include_fractions Logical. Include MW fractions if available?
#'   Default is `TRUE`.
#' @param include_purity Logical. Include purity metrics (HMWS, monomer, LMWS)
#'   if available? Default is `TRUE`.
#' @param sample_id Column name for sample identifiers.
#' @param additional_cols Character vector of additional columns to include
#'   in the summary.
#' @param digits Number of decimal places for numeric columns. Default is 2.
#'
#' @return A tibble with one row per sample containing:
#'   \describe{
#'     \item{sample_id}{Sample identifier}
#'     \item{Mn}{Number-average molecular weight}
#'     \item{Mw}{Weight-average molecular weight}
#'     \item{Mz}{Z-average molecular weight}
#'     \item{dispersity}{Polydispersity index (Mw/Mn)}
#'     \item{purity_hmws}{Percent high MW species (if available)}
#'     \item{purity_monomer}{Percent monomer (if available)}
#'     \item{purity_lmws}{Percent low MW species (if available)}
#'   }
#'
#' @details
#' This function creates a publication-ready summary table of SEC results.
#' It automatically detects and includes available metrics.
#'
#' **Typical Summary Metrics:**
#' \itemize{
#'   \item Molecular weight averages: Mn, Mw, Mz
#'   \item Dispersity (PDI): Mw/Mn
#'   \item Purity metrics: %HMWS, %Monomer, %LMWS
#'   \item MW fractions: % above/below cutoffs
#'   \item Recovery: % mass balance
#' }
#'
#' @family sec-export
#' @export
#'
#' @examples
#' \dontrun{
#' # Generate summary after SEC processing
#' summary_tbl <- measure_sec_summary_table(
#'   processed_data,
#'   sample_id = "sample_name"
#' )
#'
#' # Print formatted table
#' print(summary_tbl)
#'
#' # Export to Excel
#' writexl::write_xlsx(summary_tbl, "sec_summary.xlsx")
#' }
measure_sec_summary_table <- function(
	data,
	mw_col = NULL,
	include_mw = TRUE,
	include_fractions = TRUE,
	include_purity = TRUE,
	sample_id = NULL,
	additional_cols = NULL,
	digits = 2
) {
	if (!is.data.frame(data)) {
		cli::cli_abort("{.arg data} must be a data frame.")
	}

	n_samples <- nrow(data)

	# Initialize result with sample IDs
	if (!is.null(sample_id) && sample_id %in% names(data)) {
		result <- tibble::tibble(sample_id = data[[sample_id]])
	} else {
		result <- tibble::tibble(sample_id = seq_len(n_samples))
	}

	# Add molecular weight columns if available
	if (include_mw) {
		mw_cols <- c(
			"Mn",
			"Mw",
			"Mz",
			"dispersity",
			"mw_mn",
			"mw_mw",
			"mw_mz",
			"mw_dispersity"
		)
		for (col in mw_cols) {
			if (col %in% names(data)) {
				result[[col]] <- round(data[[col]], digits)
			}
		}
	}

	# Add purity columns if available
	if (include_purity) {
		purity_cols <- c("purity_hmws", "purity_monomer", "purity_lmws")
		for (col in purity_cols) {
			if (col %in% names(data)) {
				result[[col]] <- round(data[[col]], digits)
			}
		}
	}

	# Add fraction columns if available
	if (include_fractions) {
		# Look for fraction columns (pattern: frac_* or *_fraction)
		frac_cols <- names(data)[
			grepl(
				"frac|fraction",
				names(data),
				ignore.case = TRUE
			)
		]
		for (col in frac_cols) {
			if (is.numeric(data[[col]])) {
				result[[col]] <- round(data[[col]], digits)
			}
		}
	}

	# Add additional columns if specified
	if (!is.null(additional_cols)) {
		for (col in additional_cols) {
			if (col %in% names(data)) {
				if (is.numeric(data[[col]])) {
					result[[col]] <- round(data[[col]], digits)
				} else {
					result[[col]] <- data[[col]]
				}
			}
		}
	}

	# Add class for custom printing
	class(result) <- c("sec_summary_table", class(result))

	result
}

#' Validate sec_summary_table structure
#' @noRd
validate_sec_summary_table <- function(x) {
	if (!inherits(x, "tbl_df")) {
		cli::cli_warn("Invalid sec_summary_table: not a tibble.")
		return(FALSE)
	}
	if (!"sample_id" %in% names(x)) {
		cli::cli_warn("Invalid sec_summary_table: missing sample_id column.")
		return(FALSE)
	}
	TRUE
}

#' @export
print.sec_summary_table <- function(x, ...) {
	validate_sec_summary_table(x)
	cat("SEC Analysis Summary\n")
	cat(strrep("=", 60), "\n")
	cat("Samples:", nrow(x), "\n\n")

	# Print as tibble (removes custom class temporarily)
	print(tibble::as_tibble(x), ...)

	invisible(x)
}

# ==============================================================================
# Multi-Sample Comparison
# ==============================================================================

#' Compare Multiple SEC Samples
#'
#' Compares SEC analysis results across multiple samples, providing a summary
#' table of key metrics and optional overlay plots.
#'
#' @param ... Data frames containing SEC results. Each should be a baked recipe
#'   output with measure columns (MW, concentration, etc.).
#' @param samples Character vector of sample names. If `NULL`, uses names from
#'   `...` or generates sequential names ("Sample 1", "Sample 2", etc.).
#' @param metrics Character vector specifying which metrics to compare. Options:
#'   \itemize{
#'     \item `"mw_averages"`: Mn, Mw, Mz, dispersity
#'     \item `"mwd"`: Molecular weight distribution overlap
#'     \item `"branching"`: Branching metrics (if available)
#'   }
#'   Default includes all available metrics.
#' @param plot Logical. Generate MWD comparison plot? Default is `TRUE`.
#'   Note: Plot is only generated when `"mwd"` is included in `metrics`.
#' @param reference Integer or character. Which sample to use as reference for
#'   percent differences. Default is `1` (first sample).
#' @param digits Integer. Number of decimal places for numeric values.
#'   Default is `2`. Must be non-negative.
#'
#' @return A list of class `sec_comparison` containing:
#'   \describe{
#'     \item{summary}{Tibble with comparison metrics for all samples}
#'     \item{differences}{Tibble with absolute and percent differences vs reference}
#'     \item{plot}{ggplot2 object (if `plot = TRUE` and `"mwd"` in `metrics`),
#'       otherwise `NULL`}
#'     \item{samples}{Character vector of sample names}
#'     \item{reference}{Name of reference sample}
#'   }
#'
#' @details
#' This function is useful for:
#' \itemize{
#'   \item Batch-to-batch comparison
#'   \item Stability studies
#'   \item Process optimization
#'   \item Quality control
#' }
#'
#' **Input Data Handling:**
#'
#' When input data frames contain multiple rows (e.g., multiple injections),
#' numeric metrics are averaged across rows. For single-row data frames
#' (typical for processed SEC results), values are used directly.
#'
#' The function recognizes molecular weight columns with either naming
#' convention: prefixed (`mw_mn`, `mw_mw`, `mw_mz`, `mw_dispersity`) or
#' standard (`Mn`, `Mw`, `Mz`, `dispersity`).
#'
#' **Comparison Metrics:**
#'
#' *MW Averages:* Compares Mn, Mw, Mz, and dispersity with absolute and
#' percent differences from the reference sample.
#'
#' *MWD Overlay:* Creates overlaid MWD plots showing distribution differences.
#' Useful for detecting bimodality, tailing, or distribution shifts.
#'
#' *Branching:* Compares branching index and frequency if available from
#' triple-detection data.
#'
#' @family sec-export
#' @export
#'
#' @examples
#' \dontrun{
#' # Compare three batches
#' comparison <- measure_sec_compare(
#'   batch1_data,
#'   batch2_data,
#'   batch3_data,
#'   samples = c("Batch 1", "Batch 2", "Batch 3")
#' )
#'
#' # View summary table
#' comparison$summary
#'
#' # View differences from reference
#' comparison$differences
#'
#' # Display overlay plot
#' comparison$plot
#'
#' # Compare only MW averages without plot
#' comparison <- measure_sec_compare(
#'   sample_a, sample_b,
#'   metrics = "mw_averages",
#'   plot = FALSE
#' )
#' }
measure_sec_compare <- function(
	...,
	samples = NULL,
	metrics = c("mw_averages", "mwd", "branching"),
	plot = TRUE,
	reference = 1,
	digits = 2
) {
	# Capture input data frames
	data_list <- list(...)

	if (length(data_list) < 2) {
		cli::cli_abort("At least 2 samples are required for comparison.")
	}

	# Validate all inputs are data frames
	for (i in seq_along(data_list)) {
		if (!is.data.frame(data_list[[i]])) {
			cli::cli_abort(
				"All inputs must be data frames. Input {i} is not a data frame."
			)
		}
	}

	# Validate parameters
	if (!is.logical(plot) || length(plot) != 1 || is.na(plot)) {
		cli::cli_abort("{.arg plot} must be TRUE or FALSE.")
	}

	if (!is.numeric(digits) || length(digits) != 1 || digits < 0) {
		cli::cli_abort("{.arg digits} must be a non-negative number.")
	}
	digits <- as.integer(digits)

	# Validate metrics

	valid_metrics <- c("mw_averages", "mwd", "branching")
	invalid_metrics <- setdiff(metrics, valid_metrics)
	if (length(invalid_metrics) > 0) {
		cli::cli_warn(
			c(
				"Unrecognized metrics ignored: {.val {invalid_metrics}}",
				"i" = "Valid options: {.val {valid_metrics}}"
			)
		)
	}
	metrics <- intersect(metrics, valid_metrics)
	if (length(metrics) == 0) {
		cli::cli_abort(
			c(
				"No valid metrics specified.",
				"i" = "Valid options: {.val {valid_metrics}}"
			)
		)
	}

	# Generate sample names
	if (is.null(samples)) {
		# Try to get names from ... arguments
		call_args <- as.list(match.call())[-1]
		call_args <- call_args[
			!names(call_args) %in%
				c("samples", "metrics", "plot", "reference", "digits")
		]
		arg_names <- vapply(call_args, deparse, character(1))

		# Use argument names if they look like variable names, otherwise generate
		samples <- ifelse(
			grepl("^[a-zA-Z][a-zA-Z0-9_.]*$", arg_names),
			arg_names,
			paste("Sample", seq_along(data_list))
		)
	}

	if (length(samples) != length(data_list)) {
		cli::cli_abort(
			"{.arg samples} length ({length(samples)}) must match number of data inputs ({length(data_list)})."
		)
	}

	# Resolve reference
	if (is.character(reference)) {
		ref_idx <- match(reference, samples)
		if (is.na(ref_idx)) {
			cli::cli_abort(
				"Reference sample {.val {reference}} not found in sample names."
			)
		}
	} else {
		ref_idx <- as.integer(reference)
		if (ref_idx < 1 || ref_idx > length(data_list)) {
			cli::cli_abort(
				"{.arg reference} must be between 1 and {length(data_list)}."
			)
		}
	}

	# Build comparison summary
	summary_list <- list()

	for (i in seq_along(data_list)) {
		df <- data_list[[i]]
		row <- tibble::tibble(sample = samples[i])

		# MW averages
		if ("mw_averages" %in% metrics) {
			mw_cols <- c(
				"mw_mn",
				"mw_mw",
				"mw_mz",
				"mw_dispersity",
				"Mn",
				"Mw",
				"Mz",
				"dispersity"
			)
			for (col in mw_cols) {
				if (col %in% names(df) && is.numeric(df[[col]])) {
					# Take first value if multiple rows, or mean
					val <- if (nrow(df) == 1) df[[col]] else mean(df[[col]], na.rm = TRUE)
					row[[col]] <- round(val, digits)
				}
			}
		}

		# Branching metrics
		if ("branching" %in% metrics) {
			branch_cols <- c(
				"branching_index",
				"branching_frequency",
				"g_ratio",
				"g_prime"
			)
			for (col in branch_cols) {
				if (col %in% names(df) && is.numeric(df[[col]])) {
					val <- if (nrow(df) == 1) df[[col]] else mean(df[[col]], na.rm = TRUE)
					row[[col]] <- round(val, digits)
				}
			}
		}

		summary_list[[i]] <- row
	}

	summary_tbl <- dplyr::bind_rows(summary_list)

	# Calculate differences from reference
	diff_tbl <- calculate_differences(summary_tbl, ref_idx, digits)

	# Create overlay plot if requested
	plot_obj <- NULL
	if (plot && "mwd" %in% metrics) {
		plot_obj <- create_comparison_plot(data_list, samples)
	}

	# Build result object
	result <- list(
		summary = summary_tbl,
		differences = diff_tbl,
		plot = plot_obj,
		samples = samples,
		reference = samples[ref_idx]
	)

	class(result) <- c("sec_comparison", "list")
	result
}

#' Calculate differences from reference sample
#' @noRd
calculate_differences <- function(summary_tbl, ref_idx, digits = 2) {
	ref_row <- summary_tbl[ref_idx, ]
	numeric_cols <- names(summary_tbl)[
		vapply(
			summary_tbl,
			is.numeric,
			logical(1)
		)
	]

	if (length(numeric_cols) == 0) {
		return(tibble::tibble(sample = summary_tbl$sample))
	}

	diff_list <- list()

	for (i in seq_len(nrow(summary_tbl))) {
		row <- tibble::tibble(sample = summary_tbl$sample[i])

		for (col in numeric_cols) {
			val <- summary_tbl[[col]][i]
			ref_val <- ref_row[[col]]

			if (!is.na(val) && !is.na(ref_val) && ref_val != 0) {
				# Absolute difference
				row[[paste0(col, "_diff")]] <- round(val - ref_val, digits)
				# Percent difference
				row[[paste0(col, "_pct")]] <- round((val - ref_val) / ref_val * 100, 1)
			}
		}

		diff_list[[i]] <- row
	}

	dplyr::bind_rows(diff_list)
}

#' Create MWD overlay comparison plot
#' @noRd
create_comparison_plot <- function(data_list, samples) {
	if (!requireNamespace("ggplot2", quietly = TRUE)) {
		cli::cli_warn("Package {.pkg ggplot2} is required for plotting.")
		return(NULL)
	}

	# Combine data for plotting
	plot_data_list <- list()
	skipped_no_mw <- character()
	skipped_invalid_rows <- list()
	filtered_non_positive <- list()

	for (i in seq_along(data_list)) {
		df <- data_list[[i]]

		# Check for MW column
		mw_col <- if ("mw" %in% names(df)) "mw" else NULL

		if (is.null(mw_col)) {
			skipped_no_mw <- c(skipped_no_mw, samples[i])
			next
		}

		# Extract MW distribution data
		if (
			inherits(df[[mw_col]], "measure_list") ||
				(is.list(df[[mw_col]]) && length(df[[mw_col]]) > 0)
		) {
			# Handle measure_list format
			invalid_rows <- 0
			total_filtered <- 0

			for (j in seq_len(nrow(df))) {
				m <- df[[mw_col]][[j]]
				if (is.null(m) || is.null(m$location) || is.null(m$value)) {
					invalid_rows <- invalid_rows + 1
					next
				}

				slice_df <- tibble::tibble(
					sample = samples[i],
					location = m$location,
					mw = m$value
				)

				# Filter positive MW values and track count
				original_rows <- nrow(slice_df)
				slice_df <- slice_df[slice_df$mw > 0, ]
				total_filtered <- total_filtered + (original_rows - nrow(slice_df))

				if (nrow(slice_df) > 0) {
					plot_data_list <- c(plot_data_list, list(slice_df))
				}
			}

			if (invalid_rows > 0) {
				skipped_invalid_rows[[samples[i]]] <- invalid_rows
			}
			if (total_filtered > 0) {
				filtered_non_positive[[samples[i]]] <- total_filtered
			}
		}
	}

	# Report warnings for skipped/filtered data
	if (length(skipped_no_mw) > 0) {
		cli::cli_warn(
			c(
				"Sample{?s} excluded from plot (no {.field mw} column): {.val {skipped_no_mw}}",
				"i" = "MWD plot requires an {.field mw} measure column."
			)
		)
	}

	if (length(skipped_invalid_rows) > 0) {
		for (sample_name in names(skipped_invalid_rows)) {
			cli::cli_warn(
				"Sample {.val {sample_name}}: Skipped {skipped_invalid_rows[[sample_name]]} row{?s} with missing/invalid MW data."
			)
		}
	}

	if (length(filtered_non_positive) > 0) {
		for (sample_name in names(filtered_non_positive)) {
			n_filtered <- filtered_non_positive[[sample_name]]
			if (n_filtered > 10) {
				cli::cli_warn(
					c(
						"Sample {.val {sample_name}}: {n_filtered} points with non-positive MW values excluded.",
						"!" = "This may indicate calibration issues."
					)
				)
			}
		}
	}

	if (length(plot_data_list) == 0) {
		cli::cli_warn("No valid MW data found for plotting.")
		return(NULL)
	}

	plot_data <- dplyr::bind_rows(plot_data_list)

	# Create overlay plot
	p <- ggplot2::ggplot(
		plot_data,
		ggplot2::aes(
			x = log10(.data$mw),
			y = .data$location,
			color = .data$sample
		)
	) +
		ggplot2::geom_line(linewidth = 0.8) +
		ggplot2::labs(
			x = expression(log[10](M)),
			y = "Signal",
			color = "Sample",
			title = "MWD Comparison"
		) +
		ggplot2::theme_bw() +
		ggplot2::theme(
			legend.position = "bottom",
			plot.title = ggplot2::element_text(hjust = 0.5)
		)

	p
}

#' @export
print.sec_comparison <- function(x, ...) {
	cat("SEC Multi-Sample Comparison\n")
	cat(strrep("=", 60), "\n")
	cat("Samples:", length(x$samples), "\n")
	cat("Reference:", x$reference, "\n\n")

	cat("Summary:\n")
	print(x$summary, ...)

	if (ncol(x$differences) > 1) {
		cat("\nDifferences from reference:\n")
		print(x$differences, ...)
	}

	if (!is.null(x$plot)) {
		cat("\nPlot available in $plot\n")
	}

	invisible(x)
}
