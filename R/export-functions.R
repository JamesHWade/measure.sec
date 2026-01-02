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
        slice_df <- slice_df[, c(
          "sample_id",
          "slice",
          "location",
          "measure",
          "value"
        )]
      }

      slice_list <- c(slice_list, list(slice_df))
    }
  }

  # Warn about NULL measures if any were found
  if (length(null_warnings) > 0) {
    cli::cli_warn(c(
      "Skipped {length(null_warnings)} NULL measure value{?s}:",
      "i" = "Affected: {.val {null_warnings}}"
    ))
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
    frac_cols <- names(data)[grepl(
      "frac|fraction",
      names(data),
      ignore.case = TRUE
    )]
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
