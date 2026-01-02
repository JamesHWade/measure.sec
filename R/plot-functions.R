# ==============================================================================
# SEC Visualization Functions
#
# Publication-ready plots for SEC/GPC analysis results
# ==============================================================================

#' @importFrom tidyr pivot_wider pivot_longer
NULL


# ==============================================================================
# sec_results Class
# ==============================================================================

#' Create SEC Results Object
#'
#' Constructor for the `sec_results` class, which wraps processed SEC/GPC data
#' and enables ggplot2's `autoplot()` functionality for automatic visualization.
#'
#' @param data A data frame containing processed SEC results with measure columns.
#'   Typically the output from `bake()` on a prepped SEC recipe.
#' @param sample_id Optional. Column name containing sample identifiers.
#'   If `NULL`, auto-detection is attempted.
#'
#' @return An object of class `sec_results` (inherits from `tbl_df`).
#'
#' @details
#' The `sec_results` class provides a unified interface for SEC/GPC data that
#' enables:
#' \itemize{
#'   \item Automatic plot selection via `autoplot()`
#'   \item Integration with ggplot2 theming
#'   \item Summary statistics access
#' }
#'
#' **Expected Data Structure:**
#' The input data should contain measure columns (list columns with `location`
#' and `value` components). Common measure columns include:
#' \itemize{
#'   \item `ri`, `uv`, `mals` - Detector signals
#'   \item `mw` - Molecular weight from calibration
#'   \item `concentration` - Concentration profile
#'   \item `intrinsic_visc` - Intrinsic viscosity
#'   \item `rg` - Radius of gyration
#' }
#'
#' @family sec-visualization
#' @export
#'
#' @examples
#' \dontrun{
#' library(recipes)
#' library(measure)
#' library(ggplot2)
#'
#' # Process SEC data
#' processed <- recipe(~., data = sec_triple_detect) |>
#'   step_measure_input_long(ri_signal, location = vars(elution_time), col_name = "ri") |>
#'   step_sec_baseline(measures = "ri") |>
#'   step_sec_conventional_cal(standards = ps_standards) |>
#'   prep() |>
#'   bake(new_data = NULL)
#'
#' # Wrap as sec_results
#' results <- sec_results(processed, sample_id = "sample_id")
#'
#' # Use autoplot for automatic visualization
#' autoplot(results)
#' autoplot(results, type = "mwd")
#' autoplot(results, type = "chromatogram", normalize = TRUE)
#' }
sec_results <- function(data, sample_id = NULL) {
  if (!is.data.frame(data)) {
    cli::cli_abort("{.arg data} must be a data frame.")
  }

  # Validate that we have measure columns
  measure_cols <- find_measure_cols(data)
  if (length(measure_cols) == 0) {
    cli::cli_abort(c(
      "No measure columns found in {.arg data}.",
      "i" = "SEC results should contain processed measure columns.",
      "i" = "Use {.fn bake} on a prepped SEC recipe to create processed data."
    ))
  }

  # Auto-detect sample_id if not provided
  if (is.null(sample_id)) {
    potential_ids <- c("sample_id", "sample", "id", "sample_name", "name")
    for (col in potential_ids) {
      if (col %in% names(data)) {
        sample_id <- col
        break
      }
    }
  }

  # Construct the object
  new_sec_results(data, sample_id = sample_id, measure_cols = measure_cols)
}

#' Low-level constructor for sec_results
#' @noRd
new_sec_results <- function(
  x,
  sample_id = NULL,
  measure_cols = character()
) {
  stopifnot(is.data.frame(x))

  structure(
    tibble::as_tibble(x),
    class = c("sec_results", class(tibble::tibble())),
    sample_id = sample_id,
    measure_cols = measure_cols
  )
}

#' @export
print.sec_results <- function(x, ...) {
  sample_id_col <- attr(x, "sample_id")
  measure_cols <- attr(x, "measure_cols")

  cli::cli_h1("SEC Analysis Results")
  cli::cli_bullets(c(
    "*" = "Samples: {nrow(x)}",
    "*" = "Measure columns: {.val {measure_cols}}",
    "*" = "Sample ID column: {.val {sample_id_col %||% 'none'}}"
  ))
  cli::cli_text("")

  # Print as tibble
  NextMethod()
}


# ==============================================================================
# autoplot Method
# ==============================================================================

#' Automatic Plot for SEC Results
#'
#' Creates a ggplot2 visualization appropriate for SEC/GPC analysis results.
#' Automatically selects the best plot type based on available data, or
#' allows explicit selection via the `type` argument.
#'
#' @param object An `sec_results` object created by [sec_results()].
#' @param type Type of plot to create. One of:
#'   - `"auto"`: Automatically detect best plot type (default)
#'   - `"chromatogram"`: Basic chromatogram (signal vs time)
#'   - `"mwd"`: Molecular weight distribution
#'   - `"conformation"`: Rg-MW or eta-MW scaling plot
#'   - `"composition"`: UV/RI ratio or composition plot
#' @param overlay_mw Logical. For chromatogram plots, overlay molecular weight
#'   on secondary y-axis? Default is `TRUE` when MW data is available.
#' @param detectors Character vector of detector columns to plot for
#'   multi-detector overlays. Default is `c("ri", "uv", "mals")`.
#' @param log_scale Character. Apply log scale to axes. Options:
#'   - `"x"`: Log scale on x-axis (default for MWD plots)
#'   - `"y"`: Log scale on y-axis
#'   - `"both"`: Log scale on both axes
#'   - `"none"`: No log scaling
#' @param ... Additional arguments passed to the underlying plot function.
#'
#' @return A ggplot2 object.
#'
#' @details
#' When `type = "auto"` (default), the plot type is selected based on
#' available data:
#' \enumerate{
#'   \item If `mw` column present: MWD plot
#'   \item If multiple detectors: Multi-detector overlay
#'   \item Otherwise: Basic chromatogram
#' }
#'
#' The resulting ggplot2 object can be further customized with standard

#' ggplot2 functions like `+ theme_bw()` or `+ labs()`.
#'
#' @family sec-visualization
#' @export
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#'
#' # Create sec_results object
#' results <- sec_results(processed_sec_data)
#'
#' # Auto-detect best plot type
#' autoplot(results)
#'
#' # Specific plot types
#' autoplot(results, type = "chromatogram")
#' autoplot(results, type = "mwd", show_averages = TRUE)
#' autoplot(results, type = "conformation")
#'
#' # Customize with ggplot2
#' autoplot(results, type = "mwd") +
#'   theme_classic() +
#'   labs(title = "Molecular Weight Distribution")
#' }
autoplot.sec_results <- function(
  object,
  type = c("auto", "chromatogram", "mwd", "conformation", "composition"),
  overlay_mw = TRUE,
  detectors = c("ri", "uv", "mals"),
  log_scale = c("x", "none", "y", "both"),
  ...
) {
  check_ggplot2_available()

  type <- match.arg(type)
  log_scale <- match.arg(log_scale)

  # Get measure columns from attributes
  measure_cols <- attr(object, "measure_cols")
  sample_id <- attr(object, "sample_id")

  if (is.null(measure_cols) || length(measure_cols) == 0) {
    measure_cols <- find_measure_cols(object)
  }

  # Auto-detect plot type
  if (type == "auto") {
    type <- detect_plot_type(measure_cols)
  }

  # Dispatch to appropriate plot function
  p <- switch(
    type,
    "chromatogram" = autoplot_chromatogram(
      object,
      sample_id = sample_id,
      overlay_mw = overlay_mw,
      detectors = detectors,
      ...
    ),
    "mwd" = autoplot_mwd(
      object,
      sample_id = sample_id,
      log_scale = log_scale,
      ...
    ),
    "conformation" = autoplot_conformation(
      object,
      sample_id = sample_id,
      ...
    ),
    "composition" = autoplot_composition(
      object,
      sample_id = sample_id,
      ...
    ),
    cli::cli_abort("Unknown plot type: {.val {type}}")
  )

  p
}

#' Detect best plot type based on available measures
#' @noRd
detect_plot_type <- function(measure_cols) {
  # Check for MW-related columns -> MWD plot
  if ("mw" %in% measure_cols) {
    return("mwd")
  }

  # Check for conformation data -> conformation plot
  if ("rg" %in% measure_cols || "intrinsic_visc" %in% measure_cols) {
    return("conformation")
  }

  # Check for composition data
  if ("composition" %in% measure_cols || "uv_ri_ratio" %in% measure_cols) {
    return("composition")
  }

  # Default to chromatogram
  "chromatogram"
}

#' Internal: Chromatogram plot for autoplot
#' @noRd
autoplot_chromatogram <- function(
  object,
  sample_id = NULL,
  overlay_mw = TRUE,
  detectors = c("ri", "uv", "mals"),
  ...
) {
  measure_cols <- attr(object, "measure_cols") %||% find_measure_cols(object)

  # Find available detectors
  available_detectors <- intersect(detectors, measure_cols)

  if (length(available_detectors) == 0) {
    # Use first available measure
    available_detectors <- measure_cols[1]
  }

  if (length(available_detectors) > 1) {
    # Multi-detector overlay
    plot_sec_multidetector(
      object,
      detectors = available_detectors,
      sample_id = sample_id,
      ...
    )
  } else {
    # Single detector chromatogram
    plot_sec_chromatogram(
      object,
      measures = available_detectors,
      sample_id = sample_id,
      ...
    )
  }
}

#' Internal: MWD plot for autoplot
#' @noRd
autoplot_mwd <- function(object, sample_id = NULL, log_scale = "x", ...) {
  # Use existing MWD plot function
  p <- plot_sec_mwd(
    object,
    sample_id = sample_id,
    log_mw = log_scale %in% c("x", "both"),
    ...
  )

  # Apply y log scale if requested
  if (log_scale %in% c("y", "both")) {
    p <- p + ggplot2::scale_y_log10()
  }

  p
}

#' Internal: Conformation plot for autoplot
#' @noRd
autoplot_conformation <- function(object, sample_id = NULL, ...) {
  measure_cols <- attr(object, "measure_cols") %||% find_measure_cols(object)

  # Determine conformation plot type
  if ("rg" %in% measure_cols && "mw" %in% measure_cols) {
    type <- "rg_mw"
  } else if ("intrinsic_visc" %in% measure_cols && "mw" %in% measure_cols) {
    type <- "eta_mw"
  } else if ("rh" %in% measure_cols && "mw" %in% measure_cols) {
    type <- "rh_mw"
  } else {
    cli::cli_abort(c(
      "Cannot create conformation plot.",
      "i" = "Need both MW and conformation data (rg, intrinsic_visc, or rh).",
      "i" = "Available measures: {.val {measure_cols}}"
    ))
  }

  plot_sec_conformation(
    object,
    type = type,
    sample_id = sample_id,
    ...
  )
}

#' Internal: Composition plot for autoplot
#' @noRd
autoplot_composition <- function(object, sample_id = NULL, ...) {
  measure_cols <- attr(object, "measure_cols") %||% find_measure_cols(object)

  # Check for composition column
  comp_col <- NULL
  if ("composition" %in% measure_cols) {
    comp_col <- "composition"
  } else if ("uv_ri_ratio" %in% measure_cols) {
    comp_col <- "uv_ri_ratio"
  }

  if (is.null(comp_col)) {
    cli::cli_abort(c(
      "Cannot create composition plot.",
      "i" = "Need composition or uv_ri_ratio measure column.",
      "i" = "Available measures: {.val {measure_cols}}"
    ))
  }

  # Use chromatogram plot for composition
  plot_sec_chromatogram(
    object,
    measures = comp_col,
    sample_id = sample_id,
    y_label = if (comp_col == "composition") "Composition" else "UV/RI Ratio",
    ...
  )
}

#' Plot SEC Chromatogram
#'
#' Creates a chromatogram plot from SEC data showing detector signal vs
#' elution time/volume.
#'
#' @param data A data frame containing SEC results with measure columns,
#'   or a tibble from [measure_sec_slice_table()].
#' @param measures Character vector of measure column names to plot.
#'   If `NULL`, plots all measure columns found.
#' @param sample_id Column name containing sample identifiers. If `NULL`,
#'   attempts to auto-detect or uses row numbers.
#' @param x_label Label for x-axis. Default is "Elution Time (min)".
#' @param y_label Label for y-axis. Default is "Signal".
#' @param normalize Logical. Normalize signals to 0-1 range for comparison?
#'   Default is `FALSE`.
#' @param facet_by How to facet the plot. One of:
#'   - `"none"`: All on single plot (default)
#'   - `"measure"`: Separate panel per detector/measure
#'   - `"sample"`: Separate panel per sample
#' @param color_by What to map to color aesthetic. One of `"sample"` (default)
#'   or `"measure"`.
#' @param ... Additional arguments passed to [ggplot2::geom_line()].
#'
#' @return A ggplot2 object.
#'
#' @details
#' This is the fundamental SEC visualization showing raw or processed
#' chromatographic data. Works with both:
#' - Processed recipe output (data frames with measure_list columns)
#' - Slice tables from [measure_sec_slice_table()]
#'
#' @family sec-visualization
#' @export
#'
#' @examples
#' \dontrun{
#' library(recipes)
#' library(measure)
#'
#' # Process SEC data
#' processed <- recipe(~., data = sec_triple_detect) |>
#'   step_measure_input_long(
#'     ri_signal,
#'     location = vars(elution_time),
#'     col_name = "ri"
#'   ) |>
#'   step_sec_baseline(measures = "ri") |>
#'   prep() |>
#'   bake(new_data = NULL)
#'
#' # Basic chromatogram
#' plot_sec_chromatogram(processed, measures = "ri")
#'
#' # Normalized overlay of multiple samples
#' plot_sec_chromatogram(processed, measures = "ri", normalize = TRUE)
#'
#' # Faceted by sample
#' plot_sec_chromatogram(processed, facet_by = "sample")
#' }
plot_sec_chromatogram <- function(
  data,
  measures = NULL,

  sample_id = NULL,
  x_label = "Elution Time (min)",
  y_label = "Signal",
  normalize = FALSE,
  facet_by = c("none", "measure", "sample"),
  color_by = c("sample", "measure"),
  ...
) {
  check_ggplot2_available()

  facet_by <- match.arg(facet_by)
  color_by <- match.arg(color_by)

  # Convert to slice table if needed
  slice_data <- prepare_plot_data(
    data,
    measures = measures,
    sample_id = sample_id
  )

  if (nrow(slice_data) == 0) {
    cli::cli_abort("No data to plot after extraction.")
  }

  # Normalize if requested
  if (normalize) {
    slice_data <- slice_data |>
      dplyr::group_by(.data$sample_id, .data$measure) |>
      dplyr::mutate(
        value = (.data$value - min(.data$value, na.rm = TRUE)) /
          (max(.data$value, na.rm = TRUE) -
            min(.data$value, na.rm = TRUE) +
            1e-10)
      ) |>
      dplyr::ungroup()
    y_label <- paste0(y_label, " (normalized)")
  }

  # Build the plot
  if (color_by == "sample") {
    p <- ggplot2::ggplot(
      slice_data,
      ggplot2::aes(
        x = .data$location,
        y = .data$value,
        color = .data$sample_id,
        group = interaction(.data$sample_id, .data$measure)
      )
    )
  } else {
    p <- ggplot2::ggplot(
      slice_data,
      ggplot2::aes(
        x = .data$location,
        y = .data$value,
        color = .data$measure,
        group = interaction(.data$sample_id, .data$measure)
      )
    )
  }

  p <- p +
    ggplot2::geom_line(...) +
    ggplot2::labs(x = x_label, y = y_label) +
    ggplot2::theme_minimal()

  # Add faceting
  if (facet_by == "measure") {
    p <- p + ggplot2::facet_wrap(~measure, scales = "free_y")
  } else if (facet_by == "sample") {
    p <- p + ggplot2::facet_wrap(~sample_id, scales = "free_y")
  }

  p
}


#' Plot Molecular Weight Distribution
#'
#' Creates a molecular weight distribution (MWD) plot showing differential
#' or cumulative distribution.
#'
#' @param data A data frame containing SEC results with calibrated MW data.
#' @param mw_col Name of the molecular weight measure column. Default is `"mw"`.
#' @param concentration_col Name of the concentration/signal column used for
#'   weighting. If `NULL`, uses the first available measure column.
#' @param sample_id Column name containing sample identifiers.
#' @param type Type of distribution to plot:
#'   - `"differential"`: dW/d(log M) vs log M (default)
#'   - `"cumulative"`: Cumulative weight fraction vs log M
#'   - `"both"`: Both on faceted plot
#' @param show_averages Logical. Show vertical lines for Mn, Mw, Mz?
#'   Default is `TRUE`. Requires these columns to be present in data.
#' @param log_mw Logical. Use log10(MW) on x-axis? Default is `TRUE`.
#' @param x_label Label for x-axis. Default is auto-generated based on `log_mw`.
#' @param y_label Label for y-axis. Default is auto-generated based on `type`.
#' @param ... Additional arguments passed to [ggplot2::geom_line()].
#'
#' @return A ggplot2 object.
#'
#' @details
#' The MWD plot is the standard way to visualize polymer molecular weight
#' distributions. The differential distribution (dW/d log M) shows the
#' weight fraction of polymer at each molecular weight.
#'
#' When `show_averages = TRUE`, vertical dashed lines are added for:
#' - Mn (number-average): leftmost
#' - Mw (weight-average): middle
#' - Mz (z-average): rightmost
#'
#' @family sec-visualization
#' @export
#'
#' @examples
#' \dontrun{
#' # After SEC processing with calibration
#' plot_sec_mwd(processed_data)
#'
#' # Cumulative distribution
#' plot_sec_mwd(processed_data, type = "cumulative")
#'
#' # Without MW average lines
#' plot_sec_mwd(processed_data, show_averages = FALSE)
#' }
plot_sec_mwd <- function(
  data,
  mw_col = "mw",
  concentration_col = NULL,
  sample_id = NULL,
  type = c("differential", "cumulative", "both"),
  show_averages = TRUE,
  log_mw = TRUE,
  x_label = NULL,
  y_label = NULL,
  ...
) {
  check_ggplot2_available()

  type <- match.arg(type)

  # Set default axis labels
  if (is.null(x_label)) {
    x_label <- if (log_mw) expression(log[10](M)) else "Molecular Weight (Da)"
  }
  if (is.null(y_label)) {
    y_label <- switch(
      type,
      "differential" = expression(dW / d ~ log ~ M),
      "cumulative" = "Cumulative Weight Fraction",
      "both" = "Distribution"
    )
  }

  # Check for MW column
  if (!mw_col %in% names(data)) {
    cli::cli_abort(
      "Molecular weight column {.val {mw_col}} not found in data."
    )
  }

  # Extract slice data
  measures_to_get <- mw_col
  if (!is.null(concentration_col) && concentration_col %in% names(data)) {
    measures_to_get <- c(measures_to_get, concentration_col)
  }

  slice_data <- prepare_plot_data(
    data,
    measures = measures_to_get,
    sample_id = sample_id
  )

  # Pivot to wide if we have multiple measures
  if (length(measures_to_get) > 1) {
    slice_data <- tidyr::pivot_wider(
      slice_data,
      names_from = "measure",
      values_from = "value"
    )
  } else {
    # Rename value to mw_col name for consistency
    names(slice_data)[names(slice_data) == "value"] <- mw_col
    slice_data$measure <- NULL
  }

  # Calculate distribution
  slice_data <- slice_data |>
    dplyr::filter(!is.na(.data[[mw_col]]) & .data[[mw_col]] > 0) |>
    dplyr::mutate(log_mw = log10(.data[[mw_col]]))

  # If no concentration column, use uniform weighting
  if (is.null(concentration_col) || !concentration_col %in% names(slice_data)) {
    slice_data$weight <- 1
  } else {
    slice_data$weight <- slice_data[[concentration_col]]
  }

  # Calculate differential distribution per sample
  slice_data <- slice_data |>
    dplyr::group_by(.data$sample_id) |>
    dplyr::arrange(.data$log_mw) |>
    dplyr::mutate(
      # Normalize weights
      weight_norm = .data$weight / sum(.data$weight, na.rm = TRUE),
      # Cumulative
      cumulative = cumsum(.data$weight_norm),
      # Differential (dW/d log M approximation)
      d_log_mw = c(diff(.data$log_mw), NA),
      differential = .data$weight_norm / abs(.data$d_log_mw)
    ) |>
    dplyr::ungroup()

  # Set x values based on log_mw preference
  if (log_mw) {
    slice_data$x_val <- slice_data$log_mw
  } else {
    slice_data$x_val <- slice_data[[mw_col]]
  }

  # Build plot based on type
  if (type == "differential") {
    p <- ggplot2::ggplot(
      slice_data,
      ggplot2::aes(
        x = .data$x_val,
        y = .data$differential,
        color = .data$sample_id
      )
    ) +
      ggplot2::geom_line(...)
  } else if (type == "cumulative") {
    p <- ggplot2::ggplot(
      slice_data,
      ggplot2::aes(
        x = .data$x_val,
        y = .data$cumulative,
        color = .data$sample_id
      )
    ) +
      ggplot2::geom_line(...)
  } else {
    # Both - pivot longer and facet
    plot_data <- slice_data |>
      tidyr::pivot_longer(
        cols = c("differential", "cumulative"),
        names_to = "dist_type",
        values_to = "y_val"
      )

    p <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(
        x = .data$x_val,
        y = .data$y_val,
        color = .data$sample_id
      )
    ) +
      ggplot2::geom_line(...) +
      ggplot2::facet_wrap(~dist_type, scales = "free_y")
  }

  p <- p +
    ggplot2::labs(x = x_label, y = y_label, color = "Sample") +
    ggplot2::theme_minimal()

  # Add MW average lines if available and requested
  if (show_averages && type != "both") {
    p <- add_mw_average_lines(p, data, sample_id, log_mw)
  }

  p
}


#' Plot Multi-Detector SEC Overlay
#'
#' Creates an overlay plot of multiple SEC detectors, optionally normalized
#' and aligned.
#'
#' @param data A data frame containing SEC results with multiple detector
#'   measure columns.
#' @param detectors Character vector of detector column names to include.
#'   Common values: `c("ri", "uv", "mals", "visc")`.
#' @param sample_id Column name containing sample identifiers. If `NULL`,
#'   plots all samples or auto-detects.
#' @param samples Character vector of specific sample IDs to plot.
#'   If `NULL`, plots all samples.
#' @param normalize Logical. Normalize each detector to 0-1 range for
#'   comparison? Default is `TRUE`.
#' @param x_label Label for x-axis. Default is "Elution Time (min)".
#' @param facet Logical. Create separate panel for each sample?
#'   Default is `FALSE` (overlay).
#' @param ... Additional arguments passed to [ggplot2::geom_line()].
#'
#' @return A ggplot2 object.
#'
#' @details
#' Multi-detector overlay plots are essential for:
#' - Verifying detector alignment after delay correction
#' - Identifying composition drift in copolymers (UV/RI differences)
#' - Detecting aggregates (MALS response higher than expected from RI)
#' - Checking for baseline issues across detectors
#'
#' When `normalize = TRUE` (default), each detector signal is scaled to 0-1
#' range, making it easy to compare peak shapes and positions across
#' detectors with very different response magnitudes.
#'
#' @family sec-visualization
#' @export
#'
#' @examples
#' \dontrun{
#' # Multi-detector overlay for single sample
#' plot_sec_multidetector(
#'   processed_data,
#'   detectors = c("ri", "uv", "mals"),
#'   samples = "PMMA-1"
#' )
#'
#' # Faceted by sample
#' plot_sec_multidetector(
#'   processed_data,
#'   detectors = c("ri", "uv"),
#'   facet = TRUE
#' )
#' }
plot_sec_multidetector <- function(
  data,
  detectors,
  sample_id = NULL,
  samples = NULL,
  normalize = TRUE,
  x_label = "Elution Time (min)",
  facet = FALSE,
  ...
) {
  check_ggplot2_available()

  if (missing(detectors) || length(detectors) == 0) {
    cli::cli_abort(
      "{.arg detectors} must specify at least one detector column."
    )
  }

  # Check detectors exist
  available_measures <- find_measure_cols(data)
  missing_det <- setdiff(detectors, available_measures)
  if (length(missing_det) > 0) {
    cli::cli_warn(
      "Detector columns not found: {.val {missing_det}}. Skipping."
    )
    detectors <- intersect(detectors, available_measures)
    if (length(detectors) == 0) {
      cli::cli_abort("No valid detector columns found.")
    }
  }

  # Extract slice data
  slice_data <- prepare_plot_data(
    data,
    measures = detectors,
    sample_id = sample_id
  )

  # Filter to specific samples if requested
  if (!is.null(samples)) {
    slice_data <- dplyr::filter(slice_data, .data$sample_id %in% samples)
    if (nrow(slice_data) == 0) {
      cli::cli_abort("No data found for samples: {.val {samples}}")
    }
  }

  # Normalize if requested
  if (normalize) {
    slice_data <- slice_data |>
      dplyr::group_by(.data$sample_id, .data$measure) |>
      dplyr::mutate(
        value = (.data$value - min(.data$value, na.rm = TRUE)) /
          (max(.data$value, na.rm = TRUE) -
            min(.data$value, na.rm = TRUE) +
            1e-10)
      ) |>
      dplyr::ungroup()
  }

  # Create prettier detector labels
  detector_labels <- c(
    "ri" = "RI",
    "uv" = "UV",
    "mals" = "MALS",
    "visc" = "Viscometer",
    "viscometer" = "Viscometer",
    "dad" = "DAD",
    "lals" = "LALS",
    "rals" = "RALS",
    "dls" = "DLS"
  )

  slice_data <- slice_data |>
    dplyr::mutate(
      detector = dplyr::if_else(
        .data$measure %in% names(detector_labels),
        detector_labels[.data$measure],
        .data$measure
      )
    )

  # Build plot
  p <- ggplot2::ggplot(
    slice_data,
    ggplot2::aes(
      x = .data$location,
      y = .data$value,
      color = .data$detector,
      group = interaction(.data$sample_id, .data$measure)
    )
  ) +
    ggplot2::geom_line(...) +
    ggplot2::labs(
      x = x_label,
      y = if (normalize) "Normalized Signal" else "Signal",
      color = "Detector"
    ) +
    ggplot2::theme_minimal()

  if (facet) {
    p <- p + ggplot2::facet_wrap(~sample_id)
  }

  p
}


#' Plot SEC Conformation Data
#'
#' Creates conformation plots showing structure-MW relationships from
#' multi-detector SEC data (Mark-Houwink plots, Rg-MW scaling).
#'
#' @param data A data frame containing SEC results with MW and conformation
#'   data (Rg, intrinsic viscosity).
#' @param type Type of conformation plot:
#'   - `"rg_mw"`: Radius of gyration vs molecular weight (default)
#'   - `"eta_mw"`: Intrinsic viscosity vs molecular weight
#'   - `"rh_mw"`: Hydrodynamic radius vs molecular weight
#' @param mw_col Name of molecular weight column. Default is `"mw"`.
#' @param y_col Name of y-axis column. If `NULL`, auto-detected based on type.
#' @param sample_id Column name containing sample identifiers.
#' @param show_fit Logical. Show power-law fit line? Default is `TRUE`.
#' @param show_exponent Logical. Annotate slope/exponent on plot?
#'   Default is `TRUE`.
#' @param compare_linear Data frame with linear reference polymer for
#'   branching comparison. Should have same structure as main data.
#' @param ... Additional arguments passed to [ggplot2::geom_point()].
#'
#' @return A ggplot2 object.
#'
#' @details
#' Conformation plots reveal polymer architecture:
#'
#' **Rg-MW (radius of gyration):**
#' Log-log slope indicates conformation:
#' - 0.33: Compact sphere

#' - 0.5-0.6: Random coil (linear polymer in good solvent)
#' - 1.0: Rigid rod
#' Branched polymers show reduced Rg at same MW compared to linear.
#'
#' **eta-MW (Mark-Houwink):**
#' The relationship `[eta] = K * M^a` where:
#' - K and a are Mark-Houwink parameters
#' - a is approximately 0.5-0.8 for typical polymers
#' - Branched polymers show lower `[eta]` at same MW
#'
#' @family sec-visualization
#' @export
#'
#' @examples
#' \dontrun{
#' # Rg-MW plot from MALS data
#' plot_sec_conformation(processed_data, type = "rg_mw")
#'
#' # Mark-Houwink plot
#' plot_sec_conformation(processed_data, type = "eta_mw")
#'
#' # Compare branched to linear reference
#' plot_sec_conformation(
#'   branched_data,
#'   compare_linear = linear_reference
#' )
#' }
plot_sec_conformation <- function(
  data,
  type = c("rg_mw", "eta_mw", "rh_mw"),
  mw_col = "mw",
  y_col = NULL,
  sample_id = NULL,
  show_fit = TRUE,
  show_exponent = TRUE,
  compare_linear = NULL,
  ...
) {
  check_ggplot2_available()

  type <- match.arg(type)

  # Auto-detect y column based on type
  if (is.null(y_col)) {
    y_col <- switch(
      type,
      "rg_mw" = "rg",
      "eta_mw" = "intrinsic_visc",
      "rh_mw" = "rh"
    )
  }

  # Set axis labels
  axis_labels <- list(
    "rg_mw" = list(
      x = expression(log[10](M)),
      y = expression(log[10](R[g] ~ "(nm)"))
    ),
    "eta_mw" = list(
      x = expression(log[10](M)),
      y = expression(log[10](group("[", eta, "]") ~ "(mL/g)"))
    ),
    "rh_mw" = list(
      x = expression(log[10](M)),
      y = expression(log[10](R[h] ~ "(nm)"))
    )
  )

  # Check columns exist
  available <- find_measure_cols(data)
  if (!mw_col %in% available) {
    cli::cli_abort("MW column {.val {mw_col}} not found.")
  }
  if (!y_col %in% available) {
    cli::cli_abort("Y column {.val {y_col}} not found for type {.val {type}}.")
  }

  # Extract slice data
  slice_data <- prepare_plot_data(
    data,
    measures = c(mw_col, y_col),
    sample_id = sample_id
  )

  # Pivot to wide format
  slice_data <- tidyr::pivot_wider(
    slice_data,
    names_from = "measure",
    values_from = "value"
  )

  # Filter valid data and compute log values
  slice_data <- slice_data |>
    dplyr::filter(
      !is.na(.data[[mw_col]]) &
        .data[[mw_col]] > 0 &
        !is.na(.data[[y_col]]) &
        .data[[y_col]] > 0
    ) |>
    dplyr::mutate(
      log_mw = log10(.data[[mw_col]]),
      log_y = log10(.data[[y_col]]),
      source = "Sample"
    )

  # Add linear reference if provided
  if (!is.null(compare_linear)) {
    ref_data <- prepare_plot_data(
      compare_linear,
      measures = c(mw_col, y_col),
      sample_id = sample_id
    )
    ref_data <- tidyr::pivot_wider(
      ref_data,
      names_from = "measure",
      values_from = "value"
    )
    ref_data <- ref_data |>
      dplyr::filter(
        !is.na(.data[[mw_col]]) &
          .data[[mw_col]] > 0 &
          !is.na(.data[[y_col]]) &
          .data[[y_col]] > 0
      ) |>
      dplyr::mutate(
        log_mw = log10(.data[[mw_col]]),
        log_y = log10(.data[[y_col]]),
        source = "Linear Reference"
      )
    slice_data <- dplyr::bind_rows(slice_data, ref_data)
  }

  # Build plot
  p <- ggplot2::ggplot(
    slice_data,
    ggplot2::aes(x = .data$log_mw, y = .data$log_y, color = .data$sample_id)
  ) +
    ggplot2::geom_point(alpha = 0.5, ...) +
    ggplot2::labs(
      x = axis_labels[[type]]$x,
      y = axis_labels[[type]]$y,
      color = "Sample"
    ) +
    ggplot2::theme_minimal()

  # Add fit line
  if (show_fit) {
    p <- p +
      ggplot2::geom_smooth(
        method = "lm",
        formula = y ~ x,
        se = FALSE,
        linetype = "dashed"
      )

    # Add exponent annotation
    if (show_exponent) {
      # Calculate slope for annotation with error handling
      fit_data <- slice_data |>
        dplyr::group_by(.data$sample_id) |>
        dplyr::summarise(
          slope = tryCatch(
            stats::coef(stats::lm(.data$log_y ~ .data$log_mw))[2],
            error = function(e) NA_real_
          ),
          max_x = max(.data$log_mw, na.rm = TRUE),
          max_y = max(.data$log_y, na.rm = TRUE),
          .groups = "drop"
        ) |>
        dplyr::filter(!is.na(.data$slope))

      p <- p +
        ggplot2::geom_text(
          data = fit_data,
          ggplot2::aes(
            x = .data$max_x,
            y = .data$max_y,
            label = sprintf("slope = %.2f", .data$slope)
          ),
          hjust = 1,
          vjust = 0,
          size = 3
        )
    }
  }

  # Style for linear reference comparison
  if (!is.null(compare_linear)) {
    p <- p +
      ggplot2::aes(shape = .data$source) +
      ggplot2::scale_shape_manual(
        values = c("Sample" = 16, "Linear Reference" = 1)
      )
  }

  p
}


#' Plot SEC Calibration Curve
#'
#' Creates a calibration curve plot from SEC calibration data or a
#' prepped recipe with conventional calibration.
#'
#' @param data Either:
#'   - A data frame of calibration standards with retention and log_mw columns
#'   - A prepped recipe containing a conventional calibration step
#' @param retention_col Name of retention time/volume column.
#'   Default is `"retention_time"`.
#' @param mw_col Name of log MW column. Default is `"log_mp"`.
#' @param show_residuals Logical. Show residual plot below? Default is `FALSE`.
#' @param show_r_squared Logical. Show R-squared value? Default is `TRUE`.
#' @param ... Additional arguments passed to [ggplot2::geom_point()].
#'
#' @return A ggplot2 object (or patchwork of plots if show_residuals = TRUE).
#'
#' @family sec-visualization
#' @export
#'
#' @examples
#' \dontrun{
#' # Plot calibration from standards data
#' plot_sec_calibration(sec_ps_standards)
#'
#' # With residuals panel
#' plot_sec_calibration(sec_ps_standards, show_residuals = TRUE)
#' }
plot_sec_calibration <- function(
  data,
  retention_col = "retention_time",
  mw_col = "log_mp",
  show_residuals = FALSE,
  show_r_squared = TRUE,
  ...
) {
  check_ggplot2_available()

  # Handle recipe input
  if (inherits(data, "recipe")) {
    cli::cli_abort(
      "Recipe input not yet supported. Please provide standards data frame."
    )
  }

  # Check columns
  if (!retention_col %in% names(data)) {
    cli::cli_abort("Retention column {.val {retention_col}} not found.")
  }
  if (!mw_col %in% names(data)) {
    cli::cli_abort("MW column {.val {mw_col}} not found.")
  }

  # Validate sufficient data points for polynomial fit
  n_standards <- nrow(data)
  if (n_standards < 4) {
    cli::cli_abort(c(
      "Insufficient data for calibration curve.",
      "x" = "Found {n_standards} standard{?s}, need at least 4 for cubic fit.",
      "i" = "Provide more calibration standards or use a lower polynomial."
    ))
  }

  # Fit calibration
  fit <- stats::lm(
    stats::as.formula(paste(mw_col, "~ poly(", retention_col, ", 3)")),
    data = data
  )

  # Add predictions
  data$predicted <- stats::predict(fit)
  data$residual <- data[[mw_col]] - data$predicted

  # Main calibration plot
  p_main <- ggplot2::ggplot(
    data,
    ggplot2::aes(x = .data[[retention_col]], y = .data[[mw_col]])
  ) +
    ggplot2::geom_point(size = 3, ...) +
    ggplot2::geom_line(
      ggplot2::aes(y = .data$predicted),
      color = "blue",
      linewidth = 1
    ) +
    ggplot2::labs(
      x = "Retention Time (min)",
      y = expression(log[10](M[p])),
      title = "SEC Calibration Curve"
    ) +
    ggplot2::theme_minimal()

  # Add equation and R-squared annotations
  r_sq <- summary(fit)$r.squared
  annotations <- character()

  if (show_r_squared) {
    annotations <- c(annotations, sprintf("R\u00B2 = %.5f", r_sq))
  }

  if (length(annotations) > 0) {
    p_main <- p_main +
      ggplot2::annotate(
        "text",
        x = max(data[[retention_col]]),
        y = max(data[[mw_col]]),
        label = paste(annotations, collapse = "\n"),
        hjust = 1,
        vjust = 1,
        size = 3.5
      )
  }

  if (show_residuals) {
    p_resid <- ggplot2::ggplot(
      data,
      ggplot2::aes(x = .data[[retention_col]], y = .data$residual)
    ) +
      ggplot2::geom_hline(
        yintercept = 0,
        linetype = "dashed",
        color = "gray50"
      ) +
      ggplot2::geom_point(size = 2) +
      ggplot2::labs(x = "Retention Time (min)", y = "Residual (log M)") +
      ggplot2::theme_minimal()

    # Combine with patchwork if available
    if (requireNamespace("patchwork", quietly = TRUE)) {
      return(p_main / p_resid + patchwork::plot_layout(heights = c(3, 1)))
    } else {
      cli::cli_warn(c(
        "Install {.pkg patchwork} for combined residual plots.",
        "i" = "Returning main plot only.",
        "i" = "Install with {.code install.packages(\"patchwork\")}"
      ))
      return(p_main)
    }
  }

  p_main
}


#' Quick SEC Plot
#'
#' Creates an automatic plot appropriate for SEC analysis results.
#' Dispatches to the most appropriate plot type based on available data.
#'
#' @param data A data frame containing SEC results with measure columns.
#' @param type Type of plot to create. One of:
#'   - `"auto"`: Automatically detect (default)
#'   - `"chromatogram"`: Basic chromatogram
#'   - `"mwd"`: Molecular weight distribution
#'   - `"multidetector"`: Multi-detector overlay
#'   - `"conformation"`: Rg-MW or eta-MW plot
#' @param ... Additional arguments passed to the underlying plot function.
#'
#' @return A ggplot2 object.
#'
#' @details
#' When `type = "auto"` (default), the function chooses the plot type based on
#' available data:
#' - If MW column present: MWD plot
#' - If multiple detectors: Multi-detector overlay
#' - Otherwise: Basic chromatogram
#'
#' This is a convenience wrapper that dispatches to the specific plot
#' functions: [plot_sec_chromatogram()], [plot_sec_mwd()],
#' [plot_sec_multidetector()], or [plot_sec_conformation()].
#'
#' @family sec-visualization
#' @export
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#'
#' # Auto-detect plot type
#' plot_sec(processed_sec_data)
#'
#' # Specific plot type
#' plot_sec(processed_sec_data, type = "mwd")
#' }
plot_sec <- function(
  data,
  type = c("auto", "chromatogram", "mwd", "multidetector", "conformation"),
  ...
) {
  check_ggplot2_available()

  type <- match.arg(type)

  # Find available measure columns
  measure_cols <- find_measure_cols(data)

  if (length(measure_cols) == 0) {
    cli::cli_abort("No measure columns found in data for plotting.")
  }

  # Auto-detect appropriate plot type
  if (type == "auto") {
    if ("mw" %in% measure_cols) {
      type <- "mwd"
    } else if (length(measure_cols) >= 2) {
      type <- "multidetector"
    } else {
      type <- "chromatogram"
    }
  }

  # Dispatch to appropriate function
  switch(
    type,
    "chromatogram" = plot_sec_chromatogram(data, ...),
    "mwd" = plot_sec_mwd(data, ...),
    "multidetector" = plot_sec_multidetector(
      data,
      detectors = measure_cols,
      ...
    ),
    "conformation" = plot_sec_conformation(data, ...),
    cli::cli_abort("Unknown plot type: {.val {type}}")
  )
}


# ==============================================================================
# Helper Functions
# ==============================================================================

#' Check if ggplot2 is available
#' @noRd
check_ggplot2_available <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    cli::cli_abort(
      c(
        "Package {.pkg ggplot2} is required for plotting.",
        "i" = "Install it with {.code install.packages(\"ggplot2\")}"
      )
    )
  }
}

#' Find measure columns in a data frame
#' @noRd
detect_measure_cols <- function(data) {
  # Use measure package function if available
  if (requireNamespace("measure", quietly = TRUE)) {
    return(measure::find_measure_cols(data))
  }

  # Warn about fallback behavior
  cli::cli_warn(c(
    "Package {.pkg measure} is not available.",
    "i" = "Using fallback detection which may produce different results.",
    "i" = "Install {.pkg measure} with {.code install.packages(\"measure\")}"
  ))

  # Fallback: look for list columns with specific structure
  list_cols <- names(data)[vapply(data, is.list, logical(1))]

  # Check each list column for measure_tbl structure
  measure_cols <- character()
  for (col in list_cols) {
    if (length(data[[col]]) == 0) {
      next
    }
    first_elem <- data[[col]][[1]]
    if (is.null(first_elem)) {
      next
    }
    if (
      is.list(first_elem) &&
        all(c("location", "value") %in% names(first_elem))
    ) {
      measure_cols <- c(measure_cols, col)
    }
  }

  measure_cols
}

# Alias for backwards compatibility within package
find_measure_cols <- detect_measure_cols

#' Prepare data for plotting
#'
#' Converts data with measure columns to long-format slice table suitable
#' for ggplot2.
#' @noRd
prepare_plot_data <- function(data, measures = NULL, sample_id = NULL) {
  # If already in slice table format, return as-is
  if (all(c("sample_id", "location", "value") %in% names(data))) {
    if (!is.null(measures) && "measure" %in% names(data)) {
      data <- dplyr::filter(data, .data$measure %in% measures)
    }
    return(data)
  }

  # Otherwise, extract from measure columns
  measure_sec_slice_table(
    data,
    measures = measures,
    sample_id = sample_id,
    include_location = TRUE,
    pivot = FALSE
  )
}

#' Add MW average vertical lines to a plot
#' @noRd
add_mw_average_lines <- function(p, data, sample_id, log_mw) {
  # Look for MW average columns
  mw_cols <- c("mw_mn", "mw_mw", "mw_mz", "Mn", "Mw", "Mz")
  available <- intersect(mw_cols, names(data))

  if (length(available) == 0) {
    cli::cli_warn(c(
      "Cannot add MW average lines.",
      "i" = "No MW average columns found. Expected one of: {.val {mw_cols}}"
    ))
    return(p)
  }

  # Determine sample ID column
  if (!is.null(sample_id) && sample_id %in% names(data)) {
    id_col <- sample_id
  } else {
    id_col <- NULL
  }

  # Create data for vertical lines
  line_data <- list()
  line_colors <- c(
    "mw_mn" = "#E69F00",
    "Mn" = "#E69F00",
    "mw_mw" = "#56B4E9",
    "Mw" = "#56B4E9",
    "mw_mz" = "#009E73",
    "Mz" = "#009E73"
  )
  line_labels <- c(
    "mw_mn" = "Mn",
    "Mn" = "Mn",
    "mw_mw" = "Mw",
    "Mw" = "Mw",
    "mw_mz" = "Mz",
    "Mz" = "Mz"
  )

  for (col in available) {
    vals <- data[[col]]
    if (all(is.na(vals))) {
      next
    }

    if (log_mw) {
      vals <- log10(vals)
    }

    for (i in seq_along(vals)) {
      if (!is.na(vals[i])) {
        line_data <- c(
          line_data,
          list(
            data.frame(
              xintercept = vals[i],
              average_type = line_labels[col],
              color = line_colors[col],
              sample_id = if (!is.null(id_col)) data[[id_col]][i] else i
            )
          )
        )
      }
    }
  }

  if (length(line_data) > 0) {
    line_df <- dplyr::bind_rows(line_data)
    p <- p +
      ggplot2::geom_vline(
        data = line_df,
        ggplot2::aes(
          xintercept = .data$xintercept,
          linetype = .data$average_type
        ),
        alpha = 0.7
      ) +
      ggplot2::scale_linetype_manual(
        values = c("Mn" = "dotted", "Mw" = "dashed", "Mz" = "dotdash"),
        name = "MW Average"
      )
  }

  p
}
