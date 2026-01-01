# ==============================================================================
# SEC Quality Control Functions
#
# System suitability, column performance, and quality metrics
# ==============================================================================

#' Calculate Peak Resolution
#'
#' Calculates the resolution between two chromatographic peaks using the
#' USP or EP formula.
#'
#' @param retention_1 Retention time of the first peak (earlier eluting).
#' @param retention_2 Retention time of the second peak (later eluting).
#' @param width_1 Peak width of the first peak. See `width_type` for units.
#' @param width_2 Peak width of the second peak.
#' @param width_type Type of peak width measurement:
#'   \itemize{
#'     \item `"baseline"` (default): Width at baseline (Wb)
#'     \item `"half_height"`: Width at half height (W0.5h)
#'     \item `"tangent"`: Width from tangent lines at inflection points
#'   }
#' @param method Resolution formula to use:
#'   \itemize{
#'     \item `"usp"` (default): Rs = 2(t2 - t1) / (w1 + w2)
#'     \item `"ep"`: Rs = 1.18(t2 - t1) / (w1_0.5h + w2_0.5h)
#'   }
#'
#' @return Numeric resolution value. Rs > 1.5 indicates baseline separation.
#'
#' @details
#' Resolution quantifies the degree of separation between adjacent peaks:
#'
#' **USP Formula (baseline width):**
#' \deqn{R_s = \frac{2(t_2 - t_1)}{W_{b1} + W_{b2}}}
#'
#' **EP Formula (half-height width):**
#' \deqn{R_s = \frac{1.18(t_2 - t_1)}{W_{0.5h,1} + W_{0.5h,2}}}
#'
#' **Interpretation:**
#' \itemize{
#'   \item Rs < 1.0: Peaks overlap significantly
#'   \item Rs = 1.0: ~94% separation (4 sigma)
#'   \item Rs = 1.5: Baseline separation (~99.7%)
#'   \item Rs > 2.0: Complete separation with gap
#' }
#'
#' @family sec-qc
#' @export
#'
#' @examples
#' # Calculate resolution between monomer and dimer
#' measure_sec_resolution(
#'   retention_1 = 8.2,   # dimer (elutes first in SEC)
#'   retention_2 = 9.5,   # monomer
#'   width_1 = 0.4,
#'   width_2 = 0.5
#' )
measure_sec_resolution <- function(
    retention_1,
    retention_2,
    width_1,
    width_2,
    width_type = c("baseline", "half_height", "tangent"),
    method = c("usp", "ep")
) {
  width_type <- match.arg(width_type)
  method <- match.arg(method)

 # Validate inputs
  if (!is.numeric(retention_1) || !is.numeric(retention_2)) {
    cli::cli_abort("Retention times must be numeric.")
  }
  if (!is.numeric(width_1) || !is.numeric(width_2)) {
    cli::cli_abort("Peak widths must be numeric.")
  }
  if (width_1 <= 0 || width_2 <= 0) {
    cli::cli_abort("Peak widths must be positive.")
  }

  # Ensure proper order (peak 1 should elute before peak 2)
  if (retention_1 > retention_2) {
    temp <- retention_1
    retention_1 <- retention_2
    retention_2 <- temp
    temp <- width_1
    width_1 <- width_2
    width_2 <- temp
  }

  delta_t <- retention_2 - retention_1

  if (method == "usp") {
    # USP: Rs = 2(t2 - t1) / (w1 + w2)
    # For baseline width
    Rs <- 2 * delta_t / (width_1 + width_2)
  } else {
    # EP: Rs = 1.18(t2 - t1) / (w1 + w2)
    # For half-height width
    Rs <- 1.18 * delta_t / (width_1 + width_2)
  }

  Rs
}


#' Calculate Theoretical Plate Count
#'
#' Calculates the number of theoretical plates (N) for a chromatographic peak,
#' a measure of column efficiency.
#'
#' @param retention Retention time of the peak.
#' @param width Peak width. See `width_type` for measurement method.
#' @param width_type Type of peak width measurement:
#'   \itemize{
#'     \item `"half_height"` (default): Width at 50% height (W0.5h)
#'     \item `"baseline"`: Width at baseline from tangent lines (Wb)
#'     \item `"inflection"`: Width at inflection points (Wi)
#'   }
#' @param dead_time Column dead time (t0). If provided, calculates effective
#'   plates (N_eff) using adjusted retention time.
#'
#' @return Numeric plate count. Higher values indicate better efficiency.
#'
#' @details
#' Theoretical plate count measures column efficiency:
#'
#' **Half-height width (most common):**
#' \deqn{N = 5.54 \left(\frac{t_R}{W_{0.5h}}\right)^2}
#'
#' **Baseline width:**
#' \deqn{N = 16 \left(\frac{t_R}{W_b}\right)^2}
#'
#' **With dead time correction (effective plates):**
#' \deqn{N_{eff} = 5.54 \left(\frac{t_R - t_0}{W_{0.5h}}\right)^2}
#'
#' **Typical SEC Performance:**
#' \itemize{
#'   \item Analytical SEC columns: 10,000-40,000 plates/meter
#'   \item Preparative columns: 5,000-15,000 plates/meter
#'   \item UHPLC SEC: 50,000+ plates/meter
#' }
#'
#' @family sec-qc
#' @export
#'
#' @examples
#' # Calculate plate count for a monomer peak
#' measure_sec_plate_count(
#'   retention = 9.5,
#'   width = 0.25,
#'   width_type = "half_height"
#' )
#'
#' # With dead time for effective plates
#' measure_sec_plate_count(
#'   retention = 9.5,
#'   width = 0.25,
#'   dead_time = 3.0
#' )
measure_sec_plate_count <- function(
    retention,
    width,
    width_type = c("half_height", "baseline", "inflection"),
    dead_time = NULL
) {
  width_type <- match.arg(width_type)

  if (!is.numeric(retention) || retention <= 0) {
    cli::cli_abort("Retention time must be a positive number.")
  }
  if (!is.numeric(width) || width <= 0) {
    cli::cli_abort("Peak width must be a positive number.")
  }

  # Adjust retention if dead time provided
  if (!is.null(dead_time)) {
    if (!is.numeric(dead_time) || dead_time < 0) {
      cli::cli_abort("Dead time must be a non-negative number.")
    }
    if (dead_time >= retention) {
      cli::cli_abort("Dead time must be less than retention time.")
    }
    retention <- retention - dead_time
  }

  # Calculate N based on width type
  ratio <- retention / width

  N <- switch(
    width_type,
    "half_height" = 5.54 * ratio^2,
    "baseline" = 16 * ratio^2,
    "inflection" = 4 * ratio^2
  )

  N
}


#' Calculate Peak Asymmetry Factor
#'
#' Calculates the asymmetry factor (As) or tailing factor (Tf) for a
#' chromatographic peak.
#'
#' @param leading Width of the leading (front) half of the peak at the
#'   measurement height.
#' @param tailing Width of the tailing (back) half of the peak at the
#'   measurement height.
#' @param method Asymmetry calculation method:
#'   \itemize{
#'     \item `"usp"` (default): Tailing factor at 5% height
#'     \item `"ep"`: Asymmetry factor at 10% height
#'   }
#'
#' @return Numeric asymmetry value. Values > 1 indicate tailing, < 1 indicate
#'   fronting. Ideal value is 1.0.
#'
#' @details
#' Peak asymmetry indicates deviation from ideal Gaussian peak shape:
#'
#' **USP Tailing Factor (at 5% height):**
#' \deqn{T_f = \frac{W_{0.05}}{2f}}
#'
#' where W_0.05 is the width at 5% height and f is the leading half-width.
#'
#' **EP Asymmetry Factor (at 10% height):**
#' \deqn{A_s = \frac{b}{a}}
#'
#' where b is the tailing half-width and a is the leading half-width.
#'
#' **Interpretation:**
#' \itemize{
#'   \item As = 1.0: Symmetric (ideal)
#'   \item As < 0.9 or > 1.2: Slight asymmetry (acceptable)
#'   \item As < 0.8 or > 1.5: Significant asymmetry (investigate)
#'   \item As > 2.0: Severe tailing (column/sample issue)
#' }
#'
#' @family sec-qc
#' @export
#'
#' @examples
#' # Calculate USP tailing factor
#' measure_sec_asymmetry(
#'   leading = 0.12,
#'   tailing = 0.15,
#'   method = "usp"
#' )
measure_sec_asymmetry <- function(
    leading,
    tailing,
    method = c("usp", "ep")
) {
  method <- match.arg(method)

  if (!is.numeric(leading) || leading <= 0) {
    cli::cli_abort("Leading width must be a positive number.")
  }
  if (!is.numeric(tailing) || tailing <= 0) {
    cli::cli_abort("Tailing width must be a positive number.")
  }

  if (method == "usp") {
    # USP Tailing Factor: Tf = (a + b) / 2a
    # where a = leading half-width, b = tailing half-width
    Tf <- (leading + tailing) / (2 * leading)
    return(Tf)
  } else {
    # EP Asymmetry Factor: As = b / a
    As <- tailing / leading
    return(As)
  }
}


#' Calculate Mass Recovery
#'
#' Calculates the mass recovery (percent of injected mass detected) for
#' SEC analysis.
#'
#' @param detected_mass Mass detected by integration of the chromatogram.
#' @param injected_mass Mass injected onto the column.
#' @param units Units for mass values. Both must be in the same units.
#'
#' @return Numeric recovery percentage (0-100+).
#'
#' @details
#' Mass recovery verifies that the analytical system is detecting all of the
#' injected sample:
#'
#' \deqn{\% Recovery = \frac{m_{detected}}{m_{injected}} \times 100}
#'
#' **Interpretation:**
#' \itemize{
#'   \item 95-105%: Excellent recovery (typical acceptance)
#'   \item 90-95% or 105-110%: Acceptable (investigate if persistent)
#'   \item < 90%: Low recovery - possible column adsorption, precipitation
#'   \item > 110%: High recovery - calibration issue, interference
#' }
#'
#' **Common Causes of Low Recovery:**
#' \itemize{
#'   \item Sample adsorption to column packing
#'   \item Sample precipitation or aggregation on-column
#'   \item Detector calibration drift
#'   \item Integration baseline errors
#'   \item Sample degradation
#' }
#'
#' @family sec-qc
#' @export
#'
#' @examples
#' # Calculate recovery
#' measure_sec_recovery(
#'   detected_mass = 0.195,
#'   injected_mass = 0.200
#' )
#' # Returns 97.5%
measure_sec_recovery <- function(
    detected_mass,
    injected_mass,
    units = "mg"
) {
  if (!is.numeric(detected_mass) || detected_mass < 0) {
    cli::cli_abort("Detected mass must be a non-negative number.")
  }
  if (!is.numeric(injected_mass) || injected_mass <= 0) {
    cli::cli_abort("Injected mass must be a positive number.")
  }

  recovery <- 100 * detected_mass / injected_mass

  recovery
}


#' System Suitability Test for SEC
#'
#' Performs comprehensive system suitability testing for SEC analysis,
#' evaluating resolution, plate count, asymmetry, and other quality metrics.
#'
#' @param data A data frame or tibble containing chromatogram data with at
#'   minimum retention times and peak parameters.
#' @param peaks A data frame with peak information. Must contain columns:
#'   \itemize{
#'     \item `retention`: Peak retention time
#'     \item `width`: Peak width (at half height unless specified)
#'     \item `area`: Peak area (for recovery calculation)
#'     \item `height`: Peak height (optional, for asymmetry from raw data
#'   }
#' @param reference_peaks Character vector of peak names to use for resolution
#'   calculation (e.g., c("dimer", "monomer")).
#' @param injected_mass Injected mass for recovery calculation (optional).
#' @param criteria A list of acceptance criteria. Default uses common
#'   biopharmaceutical criteria.
#' @param column_length Column length in cm (for plates per meter calculation).
#'
#' @return A list of class `sec_suitability` containing:
#'   \describe{
#'     \item{results}{Data frame of calculated metrics and pass/fail status}
#'     \item{passed}{Logical indicating if all criteria passed}
#'     \item{summary}{Character summary of results}
#'     \item{criteria}{Criteria used for evaluation}
#'   }
#'
#' @details
#' System suitability testing (SST) verifies that the chromatographic system
#' is performing adequately before, during, and after sample analysis.
#'
#' **Standard SEC SST Parameters:**
#' \itemize{
#'   \item Resolution: Rs >= 1.5 between critical pair
#'   \item Plate count: N >= specified minimum
#'   \item Tailing factor: 0.8 <= Tf <= 1.5
#'   \item Mass recovery: 95-105%
#'   \item Retention time RSD: <= 1.0%
#'   \item Peak area RSD: <= 2.0%
#' }
#'
#' **Regulatory References:**
#' \itemize{
#'   \item USP <621> Chromatography
#'   \item ICH Q2(R1) Validation
#'   \item EP 2.2.46 Chromatographic Separation Techniques
#' }
#'
#' @family sec-qc
#' @export
#'
#' @examples
#' \dontrun{
#' # Define peaks from integration results
#' peaks <- data.frame(
#'   name = c("aggregate", "dimer", "monomer", "fragment"),
#'   retention = c(7.2, 8.5, 9.8, 11.5),
#'   width = c(0.3, 0.25, 0.28, 0.35),
#'   area = c(2.1, 5.3, 89.2, 3.4)
#' )
#'
#' # Run system suitability
#' sst <- measure_sec_suitability(
#'   peaks = peaks,
#'   reference_peaks = c("dimer", "monomer"),
#'   injected_mass = 0.200,
#'   column_length = 30
#' )
#'
#' print(sst)
#' }
measure_sec_suitability <- function(
    data = NULL,
    peaks,
    reference_peaks = NULL,
    injected_mass = NULL,
    criteria = NULL,
    column_length = NULL
) {
  # Default criteria (typical biopharmaceutical)
  if (is.null(criteria)) {
    criteria <- list(
      resolution_min = 1.5,
      plate_count_min = 5000,
      tailing_min = 0.8,
      tailing_max = 1.5,
      recovery_min = 95,
      recovery_max = 105,
      retention_rsd_max = 1.0,
      area_rsd_max = 2.0
    )
  }

  # Initialize results
  results <- list()
  all_passed <- TRUE

  # Validate peaks input
  if (!is.data.frame(peaks)) {
    cli::cli_abort("{.arg peaks} must be a data frame.")
  }

  required_cols <- c("retention", "width")
  missing_cols <- setdiff(required_cols, names(peaks))
  if (length(missing_cols) > 0) {
    cli::cli_abort(
      "Missing required columns in {.arg peaks}: {.val {missing_cols}}."
    )
  }

  # Resolution (if reference peaks specified)
  if (!is.null(reference_peaks) && length(reference_peaks) >= 2) {
    if ("name" %in% names(peaks)) {
      peak1_idx <- which(peaks$name == reference_peaks[1])
      peak2_idx <- which(peaks$name == reference_peaks[2])
    } else {
      peak1_idx <- 1
      peak2_idx <- 2
    }

    if (length(peak1_idx) > 0 && length(peak2_idx) > 0) {
      Rs <- measure_sec_resolution(
        retention_1 = peaks$retention[peak1_idx[1]],
        retention_2 = peaks$retention[peak2_idx[1]],
        width_1 = peaks$width[peak1_idx[1]],
        width_2 = peaks$width[peak2_idx[1]]
      )

      Rs_pass <- Rs >= criteria$resolution_min
      all_passed <- all_passed && Rs_pass

      results$resolution <- list(
        value = Rs,
        criterion = paste0(">= ", criteria$resolution_min),
        passed = Rs_pass,
        peaks = paste(reference_peaks[1:2], collapse = " / ")
      )
    }
  }

  # Plate count (use largest peak or first peak)
  if ("area" %in% names(peaks)) {
    main_peak_idx <- which.max(peaks$area)
  } else {
    main_peak_idx <- 1
  }

  N <- measure_sec_plate_count(
    retention = peaks$retention[main_peak_idx],
    width = peaks$width[main_peak_idx]
  )

  N_pass <- N >= criteria$plate_count_min
  all_passed <- all_passed && N_pass

  results$plate_count <- list(
    value = round(N),
    criterion = paste0(">= ", criteria$plate_count_min),
    passed = N_pass
  )

  # Plates per meter (if column length provided)
  if (!is.null(column_length) && column_length > 0) {
    plates_per_m <- N / (column_length / 100)
    results$plates_per_meter <- list(
      value = round(plates_per_m),
      criterion = "informational",
      passed = NA
    )
  }

  # Tailing factor (if leading/tailing widths available, otherwise skip)
  if ("leading" %in% names(peaks) && "tailing" %in% names(peaks)) {
    Tf <- measure_sec_asymmetry(
      leading = peaks$leading[main_peak_idx],
      tailing = peaks$tailing[main_peak_idx]
    )

    Tf_pass <- Tf >= criteria$tailing_min && Tf <= criteria$tailing_max
    all_passed <- all_passed && Tf_pass

    results$tailing_factor <- list(
      value = round(Tf, 2),
      criterion = paste0(criteria$tailing_min, "-", criteria$tailing_max),
      passed = Tf_pass
    )
  }

  # Mass recovery (if injected mass and areas available)
  if (!is.null(injected_mass) && "area" %in% names(peaks)) {
    total_area <- sum(peaks$area, na.rm = TRUE)
    # Assume area is proportional to mass (needs calibration in practice)
    # This is a placeholder - real implementation needs response factor
    recovery <- measure_sec_recovery(
      detected_mass = total_area / 100 * injected_mass,  # Simplified
      injected_mass = injected_mass
    )

    recovery_pass <- recovery >= criteria$recovery_min &&
      recovery <= criteria$recovery_max
    all_passed <- all_passed && recovery_pass

    results$recovery <- list(
      value = round(recovery, 1),
      criterion = paste0(criteria$recovery_min, "-", criteria$recovery_max, "%"),
      passed = recovery_pass
    )
  }

  # Retention time RSD (if multiple replicates in peaks)
  if ("replicate" %in% names(peaks) && "name" %in% names(peaks)) {
    # Calculate RSD per peak
    rsd_data <- peaks |>
      dplyr::group_by(.data$name) |>
      dplyr::summarize(
        mean_rt = mean(.data$retention, na.rm = TRUE),
        sd_rt = stats::sd(.data$retention, na.rm = TRUE),
        rsd_rt = 100 * .data$sd_rt / .data$mean_rt,
        .groups = "drop"
      )

    max_rsd <- max(rsd_data$rsd_rt, na.rm = TRUE)

    if (!is.na(max_rsd) && is.finite(max_rsd)) {
      rsd_pass <- max_rsd <= criteria$retention_rsd_max
      all_passed <- all_passed && rsd_pass

      results$retention_rsd <- list(
        value = round(max_rsd, 2),
        criterion = paste0("<= ", criteria$retention_rsd_max, "%"),
        passed = rsd_pass
      )
    }
  }

  # Build summary
  summary_lines <- character()
  for (metric in names(results)) {
    r <- results[[metric]]
    status <- if (is.na(r$passed)) "INFO" else if (r$passed) "PASS" else "FAIL"
    summary_lines <- c(
      summary_lines,
      sprintf("%-20s: %s (%s) [%s]",
              gsub("_", " ", metric),
              format(r$value, nsmall = 1),
              r$criterion,
              status)
    )
  }

  # Create output object
  output <- list(
    results = results,
    passed = all_passed,
    summary = paste(summary_lines, collapse = "\n"),
    criteria = criteria,
    n_peaks = nrow(peaks)
  )

  class(output) <- c("sec_suitability", "list")

  output
}


#' @export
print.sec_suitability <- function(x, ...) {
  cat("SEC System Suitability Test\n")
  cat(strrep("=", 50), "\n\n")

  cat("Overall Status: ")
  if (x$passed) {
    cat("PASSED\n\n")
  } else {
    cat("FAILED\n\n")
  }

  cat("Results:\n")
  cat(strrep("-", 50), "\n")
  cat(x$summary, "\n")
  cat(strrep("-", 50), "\n")

  invisible(x)
}


#' @export
summary.sec_suitability <- function(object, ...) {
  cat("SEC System Suitability Summary\n\n")

  n_metrics <- length(object$results)
  n_passed <- sum(vapply(object$results, function(r) {
    if (is.na(r$passed)) return(FALSE)
    r$passed
  }, logical(1)))
  n_failed <- sum(vapply(object$results, function(r) {
    if (is.na(r$passed)) return(FALSE)
    !r$passed
  }, logical(1)))
  n_info <- sum(vapply(object$results, function(r) is.na(r$passed), logical(1)))

  cat("Metrics evaluated:", n_metrics, "\n")
  cat("  Passed:", n_passed, "\n")
  cat("  Failed:", n_failed, "\n")
  if (n_info > 0) cat("  Informational:", n_info, "\n")

  cat("\nOverall:", if (object$passed) "PASSED" else "FAILED", "\n")

  invisible(object)
}
