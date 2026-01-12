# ==============================================================================
# step_sec_peaks_detect
#
# SEC/GPC optimized peak detection with finderskeepers algorithm
# ==============================================================================

#' SEC/GPC Peak Detection
#'
#' `step_sec_peaks_detect()` creates a *specification* of a recipe step that
#' detects peaks in SEC/GPC chromatography data. The default algorithm is
#' `finderskeepers`, which uses LOESS smoothing with Iterative Soft Thresholding
#' (IST) and changepoint analysis for automatic threshold detection.
#'
#' @param recipe A recipe object. The step will be added to the sequence of
#'   operations for this recipe.
#' @param measures An optional character vector of measure column names to
#'   process. If `NULL` (the default), all measure columns (columns with class
#'   `measure_list`) will be processed.
#' @param algorithm Peak detection algorithm. Currently only `"finderskeepers"`
#'   (the default) is supported. This SEC-optimized algorithm uses LOESS
#'   smoothing with Iterative Soft Thresholding and changepoint detection.
#'   For other algorithms, use [measure::step_measure_peaks_detect()] directly.
#' @param min_height Minimum peak height. For `"finderskeepers"`, this is the
#'   minimum height above baseline. For other algorithms with `snr_threshold = TRUE`,
#'   this is interpreted as a signal-to-noise ratio. Default is `1`.
#' @param min_distance Minimum distance between peaks in location units (e.g., mL).
#'   Only used for non-finderskeepers algorithms. Default is `0`.
#' @param loess_span LOESS smoothing span for `finderskeepers`. A value between 0
#'   and 1 controlling the smoothness of the fit. Default is `0.01` (minimal
#'   smoothing to preserve peak shapes).
#' @param ist_points Number of threshold levels for Iterative Soft Thresholding.
#'   Default is `50`. Higher values give finer threshold resolution.
#' @param ist_nonlinearity Nonlinearity parameter for IST threshold spacing.
#'   Default is `5`. Higher values concentrate thresholds near the baseline.
#' @param snr_threshold Logical. For non-finderskeepers algorithms, if `TRUE`,
#'   `min_height` is interpreted as a signal-to-noise ratio. Default is `FALSE`.
#' @param role Not used by this step since no new variables are created.
#' @param trained A logical to indicate if the quantities for preprocessing
#'   have been estimated.
#' @param skip A logical. Should the step be skipped when the recipe is baked?
#' @param id A character string that is unique to this step to identify it.
#'
#' @return An updated version of `recipe` with the new step added to the
#'
#'   sequence of any existing operations. A new `.peaks` column will be added
#'   containing detected peaks for each sample.
#'
#' @details
#' The `finderskeepers` algorithm is specifically designed for SEC/GPC data:
#'
#' 1
#' . **LOESS smoothing**: Applies local polynomial regression to reduce noise
#'    while preserving peak shapes (controlled by `loess_span`)
#' 2. **Iterative Soft Thresholding (IST)**: Creates a series of thresholds with
#'    nonlinear spacing to detect changes in peak structure
#' 3. **Changepoint detection**: Uses `changepoint::cpt.mean()` to automatically
#'    determine the optimal threshold for peak/baseline separation
#' 4. **Peak boundary detection**: Identifies peak start, apex, and end points
#'
#' This approach is robust to baseline drift and varying peak heights, making it
#' well-suited for polymer SEC chromatograms.
#'
#' **Peak properties stored:**
#' - `peak_id`: Integer identifier
#' - `location`: X-axis position of peak apex (elution volume)
#' - `height`: Peak height above baseline
#' - `left_base`, `right_base`: X-axis positions of peak boundaries
#' - `area`: Initially NA; use [measure::step_measure_peaks_integrate()] to calculate
#'
#' # Tidying
#'
#' When you [`tidy()`][recipes::tidy.recipe()] this step, a tibble with columns
#' `terms`, `algorithm`, `min_height`, and `id` is returned.
#'
#' @seealso [measure::step_measure_peaks_detect()] for general peak detection,
#'   [measure::step_measure_peaks_integrate()] to calculate peak areas.
#'
#' @family sec-chromatography
#' @export
#'
#' @examples
#' \dontrun{
#' library(recipes)
#' library(measure)
#'
#' # SEC peak detection with finderskeepers (default)
#' rec <- recipe(~., data = sec_triple_detect) |>
#'   step_measure_input_long(ri_signal, location = vars(elution_time)) |>
#'   step_sec_peaks_detect() |>
#'   prep()
#'
#' # Adjust sensitivity with min_height
#' rec <- recipe(~., data = sec_triple_detect) |>
#'   step_measure_input_long(ri_signal, location = vars(elution_time)) |>
#'   step_sec_peaks_detect(min_height = 5) |>
#'   prep()
#'
#' # Adjust LOESS smoothing for noisy data
#' rec <- recipe(~., data = sec_triple_detect) |>
#'   step_measure_input_long(ri_signal, location = vars(elution_time)) |>
#'   step_sec_peaks_detect(loess_span = 0.05) |>
#'   prep()
#' }
step_sec_peaks_detect <- function(
  recipe,
  measures = NULL,
  algorithm = "finderskeepers",
  min_height = 1,
  min_distance = 0,
  loess_span = 0.01,
  ist_points = 50L,
  ist_nonlinearity = 5,
  snr_threshold = FALSE,
  role = NA,
  trained = FALSE,
  skip = FALSE,
  id = recipes::rand_id("sec_peaks_detect")
) {
  # Validate algorithm - currently only finderskeepers is supported
  if (algorithm != "finderskeepers") {
    cli::cli_abort(
      c(
        "Currently only {.val finderskeepers} algorithm is supported.",
        "i" = "For other algorithms, use {.fn measure::step_measure_peaks_detect}."
      )
    )
  }

  # Validate numeric parameters
  if (!is.numeric(min_height) || length(min_height) != 1 || min_height < 0) {
    cli::cli_abort("{.arg min_height} must be a non-negative number.")
  }
  if (
    !is.numeric(min_distance) || length(min_distance) != 1 || min_distance < 0
  ) {
    cli::cli_abort("{.arg min_distance} must be a non-negative number.")
  }
  if (
    !is.numeric(loess_span) ||
      length(loess_span) != 1 ||
      loess_span <= 0 ||
      loess_span > 1
  ) {
    cli::cli_abort("{.arg loess_span} must be a number between 0 and 1.")
  }
  if (!is.numeric(ist_points) || length(ist_points) != 1 || ist_points < 10) {
    cli::cli_abort("{.arg ist_points} must be at least 10.")
  }
  if (
    !is.numeric(ist_nonlinearity) ||
      length(ist_nonlinearity) != 1 ||
      ist_nonlinearity <= 0
  ) {
    cli::cli_abort("{.arg ist_nonlinearity} must be a positive number.")
  }

  recipes::add_step(
    recipe,
    step_sec_peaks_detect_new(
      measures = measures,
      algorithm = algorithm,
      min_height = min_height,
      min_distance = min_distance,
      loess_span = loess_span,
      ist_points = as.integer(ist_points),
      ist_nonlinearity = ist_nonlinearity,
      snr_threshold = snr_threshold,
      role = role,
      trained = trained,
      skip = skip,
      id = id
    )
  )
}

step_sec_peaks_detect_new <- function(
  measures,
  algorithm,
  min_height,
  min_distance,
  loess_span,
  ist_points,
  ist_nonlinearity,
  snr_threshold,
  role,
  trained,
  skip,
  id
) {
  recipes::step(
    subclass = "sec_peaks_detect",
    measures = measures,
    algorithm = algorithm,
    min_height = min_height,
    min_distance = min_distance,
    loess_span = loess_span,
    ist_points = ist_points,
    ist_nonlinearity = ist_nonlinearity,
    snr_threshold = snr_threshold,
    role = role,
    trained = trained,
    skip = skip,
    id = id
  )
}

#' @export
prep.step_sec_peaks_detect <- function(x, training, info = NULL, ...) {
  check_for_measure(training)

  # Resolve which columns to process
  if (is.null(x$measures)) {
    measure_cols <- find_measure_cols(training)
  } else {
    measure_cols <- x$measures
  }

  step_sec_peaks_detect_new(
    measures = measure_cols,
    algorithm = x$algorithm,
    min_height = x$min_height,
    min_distance = x$min_distance,
    loess_span = x$loess_span,
    ist_points = x$ist_points,
    ist_nonlinearity = x$ist_nonlinearity,
    snr_threshold = x$snr_threshold,
    role = x$role,
    trained = TRUE,
    skip = x$skip,
    id = x$id
  )
}

#' @export
bake.step_sec_peaks_detect <- function(object, new_data, ...) {
  algorithm <- object$algorithm
  min_height <- object$min_height
  min_distance <- object$min_distance
  loess_span <- object$loess_span
  ist_points <- object$ist_points
  ist_nonlinearity <- object$ist_nonlinearity
  snr_threshold <- object$snr_threshold

  # Process each measure column
  for (col in object$measures) {
    peaks_col_name <- sub("^\\.?measures?", ".peaks", col)
    if (peaks_col_name == col) {
      peaks_col_name <- ".peaks"
    }

    peaks_list <- purrr::map(new_data[[col]], function(m) {
      loc <- m$location
      val <- m$value

      # Use SEC-specific finderskeepers algorithm
      .detect_peaks_finderskeepers(
        location = loc,
        value = val,
        min_height = min_height,
        loess_span = loess_span,
        ist_points = ist_points,
        ist_nonlinearity = ist_nonlinearity
      )
    })

    # Use internal constructor from measure (should be exported in future)
    new_data[[peaks_col_name]] <- measure:::new_peaks_list(peaks_list)
  }

  tibble::as_tibble(new_data)
}

#' @export
print.step_sec_peaks_detect <- function(
  x,
  width = max(20, options()$width - 30),
  ...
) {
  title <- glue::glue("SEC peak detection ({x$algorithm}) on ")

  if (x$trained) {
    cat(title, "<internal measurements>", sep = "")
  } else {
    cat(title)
  }
  cat("\n")

  invisible(x)
}

#' @rdname tidy.step_sec
#' @export
#' @keywords internal
tidy.step_sec_peaks_detect <- function(x, ...) {
  if (recipes::is_trained(x)) {
    terms_val <- x$measures
  } else {
    terms_val <- "<all measure columns>"
  }
  tibble::tibble(
    terms = terms_val,
    algorithm = x$algorithm,
    min_height = x$min_height,
    min_distance = x$min_distance,
    loess_span = x$loess_span,
    ist_points = x$ist_points,
    id = x$id
  )
}

#' @rdname required_pkgs.step_sec
#' @export
#' @keywords internal
required_pkgs.step_sec_peaks_detect <- function(x, ...) {
  pkgs <- c("measure.sec", "measure")
  if (!is.null(x$algorithm) && x$algorithm == "finderskeepers") {
    pkgs <- c(pkgs, "changepoint")
  }
  pkgs
}

# ==============================================================================
# Finderskeepers algorithm implementation
# ==============================================================================

#' SEC Peak Detection using Finderskeepers Algorithm
#'
#' Detects peaks in SEC/GPC chromatography data using LOESS smoothing,
#' Iterative Soft Thresholding (IST), and changepoint detection.
#'
#' @param location Numeric vector of x-axis values (elution volume/time).
#' @param value Numeric vector of y-axis values (signal intensity).
#' @param min_height Minimum peak height above baseline.
#' @param loess_span LOESS smoothing span (0 to 1).
#' @param ist_points Number of IST threshold levels.
#' @param ist_nonlinearity Nonlinearity for threshold spacing.
#'
#' @return A `peaks_tbl` object with detected peaks.
#' @noRd
.detect_peaks_finderskeepers <- function(
  location,
  value,
  min_height = 1,
  loess_span = 0.01,
  ist_points = 50L,
  ist_nonlinearity = 5
) {
  n <- length(value)

  # Handle edge cases
  if (n < 10) {
    return(measure:::new_peaks_tbl())
  }

  if (all(is.na(value))) {
    return(measure:::new_peaks_tbl())
  }

  # Remove NAs for processing (filter both location and value)
  valid_idx <- !is.na(value) & !is.na(location)
  loc_valid <- location[valid_idx]
  val_valid <- value[valid_idx]
  n_valid <- length(val_valid)

  if (n_valid < 10) {
    return(measure:::new_peaks_tbl())
  }

  # Ensure data is sorted by location first
  ord <- order(loc_valid)
  loc_sorted <- loc_valid[ord]
  val_sorted <- val_valid[ord]

  # Step 1: LOESS smoothing
  loess_fit <- tryCatch(
    stats::loess(
      val_sorted ~ loc_sorted,
      span = loess_span,
      degree = 2,
      family = "gaussian"
    ),
    error = function(e) {
      cli::cli_warn("LOESS fitting failed: {e$message}. Using raw values.")
      NULL
    }
  )

  if (!is.null(loess_fit)) {
    smoothed <- stats::predict(loess_fit, newdata = loc_sorted)
    # Handle any NAs from prediction
    smoothed[is.na(smoothed)] <- val_sorted[is.na(smoothed)]
    val_sorted <- smoothed
  }

  # Step 2: Iterative Soft Thresholding (IST) with nonlinear spacing
  # Use points with location < 0.05 as baseline region. When data starts at
  # higher values (e.g., SEC elution volumes of 5-25 mL), this returns NA,
  # resulting in all-zero IST counts and a percentile-based threshold fallback.
  early_region <- loc_sorted < 0.05

  # Baseline variance from early chromatogram region
  baseline_var <- if (sum(early_region) > 0) {
    stats::var(val_sorted[early_region], na.rm = TRUE)
  } else {
    cli::cli_inform(c(
      "i" = "No data points below x=0.05 for baseline estimation.",
      "i" = "Using percentile-based threshold for peak detection."
    ))
    NA_real_
  }

  # Generate nonlinear threshold spacing
  nonlinear_spacing <- function(n, a) {
    nlseq <- 1 / exp(seq(0, a, length.out = n))
    nlseq / max(nlseq)
  }

  # Use first value as reference (matches original finderskeepers)
  first_val <- val_sorted[1]
  max_val <- max(val_sorted, na.rm = TRUE)

  if (abs(max_val - first_val) < .Machine$double.eps * 100) {
    # Signal is essentially flat
    return(measure:::new_peaks_tbl())
  }

  thresholds <- first_val +
    (max_val - first_val) *
      nonlinear_spacing(ist_points, ist_nonlinearity)

  # Count points at each threshold level
  # When baseline_var is NA (no early region data), return 0 to trigger fallback threshold
  point_counts <- vapply(
    thresholds,
    function(thresh) {
      if (is.na(baseline_var)) {
        0L
      } else {
        as.integer(sum(abs(val_sorted - thresh) <= baseline_var, na.rm = TRUE))
      }
    },
    integer(1)
  )

  # Step 3: Changepoint detection to find optimal threshold
  rlang::check_installed(
    "changepoint",
    reason = "for finderskeepers peak detection"
  )

  threshold <- tryCatch(
    {
      cpt_result <- changepoint::cpt.mean(point_counts)
      cpts <- changepoint::cpts(cpt_result)

      if (length(cpts) == 0) {
        # No changepoint found (e.g., all counts are 0)
        # Use 10th percentile as fallback
        stats::quantile(val_sorted, 0.1, na.rm = TRUE)
      } else {
        cpt_idx <- min(cpts)
        if (cpt_idx < 1 || cpt_idx > length(thresholds)) {
          stats::quantile(val_sorted, 0.1, na.rm = TRUE)
        } else {
          thresholds[cpt_idx]
        }
      }
    },
    error = function(e) {
      cli::cli_warn(
        "Changepoint detection failed: {e$message}. Using percentile threshold."
      )
      stats::quantile(val_sorted, 0.1, na.rm = TRUE)
    }
  )

  # Step 4: Extract points above threshold
  above_thresh <- val_sorted >= threshold
  if (sum(above_thresh) < 2) {
    return(measure:::new_peaks_tbl())
  }

  # Step 5: Peak boundary detection using walk algorithm
  peaks <- .find_peak_boundaries(
    loc_sorted,
    val_sorted,
    above_thresh,
    min_height
  )

  peaks
}

#' Find peak boundaries using walk algorithm
#'
#' Walks through above-threshold points to identify peak start, apex, and end
#' positions using an upward-then-downward traversal pattern.
#'
#' @param location Sorted location vector (x-axis, e.g., elution volume).
#' @param value Sorted (smoothed) value vector (y-axis, e.g., signal).
#' @param above_thresh Logical vector indicating points above threshold.
#' @param min_height Minimum peak height to keep.
#' @return A peaks_tbl object.
#' @noRd
.find_peak_boundaries <- function(
  location,
  value,
  above_thresh,
  min_height
) {
  # Get the culled data (points above threshold)
  idx_above <- which(above_thresh)
  n_above <- length(idx_above)

  if (n_above < 2) {
    return(measure:::new_peaks_tbl())
  }

  # Extract x and y for above-threshold points
  xvar <- location[idx_above]
  yvar <- value[idx_above]

  # Initialize peak storage - pre-allocate for up to 50 peaks
  # (sufficient for typical SEC chromatograms; warns if exceeded)
  max_peaks <- 50L
  peaks_start <- rep(0, max_peaks)
  peaks_apex <- rep(0, max_peaks)
  peaks_end <- rep(0, max_peaks)

  # Walk through the data using original algorithm
  counter <- 1L
  starter <- 1L
  chaser <- starter
  scouter <- chaser + 1L

  while (scouter <= n_above) {
    # Trace upward (while signal is increasing)
    while (scouter <= n_above && yvar[chaser] < yvar[scouter]) {
      chaser <- chaser + 1L
      scouter <- scouter + 1L
    }

    # Record start and apex
    if (
      starter <= n_above && (scouter - 1L) >= 1L && (scouter - 1L) <= n_above
    ) {
      peaks_start[counter] <- xvar[starter]
      peaks_apex[counter] <- xvar[scouter - 1L]
    }

    chaser <- chaser + 1L
    scouter <- scouter + 1L

    # Trace downward (while signal is decreasing)
    while (scouter <= n_above && yvar[chaser] > yvar[scouter]) {
      chaser <- chaser + 1L
      scouter <- scouter + 1L
    }

    # Record end
    if ((scouter - 2L) >= 1L && (scouter - 2L) <= n_above) {
      peaks_end[counter] <- xvar[scouter - 2L]
    }

    # Move to next peak
    starter <- chaser + 1L
    chaser <- starter
    scouter <- chaser + 1L
    counter <- counter + 1L

    # Safety check - warn if we hit the limit
    if (counter > max_peaks) {
      cli::cli_warn(c(
        "!" = "Peak detection reached maximum limit of {max_peaks} peaks.",
        "i" = "Additional peaks beyond this limit are not reported.",
        "i" = "Consider increasing {.arg min_height} to reduce peak count."
      ))
      break
    }
  }

  # Filter to valid peaks (start > 0) and calculate heights
  valid_peaks <- peaks_start > 0
  n_valid <- sum(valid_peaks)

  if (n_valid == 0) {
    return(measure:::new_peaks_tbl())
  }

  # Extract valid peaks
  valid_start <- peaks_start[valid_peaks]
  valid_apex <- peaks_apex[valid_peaks]
  valid_end <- peaks_end[valid_peaks]

  # Calculate heights: apex_value - start_value
  heights <- vapply(
    seq_len(n_valid),
    function(i) {
      apex_val <- value[which.min(abs(location - valid_apex[i]))]
      start_val <- value[which.min(abs(location - valid_start[i]))]
      apex_val - start_val
    },
    double(1)
  )

  # Filter by min_height
  keep <- heights > min_height
  if (!any(keep)) {
    return(measure:::new_peaks_tbl())
  }

  # Build final peaks table
  final_start <- valid_start[keep]
  final_apex <- valid_apex[keep]
  final_end <- valid_end[keep]
  final_heights <- heights[keep]
  n_final <- length(final_start)

  measure:::new_peaks_tbl(
    peak_id = seq_len(n_final),
    location = final_apex,
    height = final_heights,
    left_base = final_start,
    right_base = final_end,
    area = rep(NA_real_, n_final)
  )
}
