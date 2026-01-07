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

  # Remove NAs for processing
  valid_idx <- !is.na(value)
  loc_valid <- location[valid_idx]
  val_valid <- value[valid_idx]
  n_valid <- length(val_valid)

  if (n_valid < 10) {
    return(measure:::new_peaks_tbl())
  }

  # Step 1: LOESS smoothing
  loess_fit <- tryCatch(
    stats::loess(
      val_valid ~ loc_valid,
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
    smoothed <- stats::predict(loess_fit, newdata = loc_valid)
    # Handle any NAs from prediction
    smoothed[is.na(smoothed)] <- val_valid[is.na(smoothed)]
  } else {
    smoothed <- val_valid
  }

  # Ensure data is sorted by location
  ord <- order(loc_valid)
  loc_sorted <- loc_valid[ord]
  val_sorted <- smoothed[ord]

  # Step 2: Iterative Soft Thresholding (IST) with nonlinear spacing
  # Estimate baseline variance from the lowest portion of the signal
  baseline_region <- val_sorted < stats::quantile(val_sorted, 0.1)
  if (sum(baseline_region) < 3) {
    baseline_region <- val_sorted < stats::quantile(val_sorted, 0.25)
  }
  baseline_var <- if (sum(baseline_region) > 2) {
    stats::var(val_sorted[baseline_region])
  } else {
    stats::var(val_sorted) / 100 # fallback
  }

  # Generate nonlinear threshold spacing
  # Concentrates thresholds near the baseline for better sensitivity
  nonlinear_spacing <- function(n, a) {
    nlseq <- 1 / exp(seq(0, a, length.out = n))
    nlseq / max(nlseq)
  }

  min_val <- min(val_sorted, na.rm = TRUE)
  max_val <- max(val_sorted, na.rm = TRUE)

  if (abs(max_val - min_val) < .Machine$double.eps * 100) {
    # Signal is essentially flat
    return(measure:::new_peaks_tbl())
  }

  thresholds <- min_val +
    (max_val - min_val) *
      nonlinear_spacing(ist_points, ist_nonlinearity)

  # Count points at each threshold level
  point_counts <- vapply(
    thresholds,
    function(thresh) {
      as.integer(sum(abs(val_sorted - thresh) <= sqrt(baseline_var)))
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
      cpt_idx <- min(changepoint::cpts(cpt_result))
      if (length(cpt_idx) == 0 || cpt_idx < 1 || cpt_idx > length(thresholds)) {
        # Fallback: use 10th percentile
        stats::quantile(val_sorted, 0.1)
      } else {
        thresholds[cpt_idx]
      }
    },
    error = function(e) {
      cli::cli_warn(
        "Changepoint detection failed: {e$message}. Using percentile threshold."
      )
      stats::quantile(val_sorted, 0.1)
    }
  )

  # Step 4: Extract points above threshold
  above_thresh <- val_sorted >= threshold
  if (sum(above_thresh) < 2) {
    return(measure:::new_peaks_tbl())
  }

  # Step 5: Peak detection by walking through above-threshold regions
  peaks <- .find_peak_boundaries(
    loc_sorted,
    val_sorted,
    above_thresh,
    min_height
  )

  peaks
}

#' Find peak boundaries from thresholded signal
#'
#' @param location Sorted location vector.
#' @param value Sorted (smoothed) value vector.
#' @param above_thresh Logical vector indicating points above threshold.
#' @param min_height Minimum peak height.
#' @return A peaks_tbl object.
#' @noRd
.find_peak_boundaries <- function(location, value, above_thresh, min_height) {
  # Get indices of above-threshold points
  idx_above <- which(above_thresh)
  n_above <- length(idx_above)

  if (n_above < 2) {
    return(measure:::new_peaks_tbl())
  }

  # Initialize peak storage
  peaks <- list()
  peak_count <- 0L

  # Walk through above-threshold points to find peaks
  i <- 1
  while (i < n_above) {
    # Start of potential peak region
    start_idx <- idx_above[i]

    # Find apex (highest point in this connected region)
    j <- i
    apex_idx <- start_idx
    apex_val <- value[start_idx]

    # Trace upward
    while (j < n_above && value[idx_above[j + 1]] >= value[idx_above[j]]) {
      j <- j + 1
      if (value[idx_above[j]] > apex_val) {
        apex_idx <- idx_above[j]
        apex_val <- value[idx_above[j]]
      }
    }

    # Apex found at j, now trace downward
    apex_j <- j

    while (j < n_above && value[idx_above[j + 1]] <= value[idx_above[j]]) {
      j <- j + 1
    }

    # End of peak region
    end_idx <- idx_above[j]

    # Calculate height (apex - mean of start and end values)
    base_height <- (value[start_idx] + value[end_idx]) / 2
    peak_height <- apex_val - base_height

    # Only keep if height exceeds threshold
    if (peak_height >= min_height) {
      peak_count <- peak_count + 1
      peaks[[peak_count]] <- list(
        peak_id = peak_count,
        location = location[apex_idx],
        height = peak_height,
        left_base = location[start_idx],
        right_base = location[end_idx],
        area = NA_real_
      )
    }

    # Move to next potential peak
    i <- j + 1
  }

  # Convert to peaks_tbl
  if (peak_count == 0) {
    return(measure:::new_peaks_tbl())
  }

  measure:::new_peaks_tbl(
    peak_id = vapply(peaks, function(p) as.integer(p$peak_id), integer(1)),
    location = vapply(peaks, `[[`, double(1), "location"),
    height = vapply(peaks, `[[`, double(1), "height"),
    left_base = vapply(peaks, `[[`, double(1), "left_base"),
    right_base = vapply(peaks, `[[`, double(1), "right_base"),
    area = vapply(peaks, `[[`, double(1), "area")
  )
}
