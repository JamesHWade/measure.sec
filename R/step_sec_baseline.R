# ==============================================================================
# step_sec_baseline
#
# SEC/GPC optimized baseline correction
# ==============================================================================

#' SEC/GPC Baseline Correction
#'
#' `step_sec_baseline()` creates a *specification* of a recipe step
#' that applies baseline correction optimized for Gel Permeation Chromatography
#' (GPC) or Size Exclusion Chromatography (SEC) data. This method estimates the
#' baseline by interpolating between baseline regions at the start and end of
#' the chromatogram.
#'
#' @param recipe A recipe object. The step will be added to the sequence of
#'   operations for this recipe.
#' @param measures An optional character vector of measure column names to
#'   process. If `NULL` (the default), all measure columns (columns with class
#'   `measure_list`) will be processed.
#' @param left_frac Fraction of points from the beginning to use as the left
#'   baseline region. Default is `0.05` (first 5% of data points). Not used for
#'   `method = "rf"`.
#' @param right_frac Fraction of points from the end to use as the right
#'   baseline region. Default is `0.05` (last 5% of data points). Not used for
#'   `method = "rf"`.
#' @param method Method for baseline estimation. One of:
#'   - `"linear"` (default): Linear interpolation between left and right means
#'   - `"median"`: Uses median of baseline regions (more robust to outliers)
#'   - `"spline"`: Smooth spline through baseline regions
#'   - `"rf"`: Robust local regression (IDPmisc::rfbaseline) for curved or
#'     drifting baselines. Does not use left_frac/right_frac.
#' @param rf_span Span parameter for RF baseline (only used when `method = "rf"`).
#'   A numeric value between 0 and 1 specifying the fraction of points used for
#'   local regression. Default is `2/3`. Higher values produce smoother baselines.
#' @param rf_maxit Maximum iterations for RF baseline (only used when `method = "rf"`).
#'   A numeric vector of length 2 specifying max iterations for the two stages
#'   of the robust fitting. Default is `c(20, 20)`.
#' @param role Not used by this step since no new variables are created.
#' @param trained A logical to indicate if the quantities for preprocessing
#'   have been estimated.
#' @param skip A logical. Should the step be skipped when the recipe is baked?
#' @param id A character string that is unique to this step to identify it.
#'
#' @return An updated version of `recipe` with the new step added to the
#'   sequence of any existing operations.
#'
#' @details
#' When the beginning and end of a chromatogram contain baseline-only
#' observations, this step estimates a baseline by:
#'
#' 1. Identifying baseline regions at the start and end of the chromatogram
#' 2. Computing a representative baseline value for each region (mean or median)
#' 3. Interpolating between these values to estimate the full baseline
#' 4. Subtracting the estimated baseline from the signal
#'
#' The `left_frac` and `right_frac` parameters control how much of the
#' chromatogram is considered "baseline". Choose values that:
#' - Include only the flat, signal-free regions
#' - Exclude any polymer peaks or system peaks
#' - Are large enough to average out noise
#'
#' Inspect the selected edge regions before subtraction. A peak that extends
#' into either region can contribute to the estimated baseline and lose part
#' of its response during subtraction.
#'
#' For `"linear"`, `"median"`, and `"spline"`, each selected edge region must
#' contain at least one finite signal value. Otherwise, correction stops with
#' an error identifying the unsupported edge, column, and curve. Only finite
#' values contribute to edge summaries. Linear interpolation pairs each edge
#' mean with the mean of its supported point indices; the median method pairs
#' each edge median with the median of those indices. Spline fitting uses finite
#' edge values at their original indices and requires four distinct supported
#' indices; fewer indices or a failed fit produce a warning and use the linear
#' method. Correction retains all point positions and preserves `NA`, `NaN`,
#' `Inf`, and `-Inf` signal values. These support requirements do not apply to
#' `"rf"`. For every method, curves with fewer than 10 points and curves with
#' only missing values (`NA` or `NaN`) return unchanged with a warning.
#'
#' Baseline subtraction retains a recorded `g/mL` or `mg/mL` unit on the
#' current input column. Untagged inputs remain untagged, and invalid recorded
#' units are rejected. The unit tag describes the values' units; it does not
#' qualify the selected baseline regions or establish recovery.
#'
#' **No selectors should be supplied to this step function**. The data should be
#' in the internal format produced by [measure::step_measure_input_wide()] or
#' [measure::step_measure_input_long()].
#'
#' # Tidying
#'
#' When you [`tidy()`][recipes::tidy.recipe()] this step, a tibble with columns
#' `terms`, `left_frac`, `right_frac`, `method`, and `id` is returned.
#'
#' @seealso [measure::step_measure_baseline_als()] for general-purpose baseline
#'   correction, [measure::step_measure_detrend()] for simple trend removal.
#'
#' @family sec-chromatography
#' @export
#'
#' @examples
#' \dontrun{
#' library(recipes)
#' library(measure)
#'
#' # SEC baseline correction with default settings
#' rec <- recipe(~., data = sec_triple_detect) |>
#'   step_measure_input_wide(starts_with("signal_")) |>
#'   step_sec_baseline() |>
#'   prep()
#'
#' # Using median method for robustness to outliers
#' rec <- recipe(~., data = sec_triple_detect) |>
#'   step_measure_input_wide(starts_with("signal_")) |>
#'   step_sec_baseline(left_frac = 0.1, right_frac = 0.1, method = "median") |>
#'   prep()
#'
#' # Using RF baseline for complex baseline shapes
#' rec <- recipe(~., data = sec_triple_detect) |>
#'   step_measure_input_wide(starts_with("signal_")) |>
#'   step_sec_baseline(method = "rf", rf_span = 0.5) |>
#'   prep()
#' }
step_sec_baseline <- function(
  recipe,
  measures = NULL,
  left_frac = 0.05,
  right_frac = 0.05,
  method = "linear",
  rf_span = 2 / 3,
  rf_maxit = c(20, 20),
  role = NA,
  trained = FALSE,
  skip = FALSE,
  id = recipes::rand_id("sec_baseline")
) {
  recipes::add_step(
    recipe,
    step_sec_baseline_new(
      measures = measures,
      left_frac = left_frac,
      right_frac = right_frac,
      method = method,
      rf_span = rf_span,
      rf_maxit = rf_maxit,
      role = role,
      trained = trained,
      skip = skip,
      id = id
    )
  )
}

step_sec_baseline_new <- function(
  measures,
  left_frac,
  right_frac,
  method,
  rf_span,
  rf_maxit,
  role,
  trained,
  skip,
  id
) {
  recipes::step(
    subclass = "sec_baseline",
    measures = measures,
    left_frac = left_frac,
    right_frac = right_frac,
    method = method,
    rf_span = rf_span,
    rf_maxit = rf_maxit,
    role = role,
    trained = trained,
    skip = skip,
    id = id
  )
}

#' @export
prep.step_sec_baseline <- function(x, training, info = NULL, ...) {
  check_for_measure(training)

  # Validate method first
  valid_methods <- c("linear", "median", "spline", "rf")
  if (!x$method %in% valid_methods) {
    cli::cli_abort(
      "{.arg method} must be one of {.or {.val {valid_methods}}}, not {.val {x$method}}."
    )
  }

  # Validate fraction parameters (not used for rf method)
  if (x$method != "rf") {
    if (
      !is.numeric(x$left_frac) ||
        length(x$left_frac) != 1 ||
        x$left_frac <= 0 ||
        x$left_frac >= 0.5
    ) {
      cli::cli_abort(
        "{.arg left_frac} must be a number between 0 and 0.5, not {.val {x$left_frac}}."
      )
    }

    if (
      !is.numeric(x$right_frac) ||
        length(x$right_frac) != 1 ||
        x$right_frac <= 0 ||
        x$right_frac >= 0.5
    ) {
      cli::cli_abort(
        "{.arg right_frac} must be a number between 0 and 0.5, not {.val {x$right_frac}}."
      )
    }

    if (x$left_frac + x$right_frac >= 0.5) {
      cli::cli_abort(
        "{.arg left_frac} + {.arg right_frac} must be less than 0.5."
      )
    }
  }

  # Validate RF-specific parameters
  if (x$method == "rf") {
    if (
      !is.numeric(x$rf_span) ||
        length(x$rf_span) != 1 ||
        x$rf_span <= 0 ||
        x$rf_span > 1
    ) {
      cli::cli_abort(
        "{.arg rf_span} must be a number between 0 and 1, not {.val {x$rf_span}}."
      )
    }

    if (
      !is.numeric(x$rf_maxit) || length(x$rf_maxit) != 2 || any(x$rf_maxit < 1)
    ) {
      cli::cli_abort(
        "{.arg rf_maxit} must be a numeric vector of length 2 with positive values."
      )
    }
  }

  # Resolve which columns to process
  if (is.null(x$measures)) {
    measure_cols <- find_measure_cols(training)
  } else {
    measure_cols <- x$measures
  }

  step_sec_baseline_new(
    measures = measure_cols,
    left_frac = x$left_frac,
    right_frac = x$right_frac,
    method = x$method,
    rf_span = x$rf_span,
    rf_maxit = x$rf_maxit,
    role = x$role,
    trained = TRUE,
    skip = x$skip,
    id = x$id
  )
}

#' @export
bake.step_sec_baseline <- function(object, new_data, ...) {
  for (col in object$measures) {
    concentration_unit <- attr(
      new_data[[col]],
      "sec_concentration_unit",
      exact = TRUE
    )
    .sec_validate_concentration_unit(
      concentration_unit,
      "recorded concentration unit"
    )

    result <- .compute_sec_baseline(
      new_data[[col]],
      left_frac = object$left_frac,
      right_frac = object$right_frac,
      method = object$method,
      rf_span = object$rf_span,
      rf_maxit = object$rf_maxit,
      column = col
    )
    output <- new_measure_list(result)
    if (!is.null(concentration_unit)) {
      attr(output, "sec_concentration_unit") <- concentration_unit
    }
    new_data[[col]] <- output
  }
  tibble::as_tibble(new_data)
}

#' @export
print.step_sec_baseline <- function(
  x,
  width = max(20, options()$width - 30),
  ...
) {
  title <- "SEC baseline correction on "

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
tidy.step_sec_baseline <- function(x, ...) {
  if (recipes::is_trained(x)) {
    terms_val <- x$measures
  } else {
    terms_val <- "<all measure columns>"
  }
  tibble::tibble(
    terms = terms_val,
    left_frac = x$left_frac,
    right_frac = x$right_frac,
    method = x$method,
    rf_span = x$rf_span %||% NA_real_,
    rf_maxit = list(x$rf_maxit %||% NA),
    id = x$id
  )
}

#' @rdname required_pkgs.step_sec
#' @export
#' @keywords internal
required_pkgs.step_sec_baseline <- function(x, ...) {
  pkgs <- c("measure.sec", "measure")
  if (!is.null(x$method) && x$method == "rf") {
    pkgs <- c(pkgs, "IDPmisc")
  }
  pkgs
}

# ------------------------------------------------------------------------------
# Internal computation

#' Compute SEC baseline correction for a measure_list
#'
#' @param dat A measure_list (list of measure_tbl objects).
#' @param left_frac Fraction for left baseline region.
#' @param right_frac Fraction for right baseline region.
#' @param method Baseline estimation method.
#' @param rf_span Span parameter for RF baseline.
#' @param rf_maxit Max iterations for RF baseline.
#' @param column Measure column name for error context.
#' @return A list of measure_tbl objects with baseline-corrected values.
#' @noRd
.compute_sec_baseline <- function(
  dat,
  left_frac,
  right_frac,
  method,
  rf_span = 2 / 3,
  rf_maxit = c(20, 20),
  column = NULL
) {
  purrr::map2(
    dat,
    seq_along(dat),
    function(x, i) {
      .sec_baseline_single(
        x,
        left_frac = left_frac,
        right_frac = right_frac,
        method = method,
        rf_span = rf_span,
        rf_maxit = rf_maxit,
        column = column,
        curve = i
      )
    }
  )
}

#' Apply SEC baseline correction to a single chromatogram
#'
#' @param x A measure_tbl with location and value columns.
#' @param left_frac Fraction for left baseline region.
#' @param right_frac Fraction for right baseline region.
#' @param method Baseline estimation method.
#' @param rf_span Span parameter for RF baseline.
#' @param rf_maxit Max iterations for RF baseline.
#' @param column Measure column name for error context.
#' @param curve Curve index for error context.
#' @return A measure_tbl with baseline-corrected values.
#' @noRd
.sec_baseline_single <- function(
  x,
  left_frac,
  right_frac,
  method,
  rf_span = 2 / 3,
  rf_maxit = c(20, 20),
  column = NULL,
  curve = NULL
) {
  y <- x$value
  n <- length(y)

  # Handle edge cases
  if (n < 10) {
    cli::cli_warn(
      "Chromatogram has fewer than 10 points; returning unchanged."
    )
    return(x)
  }

  if (all(is.na(y))) {
    cli::cli_warn("Chromatogram is all NA; returning unchanged.")
    return(x)
  }

  # RF method uses robust local regression - no regions needed
  if (method == "rf") {
    baseline <- .sec_baseline_rf(y, rf_span = rf_span, rf_maxit = rf_maxit)
    x$value <- y - baseline
    return(x)
  }

  # Calculate region indices for other methods
  n_left <- max(1, floor(n * left_frac))
  n_right <- max(1, floor(n * right_frac))

  left_idx <- seq_len(n_left)
  right_idx <- seq(n - n_right + 1, n)

  left_idx <- left_idx[is.finite(y[left_idx])]
  right_idx <- right_idx[is.finite(y[right_idx])]
  missing_edges <- c("left", "right")[
    c(length(left_idx), length(right_idx)) == 0
  ]
  if (length(missing_edges) > 0) {
    context <- c(
      if (!is.null(column)) paste0("column `", column, "`"),
      if (!is.null(curve)) paste("curve", curve)
    )
    context <- if (length(context) > 0) {
      paste0(" (", paste(context, collapse = ", "), ")")
    } else {
      ""
    }
    cli::cli_abort(
      "SEC baseline{context}: no finite signal values in the {.and {missing_edges}} edge region{?s}.",
      class = "measure_sec_baseline_support",
      method = method,
      missing_edges = missing_edges,
      column = column,
      curve = curve
    )
  }

  fit <- NULL
  if (method == "spline") {
    # Smooth spline through baseline regions
    baseline_idx <- c(left_idx, right_idx)
    baseline_y <- y[baseline_idx]

    if (length(unique(baseline_idx)) < 4) {
      cli::cli_warn(
        "Too few valid baseline points for spline; using linear method."
      )
    } else {
      fit <- tryCatch(
        stats::smooth.spline(
          x = baseline_idx,
          y = baseline_y,
          spar = 0.8 # Strong smoothing
        ),
        error = function(e) {
          cli::cli_warn(
            "Spline fitting failed: {e$message}; using linear method."
          )
          return(NULL)
        }
      )
    }
  }

  baseline <- if (is.null(fit)) {
    .sec_baseline_interpolate(
      y,
      left_idx,
      right_idx,
      statistic = if (method == "median") stats::median else mean,
      warn_overlap = method != "spline"
    )
  } else {
    stats::predict(fit, x = seq_len(n))$y
  }

  # Subtract baseline
  x$value <- y - baseline
  x$value[!is.finite(y)] <- y[!is.finite(y)]
  x
}

.sec_baseline_interpolate <- function(
  y,
  left_idx,
  right_idx,
  statistic = mean,
  warn_overlap = TRUE
) {
  left_value <- statistic(y[left_idx])
  right_value <- statistic(y[right_idx])
  left_center <- statistic(left_idx)
  right_center <- statistic(right_idx)

  if (abs(right_center - left_center) < .Machine$double.eps) {
    if (warn_overlap) {
      cli::cli_warn("Baseline regions overlap; using constant baseline.")
    }
    return(rep(mean(c(left_value, right_value)), length(y)))
  }

  slope <- (right_value - left_value) / (right_center - left_center)
  intercept <- left_value - slope * left_center
  intercept + slope * seq_along(y)
}

#' Compute RF (Robust Fitting) baseline using IDPmisc::rfbaseline
#'
#' This function uses robust local regression to estimate the baseline.
#' It is particularly effective for complex baselines with curvature or drift.
#'
#' @param y Numeric vector of signal values.
#' @param rf_span Span parameter (fraction of points for local regression).
#' @param rf_maxit Max iterations for the two stages of robust fitting.
#' @return Numeric vector of baseline estimates.
#' @noRd
.sec_baseline_rf <- function(y, rf_span = 2 / 3, rf_maxit = c(20, 20)) {
  rlang::check_installed("IDPmisc", reason = "for RF baseline correction")

  # Handle edge case of flat or near-flat signal
  signal_range <- diff(range(y, na.rm = TRUE))
  signal_sd <- stats::sd(y, na.rm = TRUE)

  if (signal_range < 1e-6 || signal_sd < 1e-6) {
    return(rep(min(y, na.rm = TRUE), length(y)))
  }

  # Fit RF baseline
  baseline_fit <- tryCatch(
    IDPmisc::rfbaseline(
      x = seq_along(y),
      y = y,
      span = rf_span,
      maxit = rf_maxit
    )$fit,
    error = function(e) {
      cli::cli_warn(
        c(
          "RF baseline fitting failed: {e$message}",
          "i" = "Falling back to minimum value as baseline."
        )
      )
      rep(min(y, na.rm = TRUE), length(y))
    }
  )

  baseline_fit
}
