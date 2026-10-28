# ==============================================================================
# step_sec_integration_window
#
# Define integration window for MW calculations
# ==============================================================================

#' SEC/GPC Integration Window
#'
#' `step_sec_integration_window()` creates a *specification* of a recipe step
#' that defines the integration window (start and end x-axis bounds) for
#' molecular weight calculations. This step adds an `.integration_window`
#' column that downstream steps like [step_sec_mw_averages()] can use.
#'
#' @param recipe A recipe object. The step will be added to the sequence of
#'   operations for this recipe.
#' @param measures An optional character vector of measure column names to
#'   process. If `NULL` (the default), all measure columns will be processed.
#' @param start Numeric. Start of integration window (x-axis value, e.g., mL).
#'   If `NULL` and `auto_detect = TRUE`, determined automatically from data.
#' @param end Numeric. End of integration window (x-axis value, e.g., mL).
#'   If `NULL` and `auto_detect = TRUE`, determined automatically from data.
#' @param auto_detect Logical. If `TRUE` (default), automatically determine
#'   window bounds from the data when `start` or `end` is `NULL`.
#' @param signal_threshold Numeric between 0 and 1. When auto-detecting, the
#'   fraction of maximum signal to use as threshold for defining significant
#'   signal region. Default is `0.01` (1% of max).
#' @param extend_beyond_cal Numeric. Fraction of calibration range to extend
#'   beyond the maximum calibration point for capturing low MW species.
#'   Default is `0.5` (50% extension). Only used when calibration info is
#'   available via `calibration_range`.
#' @param calibration_range Optional numeric vector `c(min, max)` specifying

#'   the calibration range in x-axis units. When provided, auto-detection
#'   respects these bounds and applies `extend_beyond_cal`.
#' @param min_window_width Numeric. Minimum window width to ensure valid
#'   integration. Default is `1.0` (1 mL for SEC).
#' @param role Not used by this step since no new variables are created.
#' @param trained A logical to indicate if the quantities for preprocessing
#'   have been estimated.
#' @param skip A logical. Should the step be skipped when the recipe is baked?
#' @param id A character string that is unique to this step to identify it.
#'
#' @return An updated version of `recipe` with the new step added to the
#'   sequence of any existing operations. An `.integration_window` column
#'   will be added containing a tibble with `start` and `end` for each sample.
#'
#' @details
#' The integration window defines the x-axis region used for molecular weight
#' calculations. For SEC/GPC, this typically corresponds to elution volume.
#'
#' **Auto-detection algorithm:**
#' 1. Find significant signal region (above `signal_threshold` of max)
#' 2. Extend slightly beyond signal boundaries
#' 3. If `calibration_range` provided, constrain start to calibration minimum
#' 4. Allow extension up to `extend_beyond_cal` beyond calibration maximum
#'    to capture low MW species
#' 5. Ensure minimum window width
#'
#' **Output format:**
#' The `.integration_window` column contains a list of tibbles, one per row,
#' each with columns:
#' - `start`: Start of integration window
#' - `end`: End of integration window
#'
#' # Tidying
#'
#' When you [`tidy()`][recipes::tidy.recipe()] this step, a tibble with columns
#' `terms`, `start`, `end`, `auto_detect`, and `id` is returned.
#'
#' @seealso [step_sec_mw_averages()] which uses the integration window,
#'   [step_sec_conventional_cal()] for calibration.
#'
#' @family sec-chromatography
#' @export
#'
#' @examples
#' \dontrun{
#' library(recipes)
#' library(measure)
#'
#' # Auto-detect integration window
#' rec <- recipe(~., data = sec_data) |>
#'   step_measure_input_long(ri_signal, location = vars(elution_volume)) |>
#'   step_sec_integration_window() |>
#'   prep()
#'
#' # Explicit window bounds
#' rec <- recipe(~., data = sec_data) |>
#'   step_measure_input_long(ri_signal, location = vars(elution_volume)) |>
#'   step_sec_integration_window(start = 8.0, end = 18.0) |>
#'   prep()
#'
#' # With calibration constraints
#' rec <- recipe(~., data = sec_data) |>
#'   step_measure_input_long(ri_signal, location = vars(elution_volume)) |>
#'   step_sec_integration_window(
#'     calibration_range = c(9.0, 16.0),
#'     extend_beyond_cal = 0.5
#'   ) |>
#'   prep()
#' }
step_sec_integration_window <- function(
  recipe,
  measures = NULL,
  start = NULL,
  end = NULL,
  auto_detect = TRUE,
  signal_threshold = 0.01,
  extend_beyond_cal = 0.5,
  calibration_range = NULL,
  min_window_width = 1.0,
  role = NA,
  trained = FALSE,
  skip = FALSE,
  id = recipes::rand_id("sec_integration_window")
) {
  # Validate parameters

  if (!is.null(start) && (!is.numeric(start) || length(start) != 1)) {
    cli::cli_abort("{.arg start} must be a single numeric value or NULL.")
  }
  if (!is.null(end) && (!is.numeric(end) || length(end) != 1)) {
    cli::cli_abort("{.arg end} must be a single numeric value or NULL.")
  }
  if (!is.null(start) && !is.null(end) && start >= end) {
    cli::cli_abort("{.arg start} must be less than {.arg end}.")
  }
  if (!is.logical(auto_detect) || length(auto_detect) != 1) {
    cli::cli_abort("{.arg auto_detect} must be TRUE or FALSE.")
  }
  if (
    !is.numeric(signal_threshold) ||
      length(signal_threshold) != 1 ||
      signal_threshold <= 0 ||
      signal_threshold >= 1
  ) {
    cli::cli_abort(
      "{.arg signal_threshold} must be a number between 0 and 1 (exclusive)."
    )
  }
  if (
    !is.numeric(extend_beyond_cal) ||
      length(extend_beyond_cal) != 1 ||
      extend_beyond_cal < 0
  ) {
    cli::cli_abort("{.arg extend_beyond_cal} must be a non-negative number.")
  }
  if (!is.null(calibration_range)) {
    if (!is.numeric(calibration_range) || length(calibration_range) != 2) {
      cli::cli_abort(
        "{.arg calibration_range} must be a numeric vector of length 2."
      )
    }
    if (calibration_range[1] >= calibration_range[2]) {
      cli::cli_abort(
        "{.arg calibration_range} first element must be less than second."
      )
    }
  }
  if (
    !is.numeric(min_window_width) ||
      length(min_window_width) != 1 ||
      min_window_width <= 0
  ) {
    cli::cli_abort("{.arg min_window_width} must be a positive number.")
  }

  # Warn if both explicit bounds and auto_detect
  if (!is.null(start) && !is.null(end) && auto_detect) {
    cli::cli_inform(
      c(
        "i" = "Both {.arg start} and {.arg end} provided.",
        "i" = "{.arg auto_detect} will be ignored."
      )
    )
  }

  recipes::add_step(
    recipe,
    step_sec_integration_window_new(
      measures = measures,
      start = start,
      end = end,
      auto_detect = auto_detect,
      signal_threshold = signal_threshold,
      extend_beyond_cal = extend_beyond_cal,
      calibration_range = calibration_range,
      min_window_width = min_window_width,
      computed_start = NULL,
      computed_end = NULL,
      role = role,
      trained = trained,
      skip = skip,
      id = id
    )
  )
}

step_sec_integration_window_new <- function(
  measures,
  start,
  end,
  auto_detect,
  signal_threshold,
  extend_beyond_cal,
  calibration_range,
  min_window_width,
  computed_start,
  computed_end,
  role,
  trained,
  skip,
  id
) {
  recipes::step(
    subclass = "sec_integration_window",
    measures = measures,
    start = start,
    end = end,
    auto_detect = auto_detect,
    signal_threshold = signal_threshold,
    extend_beyond_cal = extend_beyond_cal,
    calibration_range = calibration_range,
    min_window_width = min_window_width,
    computed_start = computed_start,
    computed_end = computed_end,
    role = role,
    trained = trained,
    skip = skip,
    id = id
  )
}

#' @export
prep.step_sec_integration_window <- function(x, training, info = NULL, ...) {
  check_for_measure(training)

  # Resolve which columns to process
  if (is.null(x$measures)) {
    measure_cols <- find_measure_cols(training)
  } else {
    measure_cols <- x$measures
  }

  # Determine final window bounds
  if (!is.null(x$start) && !is.null(x$end)) {
    # Both explicitly provided
    computed_start <- x$start
    computed_end <- x$end
  } else if (x$auto_detect) {
    # Auto-detect from training data
    bounds <- .detect_integration_window(
      training = training,
      measure_cols = measure_cols,
      signal_threshold = x$signal_threshold,
      calibration_range = x$calibration_range,
      extend_beyond_cal = x$extend_beyond_cal,
      min_window_width = x$min_window_width
    )
    computed_start <- x$start %||% bounds$start
    computed_end <- x$end %||% bounds$end
  } else {
    # No auto-detect and missing bounds
    if (is.null(x$start) || is.null(x$end)) {
      cli::cli_abort(
        c(
          "Cannot determine integration window.",
          "i" = "Either provide both {.arg start} and {.arg end}, or set {.arg auto_detect = TRUE}."
        )
      )
    }
    computed_start <- x$start
    computed_end <- x$end
  }

  # Validate computed window
  if (computed_end - computed_start < x$min_window_width) {
    center <- (computed_start + computed_end) / 2
    half_width <- x$min_window_width / 2
    computed_start <- center - half_width
    computed_end <- center + half_width
    cli::cli_inform(
      c(
        "i" = "Integration window expanded to meet minimum width.",
        "i" = "Window: {round(computed_start, 2)} to {round(computed_end, 2)}"
      )
    )
  }

  step_sec_integration_window_new(
    measures = measure_cols,
    start = x$start,
    end = x$end,
    auto_detect = x$auto_detect,
    signal_threshold = x$signal_threshold,
    extend_beyond_cal = x$extend_beyond_cal,
    calibration_range = x$calibration_range,
    min_window_width = x$min_window_width,
    computed_start = computed_start,
    computed_end = computed_end,
    role = x$role,
    trained = TRUE,
    skip = x$skip,
    id = x$id
  )
}

#' @export
bake.step_sec_integration_window <- function(object, new_data, ...) {
  # Create integration window column
  n_rows <- nrow(new_data)

  # Each row gets the same window (computed during prep)
  window_list <- lapply(seq_len(n_rows), function(i) {
    tibble::tibble(
      start = object$computed_start,
      end = object$computed_end
    )
  })

  new_data[[".integration_window"]] <- window_list

  tibble::as_tibble(new_data)
}

#' @export
print.step_sec_integration_window <- function(
  x,
  width = max(20, options()$width - 30),
  ...
) {
  if (x$trained) {
    cat(
      glue::glue(
        "SEC integration window [{round(x$computed_start, 2)}, ",
        "{round(x$computed_end, 2)}]"
      )
    )
  } else {
    if (!is.null(x$start) && !is.null(x$end)) {
      cat(
        glue::glue(
          "SEC integration window [{x$start}, {x$end}] (explicit)"
        )
      )
    } else if (x$auto_detect) {
      cat("SEC integration window [auto-detect]")
    } else {
      cat("SEC integration window [not configured]")
    }
  }
  cat("\n")

  invisible(x)
}

#' @rdname tidy.step_sec
#' @export
#' @keywords internal
tidy.step_sec_integration_window <- function(x, ...) {
  if (recipes::is_trained(x)) {
    terms_val <- x$measures
    start_val <- x$computed_start
    end_val <- x$computed_end
  } else {
    terms_val <- "<all measure columns>"
    start_val <- x$start %||% NA_real_
    end_val <- x$end %||% NA_real_
  }

  tibble::tibble(
    terms = terms_val,
    start = start_val,
    end = end_val,
    auto_detect = x$auto_detect,
    extend_beyond_cal = x$extend_beyond_cal,
    id = x$id
  )
}

#' @rdname required_pkgs.step_sec
#' @export
#' @keywords internal
required_pkgs.step_sec_integration_window <- function(x, ...) {
  c("measure.sec", "measure")
}

# ==============================================================================
# Auto-detection helper
# ==============================================================================

#' Auto-detect Integration Window from Data
#'
#' Determines integration window bounds by analyzing signal in training data.
#'
#' @param training Training data tibble.
#' @param measure_cols Character vector of measure column names.
#' @param signal_threshold Fraction of max signal for threshold.
#' @param calibration_range Optional calibration range.
#' @param extend_beyond_cal Extension fraction beyond calibration.
#' @param min_window_width Minimum window width.
#'
#' @return List with `start` and `end` values.
#' @noRd
.detect_integration_window <- function(
  training,
  measure_cols,
  signal_threshold,
  calibration_range,
  extend_beyond_cal,
  min_window_width
) {
  # Collect all location and value data from measure columns
  all_locations <- numeric(0)
  all_values <- numeric(0)

  for (col in measure_cols) {
    if (col %in% names(training)) {
      measure_data <- training[[col]]
      for (m in measure_data) {
        if (!is.null(m) && nrow(m) > 0) {
          all_locations <- c(all_locations, m$location)
          all_values <- c(all_values, m$value)
        }
      }
    }
  }

  if (length(all_values) == 0 || all(is.na(all_values))) {
    cli::cli_warn(
      c(
        "No valid signal data found for auto-detection.",
        "i" = "Using default window based on data range."
      )
    )
    loc_range <- range(all_locations, na.rm = TRUE)
    return(list(start = loc_range[1], end = loc_range[2]))
  }

  # Remove NAs
  valid_idx <- !is.na(all_values) & !is.na(all_locations)
  locations <- all_locations[valid_idx]
  values <- all_values[valid_idx]

  # Find max signal
  max_val <- max(values, na.rm = TRUE)
  min_val <- min(values, na.rm = TRUE)

  # Handle negative signals (flip for analysis)
  if (abs(min_val) > abs(max_val)) {
    values <- -values
    max_val <- max(values, na.rm = TRUE)
  }

  # Find threshold
  threshold <- max_val * signal_threshold

  # Find significant signal region
  significant_idx <- values > threshold
  if (sum(significant_idx) < 2) {
    # Fall back to full range
    loc_range <- range(locations, na.rm = TRUE)
    return(list(start = loc_range[1], end = loc_range[2]))
  }

  significant_locs <- locations[significant_idx]
  sig_start <- min(significant_locs, na.rm = TRUE)
  sig_end <- max(significant_locs, na.rm = TRUE)

  # Extend slightly beyond signal boundaries (0.5 units on each side)
  suggested_start <- sig_start - 0.5
  suggested_end <- sig_end + 0.5

  # Apply calibration constraints if provided
  if (!is.null(calibration_range)) {
    cal_min <- calibration_range[1]
    cal_max <- calibration_range[2]
    cal_range_width <- cal_max - cal_min

    # Don't go below calibration minimum
    suggested_start <- max(suggested_start, cal_min)

    # Allow extension beyond calibration maximum
    extended_max <- cal_max + (extend_beyond_cal * cal_range_width)
    suggested_end <- min(suggested_end, extended_max)
  }

  # Ensure minimum width
  if (suggested_end - suggested_start < min_window_width) {
    center <- (suggested_start + suggested_end) / 2
    half_width <- min_window_width / 2
    suggested_start <- center - half_width
    suggested_end <- center + half_width
  }

  list(
    start = round(suggested_start, 3),
    end = round(suggested_end, 3)
  )
}
