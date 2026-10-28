# ==============================================================================
# step_sec_detector_delay
#
# Correct inter-detector volume delays in multi-detector SEC
# ==============================================================================

#' Correct Inter-Detector Volume Delays
#'
#' `step_sec_detector_delay()` creates a *specification* of a recipe step that
#' corrects for volume delays between detectors in multi-detector SEC systems.
#'
#' @param recipe A recipe object.
#' @param reference Character name of the reference detector column (typically RI).
#'   All other detectors will be aligned to this reference.
#' @param targets Character vector of detector column names to shift. If `NULL`,
#'   all measure columns except the reference will be shifted.
#' @param delay_volumes Named numeric vector of delay volumes in mL.
#'   Names should match the target column names. Positive values indicate the
#'   detector sees the sample after the reference.
#' @param delay_times Named numeric vector of delay times in minutes.
#'   Alternative to `delay_volumes`. Requires `flow_rate` to be specified.
#' @param flow_rate Flow rate in mL/min. Required if using `delay_times`.
#' @param method Method for shifting signals:
#'   - `"shift"` (default): Simple index shift (fastest, slight edge effects)
#'   - `"interpolate"`: Linear interpolation (smoother, preserves signal shape)
#' @param role Not used by this step.
#' @param trained Logical indicating if the step has been trained.
#' @param skip Logical. Should the step be skipped when baking?
#' @param id Unique step identifier.
#'
#' @return An updated recipe with the new step added.
#'
#' @details
#' In multi-detector SEC systems, detectors are connected in series and separated
#' by tubing. This causes each detector to see the same analyte at different times.
#' For accurate molecular weight calculations that combine signals from multiple
#' detectors (e.g., RI + MALS for absolute MW), these delays must be corrected.
#'
#' **Typical detector order and delays:**
#' - UV detector: Often first, minimal delay

#' - RI detector: Common reference detector
#' - MALS detector: Often has 0.1-0.3 mL delay from RI
#' - Viscometer: May have 0.2-0.5 mL delay
#'
#' **Determining delay volumes:**
#' 1. Inject a narrow standard and record all detector signals
#' 2. Measure the time offset between peak maxima
#' 3. Convert to volume: delay_volume = time_offset × flow_rate
#'
#' @family sec-chromatography
#' @export
#'
#' @examples
#' \dontrun{
#' library(recipes)
#' library(measure)
#'
#' # Correct UV and MALS signals relative to RI
#' rec <- recipe(~., data = sec_triple_detect) |>
#'   step_measure_input_long(ri_signal, location = vars(elution_time), col_name = "ri") |>
#'   step_measure_input_long(uv_signal, location = vars(elution_time), col_name = "uv") |>
#'   step_measure_input_long(mals_signal, location = vars(elution_time), col_name = "mals") |>
#'   step_sec_detector_delay(
#'     reference = "ri",
#'     delay_volumes = c(uv = -0.05, mals = 0.15)
#'   ) |>
#'   prep()
#' }
step_sec_detector_delay <- function(
  recipe,
  reference = NULL,
  targets = NULL,

  delay_volumes = NULL,
  delay_times = NULL,
  flow_rate = 1.0,
  method = c("shift", "interpolate"),
  role = NA,
  trained = FALSE,
  skip = FALSE,
  id = recipes::rand_id("sec_detector_delay")
) {
  method <- match.arg(method)

  # Validate inputs

  if (is.null(reference)) {
    cli::cli_abort("{.arg reference} must be specified.")
  }

  if (is.null(delay_volumes) && is.null(delay_times)) {
    cli::cli_abort(
      "Either {.arg delay_volumes} or {.arg delay_times} must be specified."
    )
  }

  if (!is.null(delay_volumes) && !is.null(delay_times)) {
    cli::cli_abort(
      "Specify either {.arg delay_volumes} or {.arg delay_times}, not both."
    )
  }

  if (!is.null(delay_times) && is.null(flow_rate)) {
    cli::cli_abort(
      "{.arg flow_rate} is required when using {.arg delay_times}."
    )
  }

  # Convert delay_times to delay_volumes if provided
  if (!is.null(delay_times)) {
    delay_volumes <- delay_times * flow_rate
  }

  if (!is.numeric(delay_volumes) || is.null(names(delay_volumes))) {
    cli::cli_abort(
      "{.arg delay_volumes} must be a named numeric vector."
    )
  }

  recipes::add_step(
    recipe,
    step_sec_detector_delay_new(
      reference = reference,
      targets = targets,
      delay_volumes = delay_volumes,
      flow_rate = flow_rate,
      method = method,
      role = role,
      trained = trained,
      skip = skip,
      id = id
    )
  )
}

step_sec_detector_delay_new <- function(
  reference,
  targets,
  delay_volumes,
  flow_rate,
  method,
  role,
  trained,
  skip,
  id
) {
  recipes::step(
    subclass = "sec_detector_delay",
    reference = reference,
    targets = targets,
    delay_volumes = delay_volumes,
    flow_rate = flow_rate,
    method = method,
    role = role,
    trained = trained,
    skip = skip,
    id = id
  )
}

#' @export
prep.step_sec_detector_delay <- function(x, training, info = NULL, ...) {
  check_for_measure(training)

  # Find all measure columns
  measure_cols <- find_measure_cols(training)

  # Validate reference exists
  if (!x$reference %in% measure_cols) {
    cli::cli_abort(
      "Reference column {.val {x$reference}} not found in measure columns: {.val {measure_cols}}."
    )
  }

  # Determine target columns
  if (is.null(x$targets)) {
    targets <- setdiff(measure_cols, x$reference)
  } else {
    targets <- x$targets
    # Validate targets exist
    missing <- setdiff(targets, measure_cols)
    if (length(missing) > 0) {
      cli::cli_abort(
        "Target columns not found: {.val {missing}}."
      )
    }
  }

  # Validate all targets have delay volumes specified
  missing_delays <- setdiff(targets, names(x$delay_volumes))
  if (length(missing_delays) > 0) {
    cli::cli_abort(
      "Delay volumes not specified for: {.val {missing_delays}}."
    )
  }

  step_sec_detector_delay_new(
    reference = x$reference,
    targets = targets,
    delay_volumes = x$delay_volumes,
    flow_rate = x$flow_rate,
    method = x$method,
    role = x$role,
    trained = TRUE,
    skip = x$skip,
    id = x$id
  )
}

#' Shift a signal by a volume delay
#' @noRd
.shift_signal <- function(m, delay_volume, flow_rate, method) {
  location <- m$location
  value <- m$value
  n <- length(location)

  if (n < 2) {
    return(m)
  }

  # Calculate the location spacing (assume uniform)
  # For time-based data, this is time per point
  dt <- mean(diff(location))

  # Convert delay volume to delay in location units
  # delay_volume (mL) / flow_rate (mL/min) = delay (min)
  delay_loc <- delay_volume / flow_rate

  # Number of points to shift
  shift_points <- delay_loc / dt

  if (method == "shift") {
    # Simple index shift
    shift_idx <- round(shift_points)

    if (shift_idx == 0) {
      return(m)
    }

    new_value <- rep(NA_real_, n)

    if (shift_idx > 0) {
      # Positive delay: detector sees sample later
      # Shift values earlier (subtract from indices)
      if (shift_idx < n) {
        new_value[1:(n - shift_idx)] <- value[(shift_idx + 1):n]
      }
    } else {
      # Negative delay: detector sees sample earlier
      # Shift values later (add to indices)
      shift_idx <- abs(shift_idx)
      if (shift_idx < n) {
        new_value[(shift_idx + 1):n] <- value[1:(n - shift_idx)]
      }
    }

    m$value <- new_value
  } else if (method == "interpolate") {
    # Interpolation method
    # Create new location grid shifted by delay
    new_location <- location - delay_loc

    # Interpolate values at the new locations
    new_value <- stats::approx(
      x = location,
      y = value,
      xout = new_location,
      rule = 1 # NA for extrapolation
    )$y

    m$value <- new_value
  }

  m
}

#' @export
bake.step_sec_detector_delay <- function(object, new_data, ...) {
  reference <- object$reference
  targets <- object$targets
  delay_volumes <- object$delay_volumes
  flow_rate <- object$flow_rate
  method <- object$method

  # Apply delay correction to each target column
  for (target in targets) {
    delay <- delay_volumes[target]

    new_data[[target]] <- new_measure_list(
      purrr::map(new_data[[target]], function(m) {
        .shift_signal(m, delay, flow_rate, method)
      })
    )
  }

  tibble::as_tibble(new_data)
}

#' @export
print.step_sec_detector_delay <- function(
  x,
  width = max(20, options()$width - 30),
  ...
) {
  title <- paste0("SEC detector delay correction (ref: ", x$reference, ")")
  if (x$trained) {
    targets_str <- paste(x$targets, collapse = ", ")
    cat(title, " on ", targets_str, sep = "")
  } else {
    cat(title)
  }
  cat("\n")
  invisible(x)
}

#' @rdname tidy.step_sec
#' @export
#' @keywords internal
tidy.step_sec_detector_delay <- function(x, ...) {
  tibble::tibble(
    reference = x$reference,
    targets = list(x$targets),
    delay_volumes = list(x$delay_volumes),
    method = x$method,
    id = x$id
  )
}

#' @rdname required_pkgs.step_sec
#' @export
#' @keywords internal
required_pkgs.step_sec_detector_delay <- function(x, ...) {
  c("measure.sec", "measure")
}
