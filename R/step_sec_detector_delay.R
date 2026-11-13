# ==============================================================================
# step_sec_detector_delay
#
# Correct inter-detector elution delays in multi-detector SEC
# ==============================================================================

#' Correct Inter-Detector Elution Delays
#'
#' \code{step_sec_detector_delay()} creates a specification of a recipe step
#' that aligns detector signals to a reference detector in multi-detector SEC
#' systems.
#'
#' @param recipe A recipe object.
#' @param reference Character name of the reference detector column. The
#'   reference is unchanged and cannot also be a target.
#' @param targets Character vector of detector columns to shift. If \code{NULL},
#'   all measure columns except the reference are shifted. \code{character()}
#'   selects no target columns.
#' @param delay_volumes Named numeric vector of delay volumes in mL. Names must
#'   identify target columns. Positive values indicate that the target detector
#'   sees the sample after the reference.
#' @param delay_times Named numeric vector of delay times in minutes. This is
#'   an alternative to \code{delay_volumes}; the values remain in minutes in
#'   the step object.
#' @param flow_rate Flow rate in mL/min. Required only when converting between
#'   time and volume axes. Any supplied value must be positive and finite.
#' @param method Method for shifting signals: \code{"shift"} (index shift with
#'   round-to-even quantization) or \code{"interpolate"} (piecewise-linear
#'   interpolation).
#' @param role Not used by this step.
#' @param trained Logical indicating if the step has been trained.
#' @param skip Logical. Should the step be skipped when baking?
#' @param id Unique step identifier.
#' @param elution_units Units of the measure \code{location} coordinates:
#'   \code{"min"} (minutes), \code{"s"} (seconds), \code{"mL"} (millilitres),
#'   or \code{"uL"} (microlitres). The default is \code{"min"}.
#'
#' @return An updated recipe with the new step added.
#'
#' @details
#' Detector delay declarations retain their own dimensions: \code{delay_times}
#' are always minutes and \code{delay_volumes} are always millilitres,
#' regardless of \code{elution_units}. The step converts the declared delay to
#' the location axis units during baking. Time-to-time and volume-to-volume
#' conversions do not need a flow rate. Converting between time and volume
#' needs a positive finite \code{flow_rate} in mL/min.
#'
#' A corrected target value at location \code{x} is evaluated at
#' \code{x + delay}. Thus, a positive delay moves target values toward lower
#' locations. The reference and every target curve in a current bake row must
#' have exactly the same numeric location grid. This step does not resample
#' mismatched grids.
#'
#' The \code{"shift"} method requires a nonzero-delay grid to be uniformly
#' spaced; adjacent widths may differ by at most 1e-8 relative to the first
#' width. Delays are rounded to the nearest whole index using R's round-to-even
#' rule. A zero delay is an exact pass-through and does not require uniform
#' spacing. The \code{"interpolate"} method accepts irregular monotonic grids,
#' returns \code{NA} outside the sampled range, and propagates missing values
#' without bridging an internal \code{NA} gap. Neither method changes the
#' location coordinates.
#'
#' \code{elution_units} is a caller declaration. This step does not infer units
#' from record metadata or add an axis-unit attribute to measure columns.
#'
#' @family sec-chromatography
#' @export
#'
#' @examples
#' \dontrun{
#' library(recipes)
#' library(measure)
#'
#' # Correct UV and MALS signals relative to RI on a minute axis
#' rec <- recipe(~., data = sec_triple_detect) |>
#'   step_measure_input_long(ri_signal, location = vars(elution_time), col_name = "ri") |>
#'   step_measure_input_long(uv_signal, location = vars(elution_time), col_name = "uv") |>
#'   step_measure_input_long(mals_signal, location = vars(elution_time), col_name = "mals") |>
#'   step_sec_detector_delay(
#'     reference = "ri",
#'     delay_volumes = c(uv = -0.05, mals = 0.15),
#'     elution_units = "min"
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
  id = recipes::rand_id("sec_detector_delay"),
  elution_units = "min"
) {
  method <- match.arg(method)

  state <- .sec_delay_resolve_fields(
    list(
      reference = reference,
      targets = targets,
      delay_volumes = delay_volumes,
      delay_times = delay_times,
      flow_rate = flow_rate,
      method = method,
      elution_units = elution_units
    )
  )

  recipes::add_step(
    recipe,
    step_sec_detector_delay_new(
      reference = state$reference,
      targets = state$targets,
      delay_volumes = state$delay_volumes,
      delay_times = state$delay_times,
      flow_rate = state$flow_rate,
      method = state$method,
      elution_units = state$elution_units,
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
  delay_times,
  flow_rate,
  method,
  elution_units,
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
    delay_times = delay_times,
    flow_rate = flow_rate,
    method = method,
    elution_units = elution_units,
    role = role,
    trained = trained,
    skip = skip,
    id = id
  )
}

.sec_delay_named_values <- function(x, arg) {
  valid_shape <- is.numeric(x) &&
    is.atomic(x) &&
    length(x) > 0L &&
    is.null(dim(x)) &&
    !is.object(x) &&
    all(names(attributes(x)) %in% "names")

  if (!valid_shape || is.null(names(x))) {
    cli::cli_abort("{.arg {arg}} must be a nonempty named numeric vector.")
  }

  x_names <- names(x)
  if (
    length(x_names) != length(x) ||
      anyNA(x_names) ||
      !all(nzchar(trimws(x_names))) ||
      anyDuplicated(x_names)
  ) {
    cli::cli_abort(
      "Names in {.arg {arg}} must be unique, nonblank, and nonmissing."
    )
  }

  if (anyNA(x) || !all(is.finite(x))) {
    cli::cli_abort("Values in {.arg {arg}} must be finite and nonmissing.")
  }

  x
}

.sec_delay_reference <- function(x) {
  valid <- is.character(x) &&
    length(x) == 1L &&
    !is.object(x) &&
    is.null(dim(x)) &&
    !is.na(x) &&
    nzchar(trimws(x))

  if (!valid) {
    cli::cli_abort("{.arg reference} must be one nonblank column name.")
  }

  unname(x)
}

.sec_delay_targets <- function(x, reference) {
  if (is.null(x)) {
    return(NULL)
  }

  valid <- is.character(x) &&
    !is.object(x) &&
    is.null(dim(x)) &&
    !anyNA(x) &&
    all(nzchar(trimws(x))) &&
    !anyDuplicated(x)

  if (!valid) {
    cli::cli_abort(
      "{.arg targets} must contain unique, nonblank measure column names."
    )
  }
  if (reference %in% x) {
    cli::cli_abort("The {.arg reference} column cannot also be a target.")
  }

  unname(x)
}

.sec_delay_method <- function(x) {
  valid <- is.character(x) &&
    length(x) == 1L &&
    !is.object(x) &&
    is.null(dim(x)) &&
    !is.na(x) &&
    x %in% c("shift", "interpolate")

  if (!valid) {
    cli::cli_abort("{.arg method} must be {.val shift} or {.val interpolate}.")
  }

  unname(x)
}

.sec_delay_flow_required <- function(delay_times, elution_units) {
  if (!is.null(delay_times)) {
    return(elution_units %in% c("mL", "uL"))
  }

  elution_units %in% c("min", "s")
}

.sec_delay_checked_conversion <- function(value, declared, arg) {
  if (!all(is.finite(value))) {
    cli::cli_abort(
      "Converting {.arg {arg}} to {.arg elution_units} produced a nonfinite delay."
    )
  }
  if (any(declared != 0 & value == 0)) {
    cli::cli_abort(
      "Converting a nonzero {.arg {arg}} to {.arg elution_units} underflowed to zero."
    )
  }

  value
}

.sec_delay_convert <- function(
  delay_volumes,
  delay_times,
  elution_units,
  flow_rate
) {
  if (!is.null(delay_times)) {
    if (elution_units == "min") {
      return(delay_times)
    }
    if (elution_units == "s") {
      return(.sec_delay_checked_conversion(
        delay_times * 60,
        delay_times,
        "delay_times"
      ))
    }

    factor <- .sec_elution_volume_factor(elution_units, flow_rate)
    volume_delay <- .sec_delay_checked_conversion(
      delay_times * flow_rate,
      delay_times,
      "delay_times"
    )
    return(.sec_delay_checked_conversion(
      volume_delay / factor,
      delay_times,
      "delay_times"
    ))
  }

  factor <- .sec_elution_volume_factor(elution_units, flow_rate)
  .sec_delay_checked_conversion(
    delay_volumes / factor,
    delay_volumes,
    "delay_volumes"
  )
}

.sec_delay_resolve_fields <- function(x) {
  if (!is.list(x)) {
    cli::cli_abort("A SEC detector delay step must be a list.")
  }

  physics_fields <- c(
    "reference",
    "targets",
    "delay_times",
    "delay_volumes",
    "flow_rate",
    "method"
  )
  field_names <- names(x)
  duplicate_fields <- physics_fields[vapply(
    physics_fields,
    function(field) {
      sum(!is.na(field_names) & field_names == field) > 1L
    },
    logical(1)
  )]
  if (length(duplicate_fields) > 0L) {
    cli::cli_abort(
      "SEC detector delay state contains duplicate configuration fields: {.field {duplicate_fields}}."
    )
  }

  reference <- .sec_delay_reference(x[["reference", exact = TRUE]])
  targets <- .sec_delay_targets(
    x[["targets", exact = TRUE]],
    reference
  )
  delay_volumes <- x[["delay_volumes", exact = TRUE]]
  delay_times <- x[["delay_times", exact = TRUE]]

  if (is.null(delay_volumes) == is.null(delay_times)) {
    cli::cli_abort(
      "Specify exactly one nonempty {.arg delay_volumes} or {.arg delay_times} vector."
    )
  }
  if (!is.null(delay_volumes)) {
    delay_volumes <- .sec_delay_named_values(delay_volumes, "delay_volumes")
  }
  if (!is.null(delay_times)) {
    delay_times <- .sec_delay_named_values(delay_times, "delay_times")
  }

  elution_units <- .sec_elution_units_from_step(x)
  method <- .sec_delay_method(x[["method", exact = TRUE]])
  flow_rate <- x[["flow_rate", exact = TRUE]]
  flow_rate <- .sec_elution_flow(
    flow_rate,
    required = .sec_delay_flow_required(delay_times, elution_units)
  )
  delays <- .sec_delay_convert(
    delay_volumes = delay_volumes,
    delay_times = delay_times,
    elution_units = elution_units,
    flow_rate = flow_rate
  )
  if (!is.null(targets)) {
    missing_delays <- setdiff(targets, names(delays))
    if (length(missing_delays) > 0L) {
      cli::cli_abort(
        "Delays not specified for target columns: {.val {missing_delays}}."
      )
    }
  }

  list(
    reference = reference,
    targets = targets,
    delay_volumes = delay_volumes,
    delay_times = delay_times,
    delays = delays,
    flow_rate = flow_rate,
    method = method,
    elution_units = elution_units,
    role = x[["role", exact = TRUE]],
    trained = x[["trained", exact = TRUE]],
    skip = x[["skip", exact = TRUE]],
    id = x[["id", exact = TRUE]]
  )
}

#' @export
prep.step_sec_detector_delay <- function(x, training, info = NULL, ...) {
  check_for_measure(training)
  state <- .sec_delay_resolve_fields(x)

  measure_cols <- find_measure_cols(training)
  if (!state$reference %in% measure_cols) {
    cli::cli_abort(
      "Reference column {.val {state$reference}} not found in measure columns: {.val {measure_cols}}."
    )
  }

  targets <- state$targets
  if (is.null(targets)) {
    targets <- setdiff(measure_cols, state$reference)
  }

  missing <- setdiff(targets, measure_cols)
  if (length(missing) > 0L) {
    cli::cli_abort("Target columns not found: {.val {missing}}.")
  }
  missing_delays <- setdiff(targets, names(state$delays))
  if (length(missing_delays) > 0L) {
    cli::cli_abort(
      "Delays not specified for target columns: {.val {missing_delays}}."
    )
  }

  step_sec_detector_delay_new(
    reference = state$reference,
    targets = targets,
    delay_volumes = state$delay_volumes,
    delay_times = state$delay_times,
    flow_rate = state$flow_rate,
    method = state$method,
    elution_units = state$elution_units,
    role = state$role,
    trained = TRUE,
    skip = state$skip,
    id = state$id
  )
}

.sec_delay_curve <- function(m) {
  if (!is.data.frame(m)) {
    cli::cli_abort("Each detector curve must be a measure table.")
  }
  location <- m[["location", exact = TRUE]]
  value <- m[["value", exact = TRUE]]
  if (
    !is.numeric(location) ||
      !is.atomic(location) ||
      !is.null(dim(location)) ||
      !is.numeric(value) ||
      !is.atomic(value) ||
      !is.null(dim(value)) ||
      length(location) != length(value)
  ) {
    cli::cli_abort(
      "Each detector curve must have equal-length numeric {.field location} and {.field value} vectors."
    )
  }

  location <- unname(as.numeric(location))
  value <- unname(as.numeric(value))
  if (!all(is.finite(location))) {
    cli::cli_abort("Detector locations must be finite numeric values.")
  }
  if (any(is.nan(value)) || any(is.infinite(value))) {
    cli::cli_abort("Detector values must be finite or ordinary NA values.")
  }

  widths <- diff(location)
  if (
    !all(is.finite(widths)) ||
      any(widths == 0) ||
      (length(widths) > 0L && !all(widths > 0) && !all(widths < 0))
  ) {
    cli::cli_abort(
      "Detector locations must be strictly monotonic with finite adjacent widths."
    )
  }

  list(location = location, value = value)
}

.sec_delay_uniform_spacing <- function(location) {
  widths <- diff(location)
  relative_widths <- abs(widths / widths[[1L]] - 1)
  if (!all(is.finite(relative_widths)) || any(relative_widths > 1e-8)) {
    cli::cli_abort(
      "A nonzero {.val shift} delay requires uniformly spaced locations (relative width tolerance 1e-8)."
    )
  }

  widths[[1L]]
}

.sec_delay_signal <- function(m, delay, method) {
  curve <- .sec_delay_curve(m)
  location <- curve$location
  value <- curve$value
  n <- length(location)

  if (delay == 0) {
    return(m)
  }
  if (n < 2L) {
    cli::cli_abort(
      "A nonzero detector delay requires at least two sampled locations."
    )
  }

  if (method == "shift") {
    spacing <- .sec_delay_uniform_spacing(location)
    shift_points <- delay / spacing
    if (!is.finite(shift_points) || abs(shift_points) >= n) {
      m$value <- rep(NA_real_, n)
      return(m)
    }

    shift_index <- round(shift_points)
    if (shift_index == 0) {
      return(m)
    }
    if (abs(shift_index) >= n) {
      m$value <- rep(NA_real_, n)
      return(m)
    }

    shift_index <- as.integer(shift_index)
    shifted <- rep(NA_real_, n)
    if (shift_index > 0L) {
      out_index <- seq_len(n - shift_index)
      shifted[out_index] <- value[out_index + shift_index]
    } else {
      in_index <- seq_len(n + shift_index)
      out_index <- in_index - shift_index
      shifted[out_index] <- value[in_index]
    }
    m$value <- shifted
    return(m)
  }

  query <- location + delay
  if (!all(is.finite(query))) {
    cli::cli_abort(
      "Adding the converted delay to detector locations produced nonfinite interpolation coordinates."
    )
  }
  if (any(query == location)) {
    cli::cli_abort(
      "The nonzero delay is not representable at the detector location scale; interpolation coordinates did not change."
    )
  }
  query_widths <- diff(query)
  location_direction <- sign(location[[2L]] - location[[1L]])
  if (
    !all(is.finite(query_widths)) ||
      any(query_widths == 0) ||
      any(sign(query_widths) != location_direction)
  ) {
    cli::cli_abort(
      "Adding the converted delay did not preserve a strictly monotonic interpolation grid."
    )
  }

  if (location[[1L]] > location[[2L]]) {
    location <- rev(location)
    value <- rev(value)
  }

  interpolated <- stats::approx(
    x = location,
    y = value,
    xout = query,
    method = "linear",
    rule = 1,
    ties = "ordered",
    na.rm = FALSE
  )$y
  if (any(is.infinite(interpolated)) || any(is.nan(interpolated))) {
    cli::cli_abort(
      "Interpolation produced {.val Inf} or {.val NaN}; interpolated values must be finite or ordinary NA."
    )
  }

  m$value <- interpolated
  m
}

.sec_delay_validate_current_data <- function(data, reference, targets) {
  measure_cols <- find_measure_cols(data)
  needed <- unique(c(reference, targets))
  missing <- setdiff(needed, measure_cols)
  if (length(missing) > 0L) {
    cli::cli_abort(
      "Reference or target measure columns are missing from the current bake data: {.val {missing}}."
    )
  }

  n_rows <- nrow(data)
  for (column in needed) {
    if (length(data[[column]]) != n_rows) {
      cli::cli_abort(
        "Measure column {.val {column}} does not match the current bake row count."
      )
    }
  }

  invisible(measure_cols)
}

#' @export
bake.step_sec_detector_delay <- function(object, new_data, ...) {
  state <- .sec_delay_resolve_fields(object)
  if (is.null(state$targets)) {
    cli::cli_abort(
      "The detector delay step must be trained before it can be baked."
    )
  }

  .sec_delay_validate_current_data(new_data, state$reference, state$targets)
  reference_curves <- new_data[[state$reference]]

  for (target in state$targets) {
    input_curves <- new_data[[target]]
    target_delay <- state$delays[[match(target, names(state$delays))]]
    output_curves <- purrr::map2(
      input_curves,
      reference_curves,
      function(target_curve, reference_curve) {
        target_data <- .sec_delay_curve(target_curve)
        reference_data <- .sec_delay_curve(reference_curve)
        if (
          length(target_data$location) != length(reference_data$location) ||
            any(target_data$location != reference_data$location)
        ) {
          cli::cli_abort(
            "Reference and target detector locations must align exactly for each current bake row; resampling is not performed."
          )
        }
        .sec_delay_signal(target_curve, target_delay, state$method)
      }
    )
    output_curve <- new_measure_list(output_curves)
    concentration_unit <- attr(
      input_curves,
      "sec_concentration_unit",
      exact = TRUE
    )
    if (!is.null(concentration_unit)) {
      attr(output_curve, "sec_concentration_unit") <- concentration_unit
    }
    new_data[[target]] <- output_curve
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
  state <- .sec_delay_resolve_fields(x)
  tibble::tibble(
    reference = state$reference,
    targets = list(state$targets),
    delay_volumes = list(state$delay_volumes),
    delay_times = list(state$delay_times),
    elution_units = state$elution_units,
    method = state$method,
    id = state$id
  )
}

#' @rdname required_pkgs.step_sec
#' @export
#' @keywords internal
required_pkgs.step_sec_detector_delay <- function(x, ...) {
  c("measure.sec", "measure")
}
