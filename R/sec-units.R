.sec_validate_positive_scalar <- function(value, argument) {
  if (
    !is.numeric(value) ||
      length(value) != 1 ||
      !is.finite(value) ||
      value <= 0
  ) {
    cli::cli_abort("{.arg {argument}} must be a positive finite number.")
  }

  invisible(value)
}

.sec_validate_concentration_unit <- function(
  unit,
  argument = "concentration_unit",
  allow_null = TRUE
) {
  if (is.null(unit) && allow_null) {
    return(invisible(NULL))
  }

  if (
    !is.character(unit) ||
      length(unit) != 1 ||
      is.na(unit) ||
      !unit %in% c("mg/mL", "g/mL")
  ) {
    cli::cli_abort(
      "{.arg {argument}} must be {.val mg/mL} or {.val g/mL}."
    )
  }

  invisible(unit)
}

.sec_resolve_concentration_unit <- function(
  unit,
  concentration_curve,
  argument = "concentration_unit"
) {
  .sec_validate_concentration_unit(unit, argument)
  recorded_unit <- attr(
    concentration_curve,
    "sec_concentration_unit",
    exact = TRUE
  )
  .sec_validate_concentration_unit(recorded_unit, "recorded concentration unit")

  if (
    !is.null(unit) && !is.null(recorded_unit) && !identical(unit, recorded_unit)
  ) {
    cli::cli_abort(
      "{.arg {argument}} conflicts with the concentration curve's recorded units."
    )
  }

  resolved <- unit %||% recorded_unit
  if (is.null(resolved)) {
    cli::cli_abort(
      "Specify {.arg {argument}} for an untagged concentration curve."
    )
  }

  resolved
}

.sec_concentration_g_per_ml <- function(value, unit) {
  if (identical(unit, "mg/mL")) value / 1000 else value
}

.sec_assert_aligned_curves <- function(signal, concentration, detector) {
  signal_location <- signal$location
  concentration_location <- concentration$location
  signal_value <- signal$value
  concentration_value <- concentration$value

  if (
    length(signal_location) != length(signal_value) ||
      length(concentration_location) != length(concentration_value) ||
      !is.numeric(signal_location) ||
      !is.numeric(concentration_location) ||
      !is.numeric(signal_value) ||
      !is.numeric(concentration_value)
  ) {
    cli::cli_abort(
      "{.val {detector}} and concentration curves must have matching locations and values."
    )
  }
  if (
    !all(is.finite(signal_location)) ||
      !all(is.finite(concentration_location)) ||
      !isTRUE(all.equal(
        signal_location,
        concentration_location,
        tolerance = 1e-9,
        check.attributes = FALSE
      ))
  ) {
    cli::cli_abort(
      "{.val {detector}} and concentration curve locations must align."
    )
  }

  invisible(NULL)
}

.sec_single_angle_mw <- function(
  scattering,
  concentration,
  concentration_unit,
  optical_constant,
  calibration_constant,
  min_signal,
  detector
) {
  .sec_assert_aligned_curves(scattering, concentration, detector)

  concentration_g_per_ml <- .sec_concentration_g_per_ml(
    concentration$value,
    concentration_unit
  )
  scattering_value <- scattering$value
  max_concentration <- max(
    concentration_g_per_ml[is.finite(concentration_g_per_ml)],
    0
  )
  threshold <- min_signal * max_concentration

  mw <- rep(NA_real_, length(scattering_value))
  valid <- is.finite(scattering_value) &
    is.finite(concentration_g_per_ml) &
    scattering_value > 0 &
    concentration_g_per_ml > threshold

  if (any(valid)) {
    mw[valid] <- scattering_value[valid] *
      calibration_constant /
      (optical_constant * concentration_g_per_ml[valid])
  }

  measure::new_measure_tbl(
    location = scattering$location,
    value = mw
  )
}
