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

.sec_validate_nonzero_scalar <- function(value, argument) {
  if (
    !is.numeric(value) ||
      length(value) != 1 ||
      !is.finite(value) ||
      value == 0
  ) {
    cli::cli_abort("{.arg {argument}} must be a nonzero finite number.")
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

  valid_shape <- is.character(unit) &&
    length(unit) == 1L &&
    !is.object(unit) &&
    is.null(dim(unit)) &&
    all(names(attributes(unit)) %in% "names")

  if (
    !valid_shape ||
      is.na(unit) ||
      !unit %in% c("mg/mL", "g/mL")
  ) {
    cli::cli_abort(
      "{.arg {argument}} must be {.val mg/mL} or {.val g/mL}."
    )
  }

  invisible(unname(unit))
}

.sec_measure_value_metadata <- function(x) {
  metadata <- measure::measure_value_metadata(x)
  legacy_unit <- attr(x, "sec_concentration_unit", exact = TRUE)

  if (is.null(legacy_unit)) {
    return(metadata)
  }

  .sec_validate_concentration_unit(
    legacy_unit,
    "recorded concentration unit"
  )
  legacy_unit <- unname(legacy_unit)

  if (
    !is.null(metadata) && !identical(metadata$quantity, "mass_concentration")
  ) {
    cli::cli_abort(
      "Measurement value metadata quantity conflicts with the legacy {.field sec_concentration_unit} declaration."
    )
  }
  if (
    !is.null(metadata) &&
      !is.null(metadata$unit) &&
      !identical(unname(metadata$unit), legacy_unit)
  ) {
    cli::cli_abort(
      "Measurement value metadata conflicts with the legacy {.field sec_concentration_unit} declaration."
    )
  }

  list(quantity = "mass_concentration", unit = legacy_unit)
}

.sec_measure_value_basis <- function(x, argument = "measurement") {
  metadata <- .sec_measure_value_metadata(x)

  if (is.null(metadata)) {
    return("untagged")
  }
  if (identical(metadata$quantity, "mass_concentration")) {
    return("mass_concentration")
  }
  if (identical(metadata$quantity, "detector_response")) {
    return("detector_response")
  }

  cli::cli_abort(
    c(
      "{.arg {argument}} declares unsupported value quantity {.val {metadata$quantity}}.",
      "i" = "Use {.val detector_response} for native detector values or {.val mass_concentration} for calibrated concentration values."
    )
  )
}

.sec_tag_concentration <- function(output, unit) {
  .sec_validate_concentration_unit(
    unit,
    "concentration_unit",
    allow_null = FALSE
  )
  output <- measure::set_measure_value_metadata(
    output,
    quantity = "mass_concentration",
    unit = unname(unit)
  )
  attr(output, "sec_concentration_unit") <- unname(unit)
  output
}

.sec_concentration_metadata <- function(x, argument = "concentration curve") {
  metadata <- .sec_measure_value_metadata(x)

  if (
    !is.null(metadata) &&
      !identical(metadata$quantity, "mass_concentration")
  ) {
    cli::cli_abort(
      "{.arg {argument}} has value quantity {.val {metadata$quantity}}; expected {.val mass_concentration}."
    )
  }

  if (!is.null(metadata$unit)) {
    .sec_validate_concentration_unit(
      metadata$unit,
      "recorded concentration unit",
      allow_null = FALSE
    )
  }

  metadata
}

.sec_copy_measure_value_metadata <- function(x, output) {
  metadata <- .sec_measure_value_metadata(x)

  if (!is.null(metadata)) {
    output <- measure::set_measure_value_metadata(
      output,
      quantity = metadata$quantity,
      unit = metadata$unit
    )
  }
  legacy_unit <- attr(x, "sec_concentration_unit", exact = TRUE)
  if (!is.null(legacy_unit)) {
    attr(output, "sec_concentration_unit") <- legacy_unit
  }

  output
}

.sec_resolve_concentration_unit <- function(
  unit,
  concentration_curve,
  argument = "concentration_unit",
  require_recorded = FALSE
) {
  unit <- .sec_validate_concentration_unit(unit, argument)
  metadata <- .sec_concentration_metadata(concentration_curve, argument)
  recorded_unit <- if (is.null(metadata)) NULL else metadata$unit

  if (require_recorded && is.null(recorded_unit)) {
    cli::cli_abort(
      c(
        "The concentration curve no longer records its unit.",
        "i" = "The unit for {.arg {argument}} was inferred during training, so the baked curve must retain its unit tag."
      )
    )
  }

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
  unit <- .sec_validate_concentration_unit(unit, allow_null = FALSE)
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
