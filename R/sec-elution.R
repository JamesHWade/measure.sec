.sec_elution_units <- function(x) {
  valid_shape <- is.character(x) &&
    length(x) == 1L &&
    !is.object(x) &&
    is.null(dim(x)) &&
    all(names(attributes(x)) %in% "names")

  if (
    !valid_shape ||
      is.na(x) ||
      !x %in% c("min", "s", "mL", "uL")
  ) {
    cli::cli_abort(
      "{.arg elution_units} must be one plain value: {.val min}, {.val s}, {.val mL}, or {.val uL}."
    )
  }

  unname(x)
}

.sec_elution_units_from_step <- function(x) {
  if (!is.list(x)) {
    cli::cli_abort("A SEC step must be a list to resolve {.arg elution_units}.")
  }

  field <- which(names(x) == "elution_units")
  if (length(field) == 0L) {
    return("min")
  }
  if (length(field) != 1L) {
    cli::cli_abort(
      "A SEC step must contain at most one {.arg elution_units} field."
    )
  }

  .sec_elution_units(x[[field]])
}

.sec_elution_flow <- function(x, required = TRUE) {
  if (
    !is.logical(required) ||
      length(required) != 1L ||
      is.na(required) ||
      is.object(required) ||
      !is.null(dim(required))
  ) {
    cli::cli_abort("The internal flow requirement must be one logical value.")
  }

  if (is.null(x)) {
    if (!required) {
      return(NULL)
    }
    cli::cli_abort(
      "{.arg flow_rate} is required for this elution-unit conversion."
    )
  }

  valid_shape <- is.numeric(x) &&
    length(x) == 1L &&
    !is.object(x) &&
    is.null(dim(x)) &&
    all(names(attributes(x)) %in% "names")

  if (!valid_shape || !is.finite(x) || x <= 0) {
    cli::cli_abort(
      "{.arg flow_rate} must be a positive finite number in mL/min."
    )
  }

  x
}

.sec_elution_volume_factor <- function(units, flow) {
  units <- .sec_elution_units(units)
  flow <- .sec_elution_flow(
    flow,
    required = units %in% c("min", "s")
  )

  factor <- switch(
    units,
    min = flow,
    s = flow / 60,
    mL = 1,
    uL = 0.001
  )

  if (!is.finite(factor) || factor <= 0) {
    cli::cli_abort(
      "The mL-per-axis-unit conversion factor for {.arg elution_units} must be positive and finite."
    )
  }

  unname(factor)
}
