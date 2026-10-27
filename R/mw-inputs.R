# Helpers for pairing a concentration signal with a per-slice MW curve.
# The tag describes how a curve was generated; the locations are still checked
# for every sample at bake time.
.sec_tag_mw_curve <- function(curves, source, scale) {
  attr(curves, "sec_mw_curve") <- list(source = source, scale = scale)
  curves
}

.sec_resolve_mw_input <- function(
  training,
  measures,
  mw_col,
  mw_scale,
  calibration
) {
  available <- find_measure_cols(training)
  tags <- lapply(
    training[available],
    attr,
    which = "sec_mw_curve",
    exact = TRUE
  )

  if (is.null(measures)) {
    measures <- available[vapply(tags, is.null, logical(1))]
  }
  if (length(measures) != 1 || !measures %in% available) {
    cli::cli_abort(
      "Select exactly one concentration measure column with {.arg measures}."
    )
  }
  if (!is.null(mw_col) && !is.null(calibration)) {
    cli::cli_abort(
      "Specify either {.arg mw_col} or {.arg calibration}, not both."
    )
  }

  # A calibration step records both its input signal and output scale. Only an
  # unambiguous matching curve is selected automatically.
  if (is.null(mw_col) && is.null(calibration)) {
    matches <- names(tags)[vapply(
      tags,
      function(tag) {
        is.list(tag) && identical(tag$source, measures)
      },
      logical(1)
    )]
    if (length(matches) > 1) {
      cli::cli_abort(
        "Multiple molecular-weight curves match {.val {measures}}; specify {.arg mw_col}."
      )
    }
    if (length(matches) == 1) {
      mw_col <- matches[[1]]
    } else {
      uncalibrated <- vapply(
        training[available],
        function(curve) {
          identical(
            attr(curve, "sec_scattering_basis", exact = TRUE),
            "relative_uncalibrated"
          ) &&
            identical(
              attr(curve, "sec_scattering_source", exact = TRUE),
              measures
            )
        },
        logical(1)
      )
      if (any(uncalibrated)) {
        cli::cli_abort(
          "The matching scattering curve is uncalibrated and cannot provide molecular weights."
        )
      }
    }
  }

  if (!is.null(mw_col)) {
    if (
      !is.character(mw_col) || length(mw_col) != 1 || !mw_col %in% available
    ) {
      cli::cli_abort(
        "{.arg mw_col} must name one molecular-weight measure column."
      )
    }
    if (identical(mw_col, measures)) {
      cli::cli_abort(
        "{.arg mw_col} and the concentration measure must be different columns."
      )
    }
    if (
      identical(
        attr(training[[mw_col]], "sec_scattering_basis", exact = TRUE),
        "relative_uncalibrated"
      )
    ) {
      cli::cli_abort(
        "An uncalibrated scattering curve cannot provide molecular weights."
      )
    }
    tag <- attr(training[[mw_col]], "sec_mw_curve", exact = TRUE)
    if (is.null(mw_scale) && is.list(tag)) {
      mw_scale <- tag$scale
    }
    if (
      !is.character(mw_scale) ||
        length(mw_scale) != 1 ||
        is.na(mw_scale) ||
        !mw_scale %in% c("log10", "absolute")
    ) {
      cli::cli_abort(
        "Specify {.arg mw_scale} as {.val log10} or {.val absolute} for {.arg mw_col}."
      )
    }
    if (is.list(tag) && !identical(tag$scale, mw_scale)) {
      cli::cli_abort(
        "{.arg mw_scale} conflicts with the recorded scale of {.val {mw_col}}."
      )
    }
  } else if (!is.null(mw_scale)) {
    cli::cli_abort(
      "{.arg mw_scale} requires a molecular-weight curve in {.arg mw_col}."
    )
  }

  list(measures = measures, mw_col = mw_col, mw_scale = mw_scale)
}

.sec_mw_slices <- function(
  signal,
  mw_curve = NULL,
  mw_scale = NULL,
  calibration = NULL,
  integration_range = NULL
) {
  location <- signal$location
  value <- signal$value
  if (length(location) != length(value) || !all(is.finite(location))) {
    cli::cli_abort(
      "The concentration curve must have finite locations and matching values."
    )
  }
  if (
    length(location) > 1 &&
      !(all(diff(location) > 0) || all(diff(location) < 0))
  ) {
    cli::cli_abort(
      "Concentration curve locations must be strictly monotonic and unique."
    )
  }
  if (!is.null(mw_curve)) {
    if (
      !isTRUE(all.equal(
        location,
        mw_curve$location,
        tolerance = 1e-9,
        check.attributes = FALSE
      ))
    ) {
      cli::cli_abort(
        "Molecular-weight and concentration curve locations must align."
      )
    }
    mw_values <- mw_curve$value
    if (length(mw_values) != length(location)) {
      cli::cli_abort(
        "Molecular-weight and concentration curves must have the same length."
      )
    }
  } else {
    mw_values <- NULL
  }

  if (!is.null(integration_range)) {
    if (
      length(integration_range) != 2 ||
        !all(is.finite(integration_range)) ||
        integration_range[1] >= integration_range[2]
    ) {
      cli::cli_abort(
        "{.arg integration_range} must contain two increasing finite values."
      )
    }
    keep <- location >= integration_range[1] & location <= integration_range[2]
    location <- location[keep]
    value <- value[keep]
    if (!is.null(mw_values)) mw_values <- mw_values[keep]
  }
  if (!all(is.finite(value))) {
    cli::cli_abort(
      "The selected concentration signal contains missing or non-finite values."
    )
  }
  if (length(location) < 2) {
    return(list(
      location = location,
      signal = value,
      mw = numeric(),
      weight = numeric()
    ))
  }

  if (!is.null(mw_values)) {
    mw <- if (identical(mw_scale, "absolute")) mw_values else 10^mw_values
  } else if (is.null(calibration)) {
    # Backward-compatible mode: the signal's axis is already log10(MW).
    mw <- 10^location
  } else if (
    is.numeric(calibration) &&
      length(calibration) == 2 &&
      all(is.finite(calibration))
  ) {
    mw <- 10^(calibration[1] * location + calibration[2])
  } else if (identical(calibration, "auto")) {
    cli::cli_warn(
      "The {.val auto} calibration is synthetic and must not be treated as an experimental calibration."
    )
    mw <- 10^seq(7, 2, length.out = length(location))
  } else {
    cli::cli_abort(
      "{.arg calibration} must be NULL, a finite slope/intercept pair, or {.val auto}."
    )
  }
  bad_mw <- !is.finite(mw) | mw <= 0
  if (any(bad_mw & value > 0)) {
    cli::cli_abort(
      "Positive signal includes slices without a valid molecular-weight calibration; adjust the integration range or calibration."
    )
  }
  if (any(value < 0)) {
    cli::cli_warn(
      "Negative concentration values were clipped to zero for molecular-weight calculations."
    )
  }
  step <- abs(diff(location))
  width <- c(
    step[1] / 2,
    if (length(step) > 1) (step[-1] + step[-length(step)]) / 2,
    step[length(step)] / 2
  )
  list(
    location = location,
    signal = pmax(value, 0),
    mw = mw,
    weight = pmax(value, 0) * width
  )
}
