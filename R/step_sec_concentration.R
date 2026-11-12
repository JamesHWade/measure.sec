# ==============================================================================
# step_sec_concentration
#
# Convert detector signal to concentration
# ==============================================================================

#' Convert Detector Signal to Concentration
#'
#' `step_sec_concentration()` creates a *specification* of a recipe step that
#' converts detector signals to absolute concentration values using calibration
#' factors and injection parameters.
#'
#' @param recipe A recipe object.
#' @param measures Character vector of detector column names to convert.
#'   If `NULL`, will convert all measure columns.
#' @param detector Type of detector signal being converted:
#'   - `"ri"`: Refractive index detector (assumes dn/dc normalized)
#'   - `"uv"`: UV detector (assumes extinction coefficient normalized)
#'   - `"auto"`: Attempt to detect from column names
#' @param injection_volume Injection volume in uL. Required for absolute
#'   concentration calculation.
#' @param injection_mass Injected mass in mg. Alternative to using
#'   `sample_concentration` and `injection_volume`.
#' @param sample_concentration Sample concentration in mg/mL. Used with
#'   `injection_volume` to calculate injected mass.
#' @param flow_rate Flow rate in mL/min for peak area calculations.
#' @param concentration_units Output concentration units: `"mg/mL"` (default)
#'   or `"g/mL"`. Units are recorded on the output measure column when mass
#'   normalization is performed.
#' @param normalize_to_mass Logical. If `TRUE`, normalize the chromatogram
#'   so that the total area equals the injected mass. Default is `TRUE`.
#' @param role Role for generated columns.
#' @param trained Logical indicating if the step has been trained.
#' @param skip Logical. Should the step be skipped when baking?
#' @param id Unique step identifier.
#'
#' @return An updated recipe with the new step added.
#'
#' @details
#' This step converts detector response (after dn/dc or extinction coefficient
#' normalization) to absolute concentration. The conversion uses the known
#' injected mass to normalize the chromatogram area:
#'
#' \deqn{c(t) = \frac{S(t) \times m_{inj}}{\int S(t) \times F \times dt}}
#'
#' where:
#' - S(t) is the normalized detector signal
#' - m_inj is the injected mass
#' - F is the flow rate
#'
#' **Workflow for concentration determination:**
#' 1. Baseline correct the chromatogram
#' 2. Apply detector-specific normalization (dn/dc or extinction coefficient)
#' 3. Apply this step with known injection parameters
#' 4. Result: concentration at each elution point
#'
#' Mass normalization excludes negative signal points from the integrated area
#' and sets those output points to zero. When any points are excluded, the
#' bake reports their count and percentage so the baseline correction can be
#' reviewed.
#'
#' When a selected measure column contains negative points, baking emits one
#' warning with class `measure_sec_negative_signal_clipped`. The warning's
#' `column` field is the affected column name, and its `diagnostics` field is a
#' tibble with one row for every curve in the current bake data, including
#' curves without negative points. Its columns are `row` (the one-based row
#' position in the current bake data), `negative_points`, `points`,
#' `positive_area`, `negative_area`, and `excluded_fraction`. The two area
#' columns use the curve's signal units times minutes. Positive and negative
#' samples are clipped separately before trapezoidal integration on the
#' sampled axis; this uses the sampled clipped ordinates rather than
#' interpolating zero crossings. `excluded_fraction` is calculated separately
#' for each curve as `negative_area / (positive_area + negative_area)`, using a
#' scaled calculation to avoid overflow. It represents the negative share of
#' absolute sampled signal area, not a recovered-mass estimate. The warning
#' retains the total clipped-point count and adds the largest per-curve
#' negative share with its current bake row; ties report the first row. Since
#' recipes may bake training data during `prep()`, this warning can also occur
#' during preparation. No clipping diagnostics are emitted when
#' `normalize_to_mass = FALSE`.
#'
#' @family sec-detectors
#' @export
#'
#' @examples
#' \dontrun{
#' library(recipes)
#' library(measure)
#'
#' # Convert RI signal to concentration
#' rec <- recipe(~., data = sec_triple_detect) |>
#'   step_measure_input_long(ri_signal, location = vars(elution_time), col_name = "ri") |>
#'   step_sec_baseline() |>
#'   step_sec_ri(dn_dc = 0.185) |>
#'   step_sec_concentration(
#'     detector = "ri",
#'     injection_volume = 100,        # uL
#'     sample_concentration = 2.0,    # mg/mL
#'     flow_rate = 1.0                # mL/min
#'   ) |>
#'   prep()
#' }
step_sec_concentration <- function(
  recipe,
  measures = NULL,
  detector = c("ri", "uv", "auto"),
  injection_volume = NULL,
  injection_mass = NULL,
  sample_concentration = NULL,
  flow_rate = 1.0,
  concentration_units = "mg/mL",
  normalize_to_mass = TRUE,
  role = NA,
  trained = FALSE,
  skip = FALSE,
  id = recipes::rand_id("sec_concentration")
) {
  if (
    !is.logical(normalize_to_mass) ||
      length(normalize_to_mass) != 1 ||
      is.na(normalize_to_mass)
  ) {
    cli::cli_abort("{.arg normalize_to_mass} must be TRUE or FALSE.")
  }
  detector <- match.arg(detector)
  if (
    !is.character(concentration_units) ||
      length(concentration_units) != 1 ||
      !concentration_units %in% c("mg/mL", "g/mL")
  ) {
    cli::cli_abort(
      "{.arg concentration_units} must be {.val mg/mL} or {.val g/mL}."
    )
  }

  # Validate injection parameters
  if (normalize_to_mass) {
    if (
      is.null(injection_mass) &&
        (is.null(injection_volume) || is.null(sample_concentration))
    ) {
      cli::cli_abort(
        c(
          "For mass normalization, provide either:",
          "*" = "{.arg injection_mass}",
          "*" = "Both {.arg injection_volume} and {.arg sample_concentration}"
        )
      )
    }

    # Calculate injection mass if not provided
    if (is.null(injection_mass)) {
      .sec_validate_positive_scalar(injection_volume, "injection_volume")
      .sec_validate_positive_scalar(
        sample_concentration,
        "sample_concentration"
      )
      # injection_volume is in uL, sample_concentration in mg/mL
      # injection_mass in mg = volume (uL) * conc (mg/mL) / 1000
      injection_mass <- injection_volume * sample_concentration / 1000
    }
    if (
      !is.numeric(injection_mass) ||
        length(injection_mass) != 1 ||
        !is.finite(injection_mass) ||
        injection_mass <= 0
    ) {
      cli::cli_abort(
        "{.arg injection_mass} must be a positive finite mass in mg."
      )
    }
  }

  if (
    !is.numeric(flow_rate) ||
      length(flow_rate) != 1 ||
      !is.finite(flow_rate) ||
      flow_rate <= 0
  ) {
    cli::cli_abort("{.arg flow_rate} must be a positive number.")
  }

  recipes::add_step(
    recipe,
    step_sec_concentration_new(
      measures = measures,
      detector = detector,
      injection_volume = injection_volume,
      injection_mass = injection_mass,
      sample_concentration = sample_concentration,
      flow_rate = flow_rate,
      concentration_units = concentration_units,
      normalize_to_mass = normalize_to_mass,
      role = role,
      trained = trained,
      skip = skip,
      id = id
    )
  )
}

step_sec_concentration_new <- function(
  measures,
  detector,
  injection_volume,
  injection_mass,
  sample_concentration,
  flow_rate,
  concentration_units,
  normalize_to_mass,
  role,
  trained,
  skip,
  id
) {
  recipes::step(
    subclass = "sec_concentration",
    measures = measures,
    detector = detector,
    injection_volume = injection_volume,
    injection_mass = injection_mass,
    sample_concentration = sample_concentration,
    flow_rate = flow_rate,
    concentration_units = concentration_units,
    normalize_to_mass = normalize_to_mass,
    role = role,
    trained = trained,
    skip = skip,
    id = id
  )
}

.sec_normalize_concentration_measure <- function(
  m,
  injection_mass,
  flow_rate,
  concentration_units,
  row
) {
  location <- m$location
  value <- m$value
  n <- length(location)

  if (
    n < 2 ||
      !all(is.finite(location)) ||
      !all(is.finite(value)) ||
      !(all(diff(location) > 0) || all(diff(location) < 0))
  ) {
    cli::cli_abort(
      "Concentration normalization requires finite values on a strictly monotonic elution-time axis with at least two points."
    )
  }

  axis_widths <- abs(diff(location))
  if (!all(is.finite(axis_widths)) || !all(axis_widths > 0)) {
    cli::cli_abort(
      "Concentration normalization requires finite, positive axis widths."
    )
  }

  pos_value <- pmax(value, 0)
  neg_value <- pmax(-value, 0)
  positive_area <- sum(
    axis_widths *
      (utils::head(pos_value, -1) + utils::tail(pos_value, -1)) /
      2
  )
  negative_area <- sum(
    axis_widths *
      (utils::head(neg_value, -1) + utils::tail(neg_value, -1)) /
      2
  )

  if (!is.finite(positive_area)) {
    cli::cli_abort(
      "Cannot normalize a concentration signal with non-finite positive integrated area."
    )
  }
  if (positive_area <= 0) {
    cli::cli_abort(
      "Cannot normalize a concentration signal with nonpositive integrated area."
    )
  }
  if (!is.finite(negative_area) || negative_area < 0) {
    cli::cli_abort(
      "Concentration normalization requires a finite, nonnegative negative-signal area."
    )
  }

  denominator <- positive_area * flow_rate
  if (!is.finite(denominator) || denominator <= 0) {
    cli::cli_abort(
      "Concentration normalization requires a finite, positive integrated area-flow denominator."
    )
  }

  normalization_factor <- injection_mass / denominator
  if (!is.finite(normalization_factor) || normalization_factor <= 0) {
    cli::cli_abort(
      "Concentration normalization requires a finite, positive normalization factor."
    )
  }

  # Keep this evaluation order stable to preserve ordinary output values.
  m$value <- pos_value * injection_mass / denominator
  if (identical(concentration_units, "g/mL")) {
    m$value <- m$value / 1000
  }

  if (!all(is.finite(m$value))) {
    cli::cli_abort(
      "Concentration normalization produced non-finite output values."
    )
  }
  if (any(pos_value > 0) && !any(m$value > 0)) {
    cli::cli_abort(
      "Concentration normalization underflowed all positive signal values to zero."
    )
  }

  scale <- max(positive_area, negative_area)
  positive_scaled <- positive_area / scale
  negative_scaled <- negative_area / scale
  excluded_fraction <- negative_scaled / (positive_scaled + negative_scaled)

  list(
    measure = m,
    row = as.integer(row),
    negative_points = as.integer(sum(value < 0)),
    points = as.integer(n),
    positive_area = positive_area,
    negative_area = negative_area,
    excluded_fraction = excluded_fraction
  )
}

.sec_concentration_diagnostics <- function(processed) {
  tibble::tibble(
    row = purrr::map_int(processed, "row"),
    negative_points = purrr::map_int(processed, "negative_points"),
    points = purrr::map_int(processed, "points"),
    positive_area = purrr::map_dbl(processed, "positive_area"),
    negative_area = purrr::map_dbl(processed, "negative_area"),
    excluded_fraction = purrr::map_dbl(processed, "excluded_fraction")
  )
}

#' @export
prep.step_sec_concentration <- function(x, training, info = NULL, ...) {
  check_for_measure(training)

  # Find measure columns if not specified
  if (is.null(x$measures)) {
    measure_cols <- find_measure_cols(training)

    if (x$detector == "auto") {
      # Try to find RI or UV columns
      ri_cols <- measure_cols[grepl("ri", measure_cols, ignore.case = TRUE)]
      uv_cols <- measure_cols[grepl("uv", measure_cols, ignore.case = TRUE)]
      measures <- c(ri_cols, uv_cols)
      if (length(measures) == 0) {
        measures <- measure_cols
      }
    } else {
      measures <- measure_cols
    }
  } else {
    measures <- x$measures
  }

  step_sec_concentration_new(
    measures = measures,
    detector = x$detector,
    injection_volume = x$injection_volume,
    injection_mass = x$injection_mass,
    sample_concentration = x$sample_concentration,
    flow_rate = x$flow_rate,
    concentration_units = x$concentration_units,
    normalize_to_mass = x$normalize_to_mass,
    role = x$role,
    trained = TRUE,
    skip = x$skip,
    id = x$id
  )
}

#' @export
bake.step_sec_concentration <- function(object, new_data, ...) {
  measures <- object$measures
  injection_mass <- object$injection_mass
  flow_rate <- object$flow_rate
  normalize_to_mass <- object$normalize_to_mass
  concentration_units <- object$concentration_units

  if (!normalize_to_mass) {
    return(tibble::as_tibble(new_data))
  }

  for (col in measures) {
    processed <- purrr::map2(
      new_data[[col]],
      seq_along(new_data[[col]]),
      function(m, row) {
        if (!is.null(injection_mass)) {
          .sec_normalize_concentration_measure(
            m = m,
            injection_mass = injection_mass,
            flow_rate = flow_rate,
            concentration_units = concentration_units,
            row = row
          )
        } else {
          list(measure = m)
        }
      }
    )

    output <- new_measure_list(purrr::map(processed, "measure"))
    if (!is.null(injection_mass)) {
      attr(output, "sec_concentration_unit") <- concentration_units

      diagnostics <- .sec_concentration_diagnostics(processed)
      negative_count <- sum(diagnostics$negative_points)
      point_count <- sum(diagnostics$points)

      if (negative_count > 0) {
        max_row <- which.max(diagnostics$excluded_fraction)
        max_fraction <- diagnostics$excluded_fraction[[max_row]]
        cli::cli_warn(
          c(
            "Negative detector values were excluded during mass normalization.",
            "i" = "{negative_count} of {point_count} signal points ({round(100 * negative_count / point_count, 1)}%) were set to zero.",
            "i" = "Column {.val {col}}: maximum negative share of absolute sampled signal area was {round(100 * max_fraction, 1)}% in row {max_row}."
          ),
          class = "measure_sec_negative_signal_clipped",
          column = col,
          diagnostics = diagnostics
        )
      }
    }
    new_data[[col]] <- output
  }

  tibble::as_tibble(new_data)
}

#' @export
print.step_sec_concentration <- function(
  x,
  width = max(20, options()$width - 30),
  ...
) {
  title <- "SEC concentration conversion"
  if (!is.null(x$injection_mass)) {
    title <- paste0(
      title,
      " (",
      round(x$injection_mass * 1000, 1),
      " ug injected)"
    )
  }

  if (x$trained) {
    cols_str <- paste(x$measures, collapse = ", ")
    cat(title, " on ", cols_str, sep = "")
  } else {
    cat(title)
  }
  cat("\n")
  invisible(x)
}

#' @rdname tidy.step_sec
#' @export
#' @keywords internal
tidy.step_sec_concentration <- function(x, ...) {
  tibble::tibble(
    measures = list(x$measures),
    detector = x$detector,
    injection_mass = x$injection_mass %||% NA_real_,
    flow_rate = x$flow_rate,
    concentration_units = x$concentration_units,
    id = x$id
  )
}

#' @rdname required_pkgs.step_sec
#' @export
#' @keywords internal
required_pkgs.step_sec_concentration <- function(x, ...) {
  c("measure.sec", "measure")
}
