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

  for (col in measures) {
    output <- new_measure_list(
      purrr::map(new_data[[col]], function(m) {
        location <- m$location
        value <- m$value
        n <- length(location)

        if (!normalize_to_mass) {
          return(m)
        }
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

        if (
          normalize_to_mass && !is.null(injection_mass) && injection_mass > 0
        ) {
          # Ensure non-negative values for integration
          pos_value <- pmax(value, 0)
          total_area <- sum(
            abs(diff(location)) *
              (utils::head(pos_value, -1) + utils::tail(pos_value, -1)) /
              2
          )

          if (total_area > 0) {
            # Normalize so that integral of c(t) * F * dt = m_inj
            # c(t) = value * m_inj / (total_area * F)
            # But since total_area is in signal*min and F in mL/min,
            # we get c(t) in mg/mL if m_inj in mg

            m$value <- value * injection_mass / (total_area * flow_rate)
            if (identical(concentration_units, "g/mL")) {
              m$value <- m$value / 1000
            }
          } else {
            cli::cli_abort(
              "Cannot normalize a concentration signal with nonpositive integrated area."
            )
          }
        }

        m
      })
    )
    if (normalize_to_mass) {
      attr(output, "sec_concentration_unit") <- concentration_units
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
