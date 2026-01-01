# ==============================================================================
# step_sec_intrinsic_visc
#
# Calculate intrinsic viscosity from specific viscosity and concentration
# ==============================================================================

#' Calculate Intrinsic Viscosity for SEC
#'
#' `step_sec_intrinsic_visc()` creates a *specification* of a recipe step that
#' calculates intrinsic viscosity (\[eta\]) from specific viscosity and
#' concentration at each elution point.
#'
#' @param recipe A recipe object.
#' @param specific_visc_col Name of the specific viscosity measure column
#'   (from `step_sec_viscometer()`).
#' @param concentration_col Name of the concentration measure column.
#' @param output_col Name for the intrinsic viscosity output column.
#'   Default is `"intrinsic_visc"`.
#' @param min_concentration Minimum concentration threshold below which
#'   intrinsic viscosity is set to NA. Default is 1e-6 mg/mL.
#' @param units Output units for intrinsic viscosity. Default is `"dL/g"`.
#'   Common alternatives are `"mL/g"` (multiply by 100).
#' @param role Role for generated columns.
#' @param trained Logical indicating if the step has been trained.
#' @param skip Logical. Should the step be skipped when baking?
#' @param id Unique step identifier.
#'
#' @return An updated recipe with the new step added, containing intrinsic
#'   viscosity at each elution point.
#'
#' @details
#' Intrinsic viscosity is defined as the limit of reduced viscosity at
#' infinite dilution:
#'
#' \deqn{[\eta] = \lim_{c \to 0} \frac{\eta_{sp}}{c}}
#'
#' In SEC, each elution slice has very low concentration, so the approximation
#' \[eta\] = eta_sp / c is valid.
#'
#' **Applications of Intrinsic Viscosity:**
#' \itemize{
#'   \item Universal Calibration: log(\[eta\] * M) is linear with retention volume,
#'     allowing calibration transfer between polymer types
#'   \item Mark-Houwink Equation: \[eta\] = K * M^a, where K and a are polymer-
#'     and solvent-specific constants
#'   \item Branching Analysis: g' = \[eta\]_branched / \[eta\]_linear provides
#'     information about long-chain branching
#'   \item Polymer Conformation: The scaling exponent in \[eta\] vs M reveals

#'     chain conformation (coil, rod, sphere)
#' }
#'
#' **Typical Intrinsic Viscosity Values (dL/g):**
#' \itemize{
#'   \item Polystyrene in THF (MW 100k): ~0.5
#'   \item PEG in water (MW 10k): ~0.2
#'   \item Proteins: 0.03-0.05 (globular), 0.2-1.0 (denatured)
#' }
#'
#' @family sec-detectors
#' @seealso [step_sec_viscometer()] for specific viscosity calculation
#' @export
#'
#' @examples
#' \dontrun{
#' library(recipes)
#' library(measure)
#'
#' # Calculate intrinsic viscosity
#' rec <- recipe(~., data = sec_visc_data) |>
#'   step_measure_input_long(dp_signal, location = vars(elution_time), col_name = "dp") |>
#'   step_measure_input_long(ri_signal, location = vars(elution_time), col_name = "ri") |>
#'   step_sec_baseline() |>
#'   step_sec_viscometer(dp_col = "dp") |>
#'   step_sec_ri(dn_dc = 0.185) |>
#'   step_sec_concentration(detector = "ri", injection_mass = 0.2, flow_rate = 1.0) |>
#'   step_sec_intrinsic_visc(
#'     specific_visc_col = "specific_visc",
#'     concentration_col = "ri"
#'   ) |>
#'   prep()
#' }
step_sec_intrinsic_visc <- function(
  recipe,
  specific_visc_col = NULL,
  concentration_col = NULL,
  output_col = "intrinsic_visc",
  min_concentration = 1e-6,
  units = c("dL/g", "mL/g"),
  role = NA,
  trained = FALSE,
  skip = FALSE,
  id = recipes::rand_id("sec_intrinsic_visc")
) {
  units <- match.arg(units)

  if (!is.numeric(min_concentration) || min_concentration < 0) {
    cli::cli_abort("{.arg min_concentration} must be a non-negative number.")
  }

  recipes::add_step(
    recipe,
    step_sec_intrinsic_visc_new(
      specific_visc_col = specific_visc_col,
      concentration_col = concentration_col,
      output_col = output_col,
      min_concentration = min_concentration,
      units = units,
      role = role,
      trained = trained,
      skip = skip,
      id = id
    )
  )
}

step_sec_intrinsic_visc_new <- function(
  specific_visc_col,
  concentration_col,
  output_col,
  min_concentration,
  units,
  role,
  trained,
  skip,
  id
) {
  recipes::step(
    subclass = "sec_intrinsic_visc",
    specific_visc_col = specific_visc_col,
    concentration_col = concentration_col,
    output_col = output_col,
    min_concentration = min_concentration,
    units = units,
    role = role,
    trained = trained,
    skip = skip,
    id = id
  )
}

#' @export
prep.step_sec_intrinsic_visc <- function(x, training, info = NULL, ...) {
  check_for_measure(training)

  measure_cols <- find_measure_cols(training)

  # Find specific viscosity column if not specified
  if (is.null(x$specific_visc_col)) {
    visc_cols <- measure_cols[grepl(
      "visc|eta",
      measure_cols,
      ignore.case = TRUE
    )]
    if (length(visc_cols) == 0) {
      cli::cli_abort(
        "No specific viscosity column found. Specify {.arg specific_visc_col} explicitly."
      )
    }
    specific_visc_col <- visc_cols[1]
  } else {
    specific_visc_col <- x$specific_visc_col
    if (!specific_visc_col %in% measure_cols) {
      cli::cli_abort(
        "Specific viscosity column {.val {specific_visc_col}} not found in measure columns."
      )
    }
  }

  # Find concentration column if not specified
  if (is.null(x$concentration_col)) {
    conc_cols <- measure_cols[grepl(
      "ri|conc",
      measure_cols,
      ignore.case = TRUE
    )]
    if (length(conc_cols) == 0) {
      cli::cli_abort(
        "No concentration column found. Specify {.arg concentration_col} explicitly."
      )
    }
    concentration_col <- conc_cols[1]
  } else {
    concentration_col <- x$concentration_col
    if (!concentration_col %in% measure_cols) {
      cli::cli_abort(
        "Concentration column {.val {concentration_col}} not found in measure columns."
      )
    }
  }

  step_sec_intrinsic_visc_new(
    specific_visc_col = specific_visc_col,
    concentration_col = concentration_col,
    output_col = x$output_col,
    min_concentration = x$min_concentration,
    units = x$units,
    role = x$role,
    trained = TRUE,
    skip = x$skip,
    id = x$id
  )
}

#' @export
bake.step_sec_intrinsic_visc <- function(object, new_data, ...) {
  specific_visc_col <- object$specific_visc_col
  concentration_col <- object$concentration_col
  output_col <- object$output_col
  min_concentration <- object$min_concentration
  units <- object$units

  # Unit conversion factor
  # Default is dL/g; if mL/g requested, multiply by 100
  unit_factor <- if (units == "mL/g") 100 else 1

  # Calculate intrinsic viscosity for each sample
  iv_list <- purrr::map2(
    new_data[[specific_visc_col]],
    new_data[[concentration_col]],
    function(visc_m, conc_m) {
      visc_val <- visc_m$value
      conc_val <- conc_m$value
      location <- visc_m$location

      # Calculate \[eta\] = eta_sp / c
      iv <- rep(NA_real_, length(visc_val))

      valid <- !is.na(visc_val) &
        !is.na(conc_val) &
        conc_val > min_concentration &
        visc_val > 0

      if (any(valid)) {
        iv[valid] <- unit_factor * visc_val[valid] / conc_val[valid]

        # Filter unrealistic values
        iv[valid][iv[valid] > 100] <- NA_real_ # \[eta\] > 100 dL/g is unusual
        iv[valid][iv[valid] < 0] <- NA_real_
      }

      new_measure_tbl(location = location, value = iv)
    }
  )

  new_data[[output_col]] <- new_measure_list(iv_list)

  tibble::as_tibble(new_data)
}

#' @export
print.step_sec_intrinsic_visc <- function(
  x,
  width = max(20, options()$width - 30),
  ...
) {
  title <- paste0("SEC intrinsic viscosity (", x$units, ")")
  if (x$trained) {
    cat(
      title,
      " on ",
      x$specific_visc_col,
      "/",
      x$concentration_col,
      " -> ",
      x$output_col,
      sep = ""
    )
  } else {
    cat(title)
  }
  cat("\n")
  invisible(x)
}

#' @rdname tidy.step_sec
#' @export
#' @keywords internal
tidy.step_sec_intrinsic_visc <- function(x, ...) {
  tibble::tibble(
    specific_visc_col = x$specific_visc_col %||% NA_character_,
    concentration_col = x$concentration_col %||% NA_character_,
    output_col = x$output_col,
    units = x$units,
    id = x$id
  )
}

#' @rdname required_pkgs.step_sec
#' @export
#' @keywords internal
required_pkgs.step_sec_intrinsic_visc <- function(x, ...) {
  c("measure.sec", "measure")
}
