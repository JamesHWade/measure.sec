# ==============================================================================
# step_sec_rals
#
# Right-angle light scattering processing
# ==============================================================================

#' Right-Angle Light Scattering Processing for SEC
#'
#' `step_sec_rals()` creates a *specification* of a recipe step that processes
#' right-angle light scattering (RALS) signals for absolute molecular weight.
#'
#' @param recipe A recipe object.
#' @param measures Character vector of RALS measure columns. If `NULL`, the
#'   step searches for measure columns containing "rals".
#' @param concentration_col Name of the concentration measure column (from
#'   `step_sec_concentration()` or similar).
#' @param concentration_unit Unit of the concentration curve: "mg/mL" or
#'   "g/mL". When omitted, the training curve must record its unit. If prep
#'   infers the unit from that tag, baked curves must retain it. An explicit
#'   unit permits untagged baked curves, but any recorded tag must agree.
#' @param angle Detection angle in degrees. Default is 90.
#' @param laser_wavelength Laser wavelength in nm.
#' @param dn_dc Refractive index increment (mL/g). May be signed but must
#'   be finite and nonzero because the optical constant uses its square.
#'   Required unless `optical_constant` is provided.
#' @param solvent_ri Solvent refractive index. Default is 1.333 (water).
#' @param optical_constant Optional optical constant K; overrides dn/dc. K
#'   includes (dn/dc)^2 and uses units mol*cm^2/g^2.
#' @param calibration_constant RALS instrument calibration constant. If
#'   `NULL`, results are in relative units. Otherwise it converts the detector
#'   signal to excess Rayleigh ratio in cm^-1.
#' @param output_mw Name for the molecular weight output column.
#' @param min_signal Minimum signal threshold (as fraction of max) below which
#'   MW is set to NA.
#' @param role Role for generated columns.
#' @param trained Logical indicating if the step has been trained.
#' @param skip Logical. Should the step be skipped when baking?
#' @param id Unique step identifier.
#'
#' @return An updated recipe with the new step added.
#'
#' @details
#' RALS provides an estimate from a 90-degree detector. For a calibrated
#' signal, the estimate is
#' \deqn{M = R(\theta)/(Kc)}
#' K = 4 pi^2 n_0^2 (dn/dc)^2/(N_A lambda_0^4), c is converted to g/mL, and
#' lambda_0 is in cm when calculating K. This is a single-angle estimate; it
#' does not fit multiple angles, extrapolate to zero angle, or estimate Rg.
#' The particle form factor is assumed to be approximately 1, so angular
#' dependence must be negligible. Without a detector calibration constant,
#' the result is relative and must not be interpreted as a molar mass.
#'
#' @family sec-detectors
#' @export
#'
#' @examples
#' \dontrun{
#' library(recipes)
#' library(measure)
#'
#' rec <- recipe(~., data = sec_triple_detect) |>
#'   step_measure_input_long(ri_signal, location = vars(elution_time), col_name = "ri") |>
#'   step_measure_input_long(rals_signal, location = vars(elution_time), col_name = "rals") |>
#'   step_sec_concentration(detector = "ri", injection_mass = 0.2, flow_rate = 1.0) |>
#'   step_sec_rals(measures = "rals", concentration_col = "ri", dn_dc = 0.185) |>
#'   prep()
#' }
step_sec_rals <- function(
  recipe,
  measures = NULL,
  concentration_col = NULL,
  angle = 90,
  laser_wavelength = 670,
  dn_dc = NULL,
  solvent_ri = 1.333,
  optical_constant = NULL,
  calibration_constant = NULL,
  output_mw = "mw_rals",
  min_signal = 0.01,
  role = NA,
  trained = FALSE,
  skip = FALSE,
  id = recipes::rand_id("sec_rals"),
  concentration_unit = NULL
) {
  if (!is.null(optical_constant) && !is.null(dn_dc)) {
    cli::cli_warn(
      "{.arg optical_constant} takes precedence over {.arg dn_dc}."
    )
  }

  if (
    !is.numeric(angle) ||
      length(angle) != 1 ||
      !is.finite(angle) ||
      angle <= 0 ||
      angle >= 180
  ) {
    cli::cli_abort("{.arg angle} must be between 0 and 180 degrees.")
  }
  .sec_validate_positive_scalar(laser_wavelength, "laser_wavelength")
  .sec_validate_positive_scalar(solvent_ri, "solvent_ri")
  .sec_validate_concentration_unit(concentration_unit)
  if (!is.null(dn_dc)) {
    .sec_validate_nonzero_scalar(dn_dc, "dn_dc")
  }
  if (!is.null(optical_constant)) {
    .sec_validate_positive_scalar(optical_constant, "optical_constant")
  }
  if (!is.null(calibration_constant)) {
    .sec_validate_positive_scalar(calibration_constant, "calibration_constant")
  }
  if (
    !is.numeric(min_signal) ||
      length(min_signal) != 1 ||
      !is.finite(min_signal) ||
      min_signal < 0 ||
      min_signal >= 1
  ) {
    cli::cli_abort(
      "{.arg min_signal} must be a finite fraction from 0 up to but not including 1."
    )
  }

  if (is.null(optical_constant) && is.null(dn_dc)) {
    cli::cli_abort(
      "Either {.arg dn_dc} or {.arg optical_constant} must be provided."
    )
  }

  recipes::add_step(
    recipe,
    step_sec_rals_new(
      measures = measures,
      concentration_col = concentration_col,
      concentration_unit = concentration_unit,
      concentration_unit_source = if (is.null(concentration_unit)) {
        NULL
      } else {
        "explicit"
      },
      angle = angle,
      laser_wavelength = laser_wavelength,
      dn_dc = dn_dc,
      solvent_ri = solvent_ri,
      optical_constant = optical_constant,
      calibration_constant = calibration_constant,
      output_mw = output_mw,
      min_signal = min_signal,
      role = role,
      trained = trained,
      skip = skip,
      id = id
    )
  )
}

step_sec_rals_new <- function(
  measures,
  concentration_col,
  concentration_unit,
  concentration_unit_source,
  angle,
  laser_wavelength,
  dn_dc,
  solvent_ri,
  optical_constant,
  calibration_constant,
  output_mw,
  min_signal,
  role,
  trained,
  skip,
  id
) {
  recipes::step(
    subclass = "sec_rals",
    measures = measures,
    concentration_col = concentration_col,
    concentration_unit = concentration_unit,
    concentration_unit_source = concentration_unit_source,
    angle = angle,
    laser_wavelength = laser_wavelength,
    dn_dc = dn_dc,
    solvent_ri = solvent_ri,
    optical_constant = optical_constant,
    calibration_constant = calibration_constant,
    output_mw = output_mw,
    min_signal = min_signal,
    role = role,
    trained = trained,
    skip = skip,
    id = id
  )
}

#' @export
prep.step_sec_rals <- function(x, training, info = NULL, ...) {
  check_for_measure(training)

  measure_cols <- find_measure_cols(training)

  if (is.null(x$measures)) {
    rals_cols <- measure_cols[grepl("rals", measure_cols, ignore.case = TRUE)]
    if (length(rals_cols) == 0) {
      cli::cli_abort(
        "No RALS column found. Specify {.arg measures} explicitly."
      )
    }
    measures <- rals_cols[1]
  } else {
    measures <- x$measures
  }

  if (length(measures) != 1) {
    cli::cli_abort("{.arg measures} must specify a single RALS column.")
  }

  if (!measures %in% measure_cols) {
    cli::cli_abort(
      "RALS column {.val {measures}} not found in measure columns."
    )
  }

  if (is.null(x$concentration_col)) {
    conc_cols <- measure_cols[
      grepl(
        "ri|conc",
        measure_cols,
        ignore.case = TRUE
      )
    ]
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

  concentration_unit <- .sec_resolve_concentration_unit(
    x$concentration_unit,
    training[[concentration_col]]
  )
  concentration_unit_source <- if (is.null(x$concentration_unit)) {
    "inferred"
  } else {
    "explicit"
  }

  step_sec_rals_new(
    measures = measures,
    concentration_col = concentration_col,
    concentration_unit = concentration_unit,
    concentration_unit_source = concentration_unit_source,
    angle = x$angle,
    laser_wavelength = x$laser_wavelength,
    dn_dc = x$dn_dc,
    solvent_ri = x$solvent_ri,
    optical_constant = x$optical_constant,
    calibration_constant = x$calibration_constant,
    output_mw = x$output_mw,
    min_signal = x$min_signal,
    role = x$role,
    trained = TRUE,
    skip = x$skip,
    id = x$id
  )
}

#' @export
bake.step_sec_rals <- function(object, new_data, ...) {
  rals_col <- object$measures
  concentration_col <- object$concentration_col
  concentration_unit <- object$concentration_unit
  concentration_unit_source <- object$concentration_unit_source
  .sec_resolve_concentration_unit(
    concentration_unit,
    new_data[[concentration_col]],
    require_recorded = identical(concentration_unit_source, "inferred")
  )
  dn_dc <- object$dn_dc
  solvent_ri <- object$solvent_ri
  laser_wavelength <- object$laser_wavelength
  calibration_constant <- object$calibration_constant
  min_signal <- object$min_signal
  output_mw <- object$output_mw

  if (is.null(calibration_constant)) {
    cli::cli_warn(
      c(
        "No {.arg calibration_constant} provided.",
        "i" = "Results will be in relative units, not absolute MW.",
        "i" = "Provide a calibration constant for absolute values."
      )
    )
    calibration_constant <- 1.0
  }

  if (is.null(object$optical_constant)) {
    K <- .optical_constant(dn_dc, solvent_ri, laser_wavelength)
  } else {
    K <- object$optical_constant
  }

  mw_list <- purrr::pmap(
    list(new_data[[rals_col]], new_data[[concentration_col]]),
    function(rals_m, conc_m) {
      .sec_single_angle_mw(
        rals_m,
        conc_m,
        concentration_unit,
        K,
        calibration_constant,
        min_signal,
        "RALS"
      )
    }
  )

  output <- new_measure_list(mw_list)
  attr(output, "sec_scattering_basis") <- if (
    is.null(object$calibration_constant)
  ) {
    "relative_uncalibrated"
  } else {
    "calibrated_single_angle"
  }
  attr(output, "sec_scattering_source") <- concentration_col
  if (!is.null(object$calibration_constant)) {
    output <- .sec_tag_mw_curve(
      output,
      source = concentration_col,
      scale = "absolute"
    )
  }
  new_data[[output_mw]] <- output

  tibble::as_tibble(new_data)
}

#' @export
print.step_sec_rals <- function(
  x,
  width = max(20, options()$width - 30),
  ...
) {
  title <- "SEC RALS processing"
  if (x$trained) {
    cat(title, " on ", x$measures, " -> ", x$output_mw, sep = "")
  } else {
    cat(title)
  }
  cat("\n")
  invisible(x)
}

#' @rdname tidy.step_sec
#' @export
#' @keywords internal
tidy.step_sec_rals <- function(x, ...) {
  tibble::tibble(
    rals_col = x$measures %||% NA_character_,
    concentration_col = x$concentration_col %||% NA_character_,
    concentration_unit = x$concentration_unit %||% NA_character_,
    concentration_unit_source = x$concentration_unit_source %||% NA_character_,
    angle = x$angle,
    laser_wavelength = x$laser_wavelength,
    dn_dc = x$dn_dc %||% NA_real_,
    output_mw = x$output_mw,
    id = x$id
  )
}

#' @rdname required_pkgs.step_sec
#' @export
#' @keywords internal
required_pkgs.step_sec_rals <- function(x, ...) {
  c("measure.sec", "measure")
}
