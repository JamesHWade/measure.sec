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
#' @param angle Detection angle in degrees. Default is 90.
#' @param laser_wavelength Laser wavelength in nm.
#' @param dn_dc Refractive index increment (mL/g). Required unless
#'   `optical_constant` is provided.
#' @param solvent_ri Solvent refractive index. Default is 1.333 (water).
#' @param optical_constant Optional optical constant K; overrides dn/dc.
#' @param calibration_constant RALS instrument calibration constant. If `NULL`,
#'   results are in relative units.
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
#' RALS provides absolute MW using a 90-degree detector. It is most accurate
#' for smaller molecules where angular dependence is minimal (Rg << lambda/20).
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
  id = recipes::rand_id("sec_rals")
) {
  if (!is.null(optical_constant) && !is.null(dn_dc)) {
    cli::cli_warn(
      "{.arg optical_constant} takes precedence over {.arg dn_dc}."
    )
  }

  if (!is.numeric(angle) || angle <= 0 || angle >= 180) {
    cli::cli_abort("{.arg angle} must be between 0 and 180 degrees.")
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

  step_sec_rals_new(
    measures = measures,
    concentration_col = concentration_col,
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
      rals_val <- rals_m$value
      conc_val <- conc_m$value
      location <- rals_m$location

      max_conc <- max(abs(conc_val), na.rm = TRUE)
      threshold <- min_signal * max_conc

      mw <- rep(NA_real_, length(rals_val))
      valid <- abs(conc_val) > threshold &
        rals_val > 0 &
        !is.na(rals_val) &
        !is.na(conc_val)

      if (any(valid)) {
        R_theta <- rals_val[valid] * calibration_constant
        mw[valid] <- R_theta / (K * conc_val[valid])

        mw[valid][mw[valid] < 100] <- NA_real_
        mw[valid][mw[valid] > 1e10] <- NA_real_
      }

      new_measure_tbl(location = location, value = mw)
    }
  )

  new_data[[output_mw]] <- new_measure_list(mw_list)

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
