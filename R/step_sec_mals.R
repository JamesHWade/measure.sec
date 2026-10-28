# ==============================================================================
# step_sec_mals
#
# Multi-Angle Light Scattering processing for absolute MW
# ==============================================================================

#' Multi-Angle Light Scattering Processing for SEC
#'
#' `step_sec_mals()` creates a *specification* of a recipe step that processes
#' one light-scattering detector signal using a single-angle approximation.
#' Multi-angle extrapolation and radius of gyration are not yet implemented.
#'
#' @param recipe A recipe object.
#' @param mals_col Name of the MALS detector measure column. For single-angle
#'   detectors (RALS/LALS), this is the only signal needed.
#' @param concentration_col Name of the concentration measure column (from
#'   `step_sec_concentration()` or similar).
#' @param concentration_unit Unit of the concentration curve: `"mg/mL"` or
#'   `"g/mL"`. When omitted, the unit must be recorded by a preceding
#'   [step_sec_concentration()] step; an untagged curve requires an explicit unit.
#' @param dn_dc Refractive index increment (mL/g). May be signed but must
#'   be finite and nonzero because the optical constant uses its square.
#'   Required for absolute MW unless `dn_dc_column` is provided.
#' @param dn_dc_column Column containing sample-specific finite, nonzero dn/dc
#'   values; either sign is valid.
#' @param wavelength Laser wavelength in nm. Default is 658 (common for MALS).
#' @param solvent_ri Solvent refractive index. Default is 1.333 (water).
#'   Common values: water = 1.333, THF = 1.407, toluene = 1.497.
#' @param angles One detection angle in degrees, default 90. The calculation
#'   assumes the particle form factor is approximately 1 at this angle.
#' @param formalism Reserved for future angular extrapolation. Only the default
#'   `"zimm"` is accepted; no angular fit is performed by this approximation.
#' @param calibration_constant Factor converting the detector signal to excess
#'   Rayleigh ratio (cm^-1). If `NULL`, results are relative and uncalibrated.
#'   The appropriate factor depends on the instrument's signal convention.
#' @param output_mw Name for the molecular weight output column. Default is `"mw_mals"`.
#' @param output_rg Must be `NULL`; radius of gyration is not yet implemented.
#' @param min_signal Minimum signal threshold (as fraction of max) below which
#'   MW is set to NA. Default is 0.01.
#' @param role Role for generated columns.
#' @param trained Logical indicating if the step has been trained.
#' @param skip Logical. Should the step be skipped when baking?
#' @param id Unique step identifier.
#'
#' @return An updated recipe with a single-angle MW estimate at each elution
#'   point. With no calibration constant the output is explicitly uncalibrated.
#'
#' @details
#' Calibrated light scattering can estimate molecular weight without polymer
#' standards. This approximation also requires an instrument calibration
#' constant and a concentration signal with known units. The relationship is:
#'
#' \deqn{\frac{K \cdot c}{R(\theta)} = \frac{1}{M_w} \cdot P(\theta) + 2A_2 c}
#'
#' where:
#' \itemize{
#'   \item K is the optical constant
#'   \item c is the concentration
#'   \item R(theta) is the excess Rayleigh ratio
#'   \item Mw is the weight-average molecular weight
#'   \item P(theta) is the particle scattering function
#'   \item A2 is the second virial coefficient
#' }
#'
#' The optical constant K is calculated as:
#' \deqn{K = \frac{4\pi^2 n_0^2 (dn/dc)^2}{N_A \lambda^4}}
#'
#' The implementation uses a single-angle approximation with P(theta) ~ 1;
#' it does not calculate an angular extrapolation or Rg.
#'
#' @note
#' Accurate results require:
#' \itemize{
#'   \item Known and accurate dn/dc value
#'   \item Calibrated instrument (calibration_constant from toluene)
#'   \item Accurate concentration from RI detector
#'   \item Clean baseline and aligned detectors
#' }
#'
#' @family sec-detectors
#' @export
#'
#' @examples
#' \dontrun{
#' library(recipes)
#' library(measure)
#'
#' # Single-angle (90 degree) light scattering
#' rec <- recipe(~., data = sec_triple_detect) |>
#'   step_measure_input_long(ri_signal, location = vars(elution_time), col_name = "ri") |>
#'   step_measure_input_long(mals_signal, location = vars(elution_time), col_name = "mals") |>
#'   step_sec_baseline() |>
#'   step_sec_ri(dn_dc = 0.185) |>
#'   step_sec_concentration(detector = "ri", injection_mass = 0.2, flow_rate = 1.0) |>
#'   step_sec_mals(
#'     mals_col = "mals",
#'     concentration_col = "ri",
#'     dn_dc = 0.185,
#'     wavelength = 658,
#'     angles = 90
#'   ) |>
#'   prep()
#' }
step_sec_mals <- function(
  recipe,
  mals_col = NULL,
  concentration_col = NULL,
  dn_dc = NULL,
  dn_dc_column = NULL,
  wavelength = 658,
  solvent_ri = 1.333,
  angles = 90,
  formalism = c("zimm", "debye", "berry"),
  calibration_constant = NULL,
  output_mw = "mw_mals",
  output_rg = NULL,
  min_signal = 0.01,
  role = NA,
  trained = FALSE,
  skip = FALSE,
  id = recipes::rand_id("sec_mals"),
  concentration_unit = NULL
) {
  formalism <- match.arg(formalism)

  # Validate inputs
  if (is.null(dn_dc) && is.null(dn_dc_column)) {
    cli::cli_abort(
      "Either {.arg dn_dc} or {.arg dn_dc_column} must be specified for absolute MW."
    )
  }

  if (!is.null(dn_dc)) {
    .sec_validate_nonzero_scalar(dn_dc, "dn_dc")
  }

  if (
    !is.numeric(wavelength) ||
      length(wavelength) != 1 ||
      !is.finite(wavelength) ||
      wavelength <= 0
  ) {
    cli::cli_abort("{.arg wavelength} must be a positive number (in nm).")
  }

  if (
    !is.numeric(solvent_ri) ||
      length(solvent_ri) != 1 ||
      !is.finite(solvent_ri) ||
      solvent_ri <= 0
  ) {
    cli::cli_abort("{.arg solvent_ri} must be a positive number.")
  }

  if (
    !is.numeric(angles) ||
      !all(is.finite(angles)) ||
      any(angles <= 0) ||
      any(angles >= 180)
  ) {
    cli::cli_abort("{.arg angles} must be between 0 and 180 degrees.")
  }
  if (length(angles) != 1) {
    cli::cli_abort(
      "Multi-angle extrapolation is not implemented; provide exactly one angle."
    )
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
  if (!identical(formalism, "zimm")) {
    cli::cli_abort(
      "Angular formalisms are not implemented for the single-angle approximation."
    )
  }
  if (!is.null(output_rg)) {
    cli::cli_abort("Radius of gyration (Rg) is not implemented by this step.")
  }
  if (
    !is.null(concentration_unit) &&
      (!is.character(concentration_unit) ||
        length(concentration_unit) != 1 ||
        !concentration_unit %in% c("mg/mL", "g/mL"))
  ) {
    cli::cli_abort(
      "{.arg concentration_unit} must be {.val mg/mL} or {.val g/mL}."
    )
  }
  if (is.null(calibration_constant)) {
    cli::cli_warn(
      "No scattering calibration constant supplied; output is relative and uncalibrated."
    )
  } else if (
    !is.numeric(calibration_constant) ||
      length(calibration_constant) != 1 ||
      !is.finite(calibration_constant) ||
      calibration_constant <= 0
  ) {
    cli::cli_abort(
      "{.arg calibration_constant} must be a positive finite number."
    )
  }

  recipes::add_step(
    recipe,
    step_sec_mals_new(
      mals_col = mals_col,
      concentration_col = concentration_col,
      concentration_unit = concentration_unit,
      dn_dc = dn_dc,
      dn_dc_column = dn_dc_column,
      wavelength = wavelength,
      solvent_ri = solvent_ri,
      angles = angles,
      formalism = formalism,
      calibration_constant = calibration_constant,
      output_mw = output_mw,
      output_rg = output_rg,
      min_signal = min_signal,
      role = role,
      trained = trained,
      skip = skip,
      id = id
    )
  )
}

step_sec_mals_new <- function(
  mals_col,
  concentration_col,
  concentration_unit,
  dn_dc,
  dn_dc_column,
  wavelength,
  solvent_ri,
  angles,
  formalism,
  calibration_constant,
  output_mw,
  output_rg,
  min_signal,
  role,
  trained,
  skip,
  id
) {
  recipes::step(
    subclass = "sec_mals",
    mals_col = mals_col,
    concentration_col = concentration_col,
    concentration_unit = concentration_unit,
    dn_dc = dn_dc,
    dn_dc_column = dn_dc_column,
    wavelength = wavelength,
    solvent_ri = solvent_ri,
    angles = angles,
    formalism = formalism,
    calibration_constant = calibration_constant,
    output_mw = output_mw,
    output_rg = output_rg,
    min_signal = min_signal,
    role = role,
    trained = trained,
    skip = skip,
    id = id
  )
}

#' Calculate optical constant K
#' @noRd
.optical_constant <- function(dn_dc, solvent_ri, wavelength_nm) {
  .sec_validate_nonzero_scalar(dn_dc, "dn_dc")
  .sec_validate_positive_scalar(solvent_ri, "solvent_ri")
  .sec_validate_positive_scalar(wavelength_nm, "wavelength")

  # K = 4 * pi^2 * n0^2 * (dn/dc)^2 / (NA * lambda^4)
  # wavelength in nm -> cm
  # dn/dc in mL/g = cm^3/g
  # K has units mol*cm^2/g^2 when lambda is in cm.

  lambda_cm <- wavelength_nm * 1e-7 # nm to cm
  NA_val <- 6.022e23 # Avogadro's number

  K <- 4 * pi^2 * solvent_ri^2 * dn_dc^2 / (NA_val * lambda_cm^4)
  K
}

#' @export
prep.step_sec_mals <- function(x, training, info = NULL, ...) {
  check_for_measure(training)

  measure_cols <- find_measure_cols(training)

  # Find MALS column if not specified
  if (is.null(x$mals_col)) {
    mals_cols <- measure_cols[
      grepl(
        "mals|ls|light",
        measure_cols,
        ignore.case = TRUE
      )
    ]
    if (length(mals_cols) == 0) {
      cli::cli_abort(
        "No MALS column found. Specify {.arg mals_col} explicitly."
      )
    }
    mals_col <- mals_cols[1]
  } else {
    mals_col <- x$mals_col
    if (!mals_col %in% measure_cols) {
      cli::cli_abort(
        "MALS column {.val {mals_col}} not found in measure columns."
      )
    }
  }

  # Find concentration column if not specified
  if (is.null(x$concentration_col)) {
    # Try to find RI or concentration column
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

  # Validate dn_dc_column exists if specified
  if (!is.null(x$dn_dc_column)) {
    if (!x$dn_dc_column %in% names(training)) {
      cli::cli_abort(
        "Column {.val {x$dn_dc_column}} not found in training data."
      )
    }
  }

  step_sec_mals_new(
    mals_col = mals_col,
    concentration_col = concentration_col,
    concentration_unit = concentration_unit,
    dn_dc = x$dn_dc,
    dn_dc_column = x$dn_dc_column,
    wavelength = x$wavelength,
    solvent_ri = x$solvent_ri,
    angles = x$angles,
    formalism = x$formalism,
    calibration_constant = x$calibration_constant,
    output_mw = x$output_mw,
    output_rg = x$output_rg,
    min_signal = x$min_signal,
    role = x$role,
    trained = TRUE,
    skip = x$skip,
    id = x$id
  )
}

#' @export
bake.step_sec_mals <- function(object, new_data, ...) {
  mals_col <- object$mals_col
  concentration_col <- object$concentration_col
  concentration_unit <- object$concentration_unit
  .sec_resolve_concentration_unit(
    concentration_unit,
    new_data[[concentration_col]]
  )
  dn_dc <- object$dn_dc
  dn_dc_column <- object$dn_dc_column
  wavelength <- object$wavelength
  solvent_ri <- object$solvent_ri
  angles <- object$angles
  calibration_constant <- object$calibration_constant
  output_mw <- object$output_mw
  min_signal <- object$min_signal

  # Get dn/dc values
  if (!is.null(dn_dc_column)) {
    dn_dc_values <- new_data[[dn_dc_column]]
  } else {
    dn_dc_values <- rep(dn_dc, nrow(new_data))
  }
  if (
    !is.numeric(dn_dc_values) ||
      length(dn_dc_values) != nrow(new_data) ||
      !all(is.finite(dn_dc_values)) ||
      any(dn_dc_values == 0)
  ) {
    cli::cli_abort("Each dn/dc value must be nonzero and finite.")
  }

  # Default calibration constant if not provided
  if (is.null(calibration_constant)) {
    calibration_constant <- 1.0
  }

  # Calculate MW for each sample
  mw_list <- purrr::pmap(
    list(
      new_data[[mals_col]],
      new_data[[concentration_col]],
      dn_dc_values
    ),
    function(mals_m, conc_m, dndc) {
      mals_val <- mals_m$value
      conc_val <- conc_m$value
      .sec_assert_aligned_curves(mals_m, conc_m, "Scattering")
      location <- mals_m$location
      if (length(dndc) != 1 || !is.finite(dndc) || dndc == 0) {
        cli::cli_abort("Each dn/dc value must be nonzero and finite.")
      }

      # Calculate optical constant
      K <- .optical_constant(dndc, solvent_ri, wavelength)

      # Determine signal threshold
      max_conc <- max(conc_val[is.finite(conc_val)], 0)
      threshold <- min_signal * max_conc

      # Calculate MW at each point
      mw <- rep(NA_real_, length(mals_val))
      valid <- conc_val > threshold &
        is.finite(mals_val) &
        is.finite(conc_val) &
        mals_val > 0

      if (any(valid)) {
        # For single angle at 90 degrees, P(theta) ~ 1 for small particles
        # Mw = R(theta) / (K * c)
        # With calibration: R = signal * calibration_constant
        R_theta <- mals_val[valid] * calibration_constant
        conc_g_per_ml <- .sec_concentration_g_per_ml(
          conc_val[valid],
          concentration_unit
        )
        mw[valid] <- R_theta / (K * conc_g_per_ml)
      }

      # Return as measure object
      new_measure_tbl(location = location, value = mw)
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
print.step_sec_mals <- function(
  x,
  width = max(20, options()$width - 30),
  ...
) {
  title <- paste0(
    "SEC MALS processing (",
    length(x$angles),
    " angle",
    if (length(x$angles) > 1) "s" else "",
    ")"
  )
  if (x$trained) {
    cat(title, " on ", x$mals_col, " -> ", x$output_mw, sep = "")
  } else {
    cat(title)
  }
  cat("\n")
  invisible(x)
}

#' @rdname tidy.step_sec
#' @export
#' @keywords internal
tidy.step_sec_mals <- function(x, ...) {
  tibble::tibble(
    mals_col = x$mals_col %||% NA_character_,
    concentration_col = x$concentration_col %||% NA_character_,
    concentration_unit = x$concentration_unit %||% NA_character_,
    dn_dc = x$dn_dc %||% NA_real_,
    wavelength = x$wavelength,
    angles = list(x$angles),
    formalism = x$formalism,
    output_mw = x$output_mw,
    id = x$id
  )
}

#' @rdname required_pkgs.step_sec
#' @export
#' @keywords internal
required_pkgs.step_sec_mals <- function(x, ...) {
  c("measure.sec", "measure")
}
