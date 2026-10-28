# ==============================================================================
# step_sec_composition
#
# Calculate copolymer composition from UV/RI signals
# ==============================================================================

#' Calculate Copolymer Composition from Detector Signals
#'
#' `step_sec_composition()` creates a *specification* of a recipe step that
#' calculates the weight fraction of components in a copolymer or blend using
#' UV and RI detector signals with known response factors.
#'
#' @param recipe A recipe object.
#' @param uv_col Name of the UV detector measure column.
#' @param ri_col Name of the RI detector measure column.
#' @param component_a_uv,component_b_uv UV response coefficients for components
#'   A and B, calibrated in the selected input-channel basis per g/mL.
#' @param component_a_ri,component_b_ri RI response coefficients for components
#'   A and B, calibrated in the selected input-channel basis per g/mL.
#' @param input_basis Input interpretation: `"detector_response"` (default)
#'   uses untagged raw UV and RI detector responses; `"normalized_concentration"`
#'   uses apparent concentration channels and converts them to g/mL before
#'   calculating composition.
#' @param uv_concentration_unit,ri_concentration_unit Units of the respective
#'   curves when `input_basis = "normalized_concentration"`. Each can be
#'   omitted when that curve records `sec_concentration_unit`; an explicit
#'   unit must agree with any recorded unit. Both are required for untagged
#'   curves. These arguments must be `NULL` for raw detector responses.
#' @param output_col Name for the output composition column. Default is
#'   `"composition_a"`. Contains weight fraction of component A (0-1).
#' @param min_signal Minimum signal threshold (as fraction of max) below which
#'   composition is set to `NA`. Default is 0.01. The threshold is a fraction of
#'   the maximum absolute RI value; its effective units are native RI-response
#'   units for raw mode and g/mL after conversion for normalized mode.
#' @param clip Logical. Clip composition values to \[0, 1\] range? Default is `TRUE`.
#' @param role Role for generated columns.
#' @param trained Logical indicating if the step has been trained.
#' @param skip Logical. Should the step be skipped when baking?
#' @param id Unique step identifier.
#'
#' @return An updated recipe with the new step added.
#'
#' @details
#' For a two-component system, the detector signals are:
#'
#' \deqn{UV = u_A c_A + u_B c_B, \qquad RI = r_A c_A + r_B c_B}
#'
#' where `c_A` and `c_B` are the component mass concentrations in g/mL. The
#' coefficients `u_A`, `u_B`, `r_A`, and `r_B` are the response coefficients
#' for the selected `input_basis`, measured in the same input-channel basis
#' per g/mL. In raw detector-response mode these are calibrated detector
#' responses per g/mL; in normalized-concentration mode they are calibrated
#' apparent channel responses per g/mL. They are not inferred from physical
#' extinction coefficients or `dn/dc` values.
#'
#' Composition is calculated by solving the two-response linear system and
#' converting the recovered component concentrations to a mass fraction:
#'
#' \deqn{c_A = \frac{UV r_B - u_B RI}{u_A r_B - u_B r_A}, \qquad
#' c_B = \frac{u_A RI - UV r_A}{u_A r_B - u_B r_A}, \qquad
#' w_A = \frac{c_A}{c_A + c_B}}
#'
#' The response matrix must be identifiable: the component response vectors
#' cannot be proportional. A zero coefficient is allowed when the full matrix
#' remains identifiable.
#'
#' `min_signal` is a fraction of the maximum absolute RI value. The effective
#' cutoff is therefore in raw RI detector units for `"detector_response"`, or
#' in g/mL after unit conversion for `"normalized_concentration"`.
#'
#' **Common Applications:**
#' \itemize{
#'   \item Styrene-acrylate copolymers (styrene is UV-active)
#'   \item Block copolymers with different chromophore content
#'   \item PEGylated proteins (protein at 280nm, PEG is UV-transparent)
#'   \item Polymer blends with known compositions
#' }
#'
#' Calibrate each response coefficient from the measured channel output of a
#' pure-component standard at known concentrations in g/mL. A literature
#' extinction coefficient or `dn/dc` value alone is not a substitute; account
#' for path length, instrument calibration, and signal scaling when relating it
#' to the actual input-channel values.
#'
#' @note
#' The response coefficients must be calibrated for the selected input basis
#' and the same g/mL concentration denominator. Raw detector responses and
#' normalized concentration channels are separate models; switching between
#' them requires recalibrating all four coefficients. If both channels have
#' been normalized to the same mass concentration response for both
#' components, the response matrix is singular and composition is not
#' identifiable from those channels.
#'
#' @family sec-composition
#' @export
#'
#' @examples
#' \dontrun{
#' library(recipes)
#' library(measure)
#'
#' # Calibrate the response coefficients with pure-component standards.
#' rec <- recipe(~., data = copolymer_data) |>
#'   step_measure_input_long(ri_signal, location = vars(elution_time), col_name = "ri") |>
#'   step_measure_input_long(uv_signal, location = vars(elution_time), col_name = "uv") |>
#'   step_sec_baseline() |>
#'   step_sec_composition(
#'     uv_col = "uv",
#'     ri_col = "ri",
#'     component_a_uv = 2.0,    # measured UV response per g/mL
#'     component_a_ri = 1.0,    # measured RI response per g/mL
#'     component_b_uv = 0.5,
#'     component_b_ri = 3.0,
#'     output_col = "styrene_fraction"
#'   ) |>
#'   prep()
#' }
step_sec_composition <- function(
  recipe,
  uv_col = NULL,
  ri_col = NULL,
  component_a_uv,
  component_a_ri,
  component_b_uv,
  component_b_ri,
  output_col = "composition_a",
  min_signal = 0.01,
  clip = TRUE,
  role = NA,
  trained = FALSE,
  skip = FALSE,
  id = recipes::rand_id("sec_composition"),
  input_basis = c("detector_response", "normalized_concentration"),
  uv_concentration_unit = NULL,
  ri_concentration_unit = NULL
) {
  input_basis <- match.arg(input_basis)

  # Validate response factors
  if (
    missing(component_a_uv) ||
      missing(component_a_ri) ||
      missing(component_b_uv) ||
      missing(component_b_ri)
  ) {
    cli::cli_abort(
      "All response factors must be specified: {.arg component_a_uv}, {.arg component_a_ri}, {.arg component_b_uv}, {.arg component_b_ri}."
    )
  }

  response_factors <- list(
    component_a_uv = component_a_uv,
    component_a_ri = component_a_ri,
    component_b_uv = component_b_uv,
    component_b_ri = component_b_ri
  )
  for (argument in names(response_factors)) {
    value <- response_factors[[argument]]
    if (
      !is.numeric(value) ||
        length(value) != 1 ||
        !is.finite(value)
    ) {
      cli::cli_abort(
        "{.arg {argument}} must be a finite numeric response coefficient."
      )
    }
  }

  .sec_validate_concentration_unit(
    uv_concentration_unit,
    "uv_concentration_unit"
  )
  .sec_validate_concentration_unit(
    ri_concentration_unit,
    "ri_concentration_unit"
  )

  if (
    identical(input_basis, "detector_response") &&
      (!is.null(uv_concentration_unit) || !is.null(ri_concentration_unit))
  ) {
    cli::cli_abort(
      c(
        "Concentration units cannot be supplied for raw detector responses.",
        "i" = "Use {.arg input_basis = 'normalized_concentration'} for concentration channels."
      )
    )
  }

  determinant <- component_a_uv *
    component_b_ri -
    component_b_uv * component_a_ri
  determinant_scale <- max(
    abs(component_a_uv * component_b_ri),
    abs(component_b_uv * component_a_ri)
  )
  if (
    !is.finite(determinant) ||
      !is.finite(determinant_scale) ||
      determinant_scale == 0 ||
      abs(determinant) <= sqrt(.Machine$double.eps) * determinant_scale
  ) {
    cli::cli_abort(
      c(
        "The component response matrix cannot identify composition.",
        "i" = "The UV and RI response vectors must not be proportional."
      )
    )
  }

  if (
    !is.numeric(min_signal) ||
      length(min_signal) != 1 ||
      !is.finite(min_signal) ||
      min_signal < 0 ||
      min_signal > 1
  ) {
    cli::cli_abort(
      "{.arg min_signal} must be a finite fraction between 0 and 1."
    )
  }
  if (!is.logical(clip) || length(clip) != 1 || is.na(clip)) {
    cli::cli_abort("{.arg clip} must be a single `TRUE` or `FALSE` value.")
  }

  ratio_a <- if (component_a_ri == 0) {
    NA_real_
  } else {
    component_a_uv / component_a_ri
  }
  ratio_b <- if (component_b_ri == 0) {
    NA_real_
  } else {
    component_b_uv / component_b_ri
  }

  recipes::add_step(
    recipe,
    step_sec_composition_new(
      uv_col = uv_col,
      ri_col = ri_col,
      component_a_uv = component_a_uv,
      component_a_ri = component_a_ri,
      component_b_uv = component_b_uv,
      component_b_ri = component_b_ri,
      ratio_a = ratio_a,
      ratio_b = ratio_b,
      determinant = determinant,
      output_col = output_col,
      min_signal = min_signal,
      clip = clip,
      input_basis = input_basis,
      uv_concentration_unit = uv_concentration_unit,
      ri_concentration_unit = ri_concentration_unit,
      role = role,
      trained = trained,
      skip = skip,
      id = id
    )
  )
}

step_sec_composition_new <- function(
  uv_col,
  ri_col,
  component_a_uv,
  component_a_ri,
  component_b_uv,
  component_b_ri,
  ratio_a,
  ratio_b,
  determinant,
  output_col,
  min_signal,
  clip,
  input_basis,
  uv_concentration_unit,
  ri_concentration_unit,
  role,
  trained,
  skip,
  id
) {
  recipes::step(
    subclass = "sec_composition",
    uv_col = uv_col,
    ri_col = ri_col,
    component_a_uv = component_a_uv,
    component_a_ri = component_a_ri,
    component_b_uv = component_b_uv,
    component_b_ri = component_b_ri,
    ratio_a = ratio_a,
    ratio_b = ratio_b,
    determinant = determinant,
    output_col = output_col,
    min_signal = min_signal,
    clip = clip,
    input_basis = input_basis,
    uv_concentration_unit = uv_concentration_unit,
    ri_concentration_unit = ri_concentration_unit,
    role = role,
    trained = trained,
    skip = skip,
    id = id
  )
}

#' @export
prep.step_sec_composition <- function(x, training, info = NULL, ...) {
  check_for_measure(training)

  measure_cols <- find_measure_cols(training)

  # Find UV column if not specified
  if (is.null(x$uv_col)) {
    uv_cols <- measure_cols[grepl("uv", measure_cols, ignore.case = TRUE)]
    if (length(uv_cols) == 0) {
      cli::cli_abort("No UV column found. Specify {.arg uv_col} explicitly.")
    }
    uv_col <- uv_cols[1]
  } else {
    uv_col <- x$uv_col
    if (!uv_col %in% measure_cols) {
      cli::cli_abort("UV column {.val {uv_col}} not found in measure columns.")
    }
  }

  # Find RI column if not specified
  if (is.null(x$ri_col)) {
    ri_cols <- measure_cols[grepl("ri", measure_cols, ignore.case = TRUE)]
    if (length(ri_cols) == 0) {
      cli::cli_abort("No RI column found. Specify {.arg ri_col} explicitly.")
    }
    ri_col <- ri_cols[1]
  } else {
    ri_col <- x$ri_col
    if (!ri_col %in% measure_cols) {
      cli::cli_abort("RI column {.val {ri_col}} not found in measure columns.")
    }
  }

  input_basis <- x$input_basis
  uv_concentration_unit <- x$uv_concentration_unit
  ri_concentration_unit <- x$ri_concentration_unit
  if (identical(input_basis, "detector_response")) {
    if (
      !is.null(attr(
        training[[uv_col]],
        "sec_concentration_unit",
        exact = TRUE
      )) ||
        !is.null(attr(
          training[[ri_col]],
          "sec_concentration_unit",
          exact = TRUE
        ))
    ) {
      cli::cli_abort(
        c(
          "Raw detector-response mode cannot use concentration-tagged curves.",
          "i" = "Use {.arg input_basis = 'normalized_concentration'} for normalized channels."
        )
      )
    }
  } else {
    uv_concentration_unit <- .sec_resolve_concentration_unit(
      uv_concentration_unit,
      training[[uv_col]],
      "uv_concentration_unit"
    )
    ri_concentration_unit <- .sec_resolve_concentration_unit(
      ri_concentration_unit,
      training[[ri_col]],
      "ri_concentration_unit"
    )
  }

  step_sec_composition_new(
    uv_col = uv_col,
    ri_col = ri_col,
    component_a_uv = x$component_a_uv,
    component_a_ri = x$component_a_ri,
    component_b_uv = x$component_b_uv,
    component_b_ri = x$component_b_ri,
    ratio_a = x$ratio_a,
    ratio_b = x$ratio_b,
    determinant = x$determinant,
    output_col = x$output_col,
    min_signal = x$min_signal,
    clip = x$clip,
    input_basis = input_basis,
    uv_concentration_unit = uv_concentration_unit,
    ri_concentration_unit = ri_concentration_unit,
    role = x$role,
    trained = TRUE,
    skip = x$skip,
    id = x$id
  )
}

#' @export
bake.step_sec_composition <- function(object, new_data, ...) {
  uv_col <- object$uv_col
  ri_col <- object$ri_col
  ratio_a <- object$ratio_a
  ratio_b <- object$ratio_b
  output_col <- object$output_col
  min_signal <- object$min_signal
  clip <- object$clip
  input_basis <- object$input_basis
  uv_concentration_unit <- object$uv_concentration_unit
  ri_concentration_unit <- object$ri_concentration_unit

  if (identical(input_basis, "detector_response")) {
    if (
      !is.null(attr(
        new_data[[uv_col]],
        "sec_concentration_unit",
        exact = TRUE
      )) ||
        !is.null(attr(
          new_data[[ri_col]],
          "sec_concentration_unit",
          exact = TRUE
        ))
    ) {
      cli::cli_abort(
        "Raw detector-response mode cannot bake concentration-tagged curves."
      )
    }
  } else {
    .sec_resolve_concentration_unit(
      uv_concentration_unit,
      new_data[[uv_col]],
      "uv_concentration_unit"
    )
    .sec_resolve_concentration_unit(
      ri_concentration_unit,
      new_data[[ri_col]],
      "ri_concentration_unit"
    )
  }

  # Calculate composition for each sample
  comp_list <- purrr::map2(
    new_data[[uv_col]],
    new_data[[ri_col]],
    function(uv_m, ri_m) {
      uv_val <- uv_m$value
      ri_val <- ri_m$value
      location <- uv_m$location

      .sec_assert_aligned_curves(uv_m, ri_m, "UV and RI")
      if (identical(input_basis, "normalized_concentration")) {
        uv_val <- .sec_concentration_g_per_ml(
          uv_val,
          uv_concentration_unit
        )
        ri_val <- .sec_concentration_g_per_ml(
          ri_val,
          ri_concentration_unit
        )
      }

      # Determine signal threshold
      max_ri <- max(abs(ri_val[is.finite(ri_val)]), 0)
      threshold <- min_signal * max_ri

      # Solve the response matrix and calculate the component mass fraction
      composition <- rep(NA_real_, length(uv_val))
      valid <- abs(ri_val) > threshold &
        is.finite(uv_val) &
        is.finite(ri_val)

      if (any(valid)) {
        component_a_conc <- (uv_val[valid] *
          object$component_b_ri -
          object$component_b_uv * ri_val[valid]) /
          object$determinant
        component_b_conc <- (object$component_a_uv *
          ri_val[valid] -
          uv_val[valid] * object$component_a_ri) /
          object$determinant
        total_conc <- component_a_conc + component_b_conc
        valid_total <- is.finite(total_conc) & total_conc > 0
        valid_indices <- which(valid)
        composition[valid_indices[valid_total]] <-
          component_a_conc[valid_total] / total_conc[valid_total]

        # Clip to \[0, 1\] if requested
        if (clip) {
          composition[valid_indices[valid_total]] <- pmax(
            0,
            pmin(1, composition[valid_indices[valid_total]])
          )
        }
      }

      # Return as measure object
      new_measure_tbl(location = location, value = composition)
    }
  )

  new_data[[output_col]] <- new_measure_list(comp_list)

  tibble::as_tibble(new_data)
}

#' @export
print.step_sec_composition <- function(
  x,
  width = max(20, options()$width - 30),
  ...
) {
  title <- paste0("SEC composition analysis (", x$input_basis, ")")
  if (x$trained) {
    cat(
      title,
      " (",
      x$uv_col,
      "/",
      x$ri_col,
      " -> ",
      x$output_col,
      ")",
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
tidy.step_sec_composition <- function(x, ...) {
  tibble::tibble(
    uv_col = x$uv_col %||% NA_character_,
    ri_col = x$ri_col %||% NA_character_,
    output_col = x$output_col,
    ratio_a = x$ratio_a,
    ratio_b = x$ratio_b,
    determinant = x$determinant,
    input_basis = x$input_basis,
    uv_concentration_unit = x$uv_concentration_unit %||% NA_character_,
    ri_concentration_unit = x$ri_concentration_unit %||% NA_character_,
    min_signal = x$min_signal,
    id = x$id
  )
}

#' @rdname required_pkgs.step_sec
#' @export
#' @keywords internal
required_pkgs.step_sec_composition <- function(x, ...) {
  c("measure.sec", "measure")
}
