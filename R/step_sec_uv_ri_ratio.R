# ==============================================================================
# step_sec_uv_ri_ratio
#
# Calculate UV/RI ratio for composition analysis
# ==============================================================================

#' Calculate UV/RI Ratio for Composition Analysis
#'
#' `step_sec_uv_ri_ratio()` creates a *specification* of a recipe step that
#' calculates the ratio of UV to RI detector signals at each elution point.
#' This ratio is useful for detecting compositional heterogeneity in copolymers
#' and conjugates.
#'
#' @param recipe A recipe object.
#' @param uv_col Name of the UV detector measure column.
#' @param ri_col Name of the RI detector measure column.
#' @param input_basis Input interpretation: `"auto"` (default) selects raw
#'   detector responses when both curves are untagged and normalized
#'   concentrations when both curves carry units; `"detector_response"` and
#'   `"normalized_concentration"` select a basis explicitly. A tagged/raw mix
#'   is rejected.
#' @param uv_concentration_unit,ri_concentration_unit Units of the respective
#'   curves for `"normalized_concentration"`. Each can be omitted when the
#'   training curve records `sec_concentration_unit`; an explicit unit must
#'   agree with any recorded unit. When prep infers a unit from its training
#'   tag, the corresponding baked curve must retain that tag. An explicit unit
#'   permits untagged baked curves, but any recorded tag must agree. Both units
#'   are required for untagged normalized training curves. These arguments must
#'   be `NULL` for raw detector responses and with `input_basis = "auto"`.
#' @param output_col Name for the output ratio column. Default is `"uv_ri_ratio"`.
#' @param min_signal Minimum signal threshold (as fraction of max) below which
#'   ratio is set to `NA`. Default is 0.01 (1%). The threshold is a fraction of
#'   the maximum absolute RI value; its effective units are native RI-response
#'   units for raw mode and g/mL after conversion for normalized mode.
#' @param smooth Logical. Apply smoothing to the ratio? Default is `TRUE`.
#' @param smooth_window Window size for smoothing (number of points). Default is 5.
#' @param role Role for generated columns.
#' @param trained Logical indicating if the step has been trained.
#' @param skip Logical. Should the step be skipped when baking?
#' @param id Unique step identifier.
#'
#' @return An updated recipe with the new step added.
#'
#' @details
#' The UV/RI ratio provides a relative detector-response profile or an
#' apparent-concentration ratio across the molecular weight distribution:
#'
#' \deqn{Ratio(t) = \frac{UV(t)}{RI(t)}}
#'
#' In raw detector-response mode this is the ratio of the input detector
#' responses in their native units. In normalized-concentration mode both
#' curves are first converted to g/mL, so the ratio is dimensionless. In
#' neither mode is this ratio by itself a mass fraction; compositional
#' interpretation requires calibration of detector responses for the
#' components of interest.
#'
#' `min_signal` is a fraction of the maximum absolute RI value. The effective
#' cutoff is therefore in raw RI detector units in raw-response mode, or in
#' g/mL after unit conversion in normalized-concentration mode.
#'
#' **Applications:**
#' \itemize{
#'   \item Copolymer composition drift with molecular weight
#'   \item Block copolymer characterization
#'   \item PEGylation analysis (protein-PEG conjugates)
#'   \item Detection of chemical heterogeneity
#'   \item End-group analysis with UV labels
#' }
#'
#' **Interpretation:**
#' \itemize{
#'   \item Constant ratio: Uniform composition across MW
#'   \item Increasing ratio with MW: More chromophore in higher MW species
#'   \item Decreasing ratio with MW: Less chromophore in higher MW species
#' }
#'
#' @note
#' Both UV and RI signals should be baseline-corrected and properly aligned
#' (using `step_sec_detector_delay()`) before calculating the ratio.
#'
#' @family sec-composition
#' @export
#'
#' @examples
#' \dontrun{
#' library(recipes)
#' library(measure)
#'
#' # Calculate UV/RI ratio for copolymer
#' rec <- recipe(~., data = sec_copolymer) |>
#'   step_measure_input_long(ri_signal, location = vars(elution_time), col_name = "ri") |>
#'   step_measure_input_long(uv_signal, location = vars(elution_time), col_name = "uv") |>
#'   step_sec_detector_delay(reference = "ri", delay_volumes = c(uv = 0.05)) |>
#'   step_sec_baseline() |>
#'   step_sec_uv_ri_ratio(uv_col = "uv", ri_col = "ri") |>
#'   prep()
#' }
step_sec_uv_ri_ratio <- function(
  recipe,
  uv_col = NULL,
  ri_col = NULL,
  output_col = "uv_ri_ratio",
  min_signal = 0.01,
  smooth = TRUE,
  smooth_window = 5,
  role = NA,
  trained = FALSE,
  skip = FALSE,
  id = recipes::rand_id("sec_uv_ri_ratio"),
  input_basis = c("auto", "detector_response", "normalized_concentration"),
  uv_concentration_unit = NULL,
  ri_concentration_unit = NULL
) {
  input_basis <- match.arg(input_basis)

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

  if (!is.logical(smooth) || length(smooth) != 1 || is.na(smooth)) {
    cli::cli_abort("{.arg smooth} must be a single `TRUE` or `FALSE` value.")
  }

  if (
    !is.numeric(smooth_window) ||
      length(smooth_window) != 1 ||
      !is.finite(smooth_window) ||
      smooth_window < 1
  ) {
    cli::cli_abort("{.arg smooth_window} must be a positive integer.")
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
    input_basis %in%
      c("auto", "detector_response") &&
      (!is.null(uv_concentration_unit) || !is.null(ri_concentration_unit))
  ) {
    cli::cli_abort(
      "Concentration units require {.arg input_basis = 'normalized_concentration'}."
    )
  }

  recipes::add_step(
    recipe,
    step_sec_uv_ri_ratio_new(
      uv_col = uv_col,
      ri_col = ri_col,
      input_basis = input_basis,
      uv_concentration_unit = uv_concentration_unit,
      ri_concentration_unit = ri_concentration_unit,
      uv_concentration_unit_source = if (is.null(uv_concentration_unit)) {
        NULL
      } else {
        "explicit"
      },
      ri_concentration_unit_source = if (is.null(ri_concentration_unit)) {
        NULL
      } else {
        "explicit"
      },
      output_col = output_col,
      min_signal = min_signal,
      smooth = smooth,
      smooth_window = as.integer(smooth_window),
      role = role,
      trained = trained,
      skip = skip,
      id = id
    )
  )
}

step_sec_uv_ri_ratio_new <- function(
  uv_col,
  ri_col,
  input_basis,
  uv_concentration_unit,
  ri_concentration_unit,
  uv_concentration_unit_source,
  ri_concentration_unit_source,
  output_col,
  min_signal,
  smooth,
  smooth_window,
  role,
  trained,
  skip,
  id
) {
  recipes::step(
    subclass = "sec_uv_ri_ratio",
    uv_col = uv_col,
    ri_col = ri_col,
    input_basis = input_basis,
    uv_concentration_unit = uv_concentration_unit,
    ri_concentration_unit = ri_concentration_unit,
    uv_concentration_unit_source = uv_concentration_unit_source,
    ri_concentration_unit_source = ri_concentration_unit_source,
    output_col = output_col,
    min_signal = min_signal,
    smooth = smooth,
    smooth_window = smooth_window,
    role = role,
    trained = trained,
    skip = skip,
    id = id
  )
}

#' @export
prep.step_sec_uv_ri_ratio <- function(x, training, info = NULL, ...) {
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
  uv_concentration_unit_source <- NULL
  ri_concentration_unit_source <- NULL
  uv_recorded_unit <- attr(
    training[[uv_col]],
    "sec_concentration_unit",
    exact = TRUE
  )
  ri_recorded_unit <- attr(
    training[[ri_col]],
    "sec_concentration_unit",
    exact = TRUE
  )
  .sec_validate_concentration_unit(
    uv_recorded_unit,
    "recorded UV concentration unit"
  )
  .sec_validate_concentration_unit(
    ri_recorded_unit,
    "recorded RI concentration unit"
  )

  if (identical(input_basis, "auto")) {
    if (is.null(uv_recorded_unit) && is.null(ri_recorded_unit)) {
      input_basis <- "detector_response"
    } else if (!is.null(uv_recorded_unit) && !is.null(ri_recorded_unit)) {
      input_basis <- "normalized_concentration"
    } else {
      cli::cli_abort(
        "UV/RI ratio cannot mix a raw detector response with a concentration-tagged curve."
      )
    }
  }

  if (identical(input_basis, "detector_response")) {
    if (!is.null(uv_recorded_unit) || !is.null(ri_recorded_unit)) {
      cli::cli_abort(
        "Raw detector-response mode cannot use concentration-tagged curves."
      )
    }
  } else {
    uv_concentration_unit_source <- if (is.null(uv_concentration_unit)) {
      "inferred"
    } else {
      "explicit"
    }
    ri_concentration_unit_source <- if (is.null(ri_concentration_unit)) {
      "inferred"
    } else {
      "explicit"
    }
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

  step_sec_uv_ri_ratio_new(
    uv_col = uv_col,
    ri_col = ri_col,
    input_basis = input_basis,
    uv_concentration_unit = uv_concentration_unit,
    ri_concentration_unit = ri_concentration_unit,
    uv_concentration_unit_source = uv_concentration_unit_source,
    ri_concentration_unit_source = ri_concentration_unit_source,
    output_col = x$output_col,
    min_signal = x$min_signal,
    smooth = x$smooth,
    smooth_window = x$smooth_window,
    role = x$role,
    trained = TRUE,
    skip = x$skip,
    id = x$id
  )
}

#' Simple moving average smoother
#' @noRd
.moving_average <- function(x, window) {
  n <- length(x)
  if (n < window) {
    return(x)
  }

  result <- x
  half_window <- floor(window / 2)

  for (i in seq_len(n)) {
    start_idx <- max(1, i - half_window)
    end_idx <- min(n, i + half_window)
    # Only average non-NA values
    vals <- x[start_idx:end_idx]
    vals <- vals[!is.na(vals)]
    if (length(vals) > 0) {
      result[i] <- mean(vals)
    }
  }

  result
}

#' @export
bake.step_sec_uv_ri_ratio <- function(object, new_data, ...) {
  uv_col <- object$uv_col
  ri_col <- object$ri_col
  input_basis <- object$input_basis
  uv_concentration_unit <- object$uv_concentration_unit
  ri_concentration_unit <- object$ri_concentration_unit
  uv_concentration_unit_source <- object$uv_concentration_unit_source
  ri_concentration_unit_source <- object$ri_concentration_unit_source
  output_col <- object$output_col
  min_signal <- object$min_signal
  smooth <- object$smooth
  smooth_window <- object$smooth_window

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
      "uv_concentration_unit",
      require_recorded = identical(
        uv_concentration_unit_source,
        "inferred"
      )
    )
    .sec_resolve_concentration_unit(
      ri_concentration_unit,
      new_data[[ri_col]],
      "ri_concentration_unit",
      require_recorded = identical(
        ri_concentration_unit_source,
        "inferred"
      )
    )
  }

  # Calculate ratio for each sample
  ratio_list <- purrr::map2(
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

      # Calculate ratio where signal is above threshold
      ratio <- rep(NA_real_, length(uv_val))
      valid <- abs(ri_val) > threshold &
        is.finite(uv_val) &
        is.finite(ri_val)
      ratio[valid] <- uv_val[valid] / ri_val[valid]

      # Apply smoothing if requested
      if (smooth && sum(valid) > smooth_window) {
        ratio <- .moving_average(ratio, smooth_window)
      }
      ratio[!valid] <- NA_real_

      # Return as measure object
      new_measure_tbl(location = location, value = ratio)
    }
  )

  output <- new_measure_list(ratio_list)
  attr(output, "sec_uv_ri_ratio_basis") <- input_basis
  new_data[[output_col]] <- output

  tibble::as_tibble(new_data)
}

#' @export
print.step_sec_uv_ri_ratio <- function(
  x,
  width = max(20, options()$width - 30),
  ...
) {
  title <- "SEC UV/RI ratio"
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
tidy.step_sec_uv_ri_ratio <- function(x, ...) {
  tibble::tibble(
    uv_col = x$uv_col %||% NA_character_,
    ri_col = x$ri_col %||% NA_character_,
    input_basis = x$input_basis,
    uv_concentration_unit = x$uv_concentration_unit %||% NA_character_,
    ri_concentration_unit = x$ri_concentration_unit %||% NA_character_,
    uv_concentration_unit_source = x$uv_concentration_unit_source %||%
      NA_character_,
    ri_concentration_unit_source = x$ri_concentration_unit_source %||%
      NA_character_,
    output_col = x$output_col,
    min_signal = x$min_signal,
    smooth = x$smooth,
    id = x$id
  )
}

#' @rdname required_pkgs.step_sec
#' @export
#' @keywords internal
required_pkgs.step_sec_uv_ri_ratio <- function(x, ...) {
  c("measure.sec", "measure")
}
