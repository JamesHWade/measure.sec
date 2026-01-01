# ==============================================================================
# step_sec_universal_cal
#
# Universal calibration for SEC using hydrodynamic volume
# ==============================================================================

#' Universal Calibration for SEC
#'
#' `step_sec_universal_cal()` creates a *specification* of a recipe step that
#' applies universal calibration to determine molecular weight from intrinsic
#' viscosity and retention data.
#'
#' @param recipe A recipe object.
#' @param measures Character vector of measure columns to apply calibration to.
#'   If `NULL`, uses all measure columns.
#' @param calibration A calibration object or data frame containing the
#'   universal calibration curve (log(\[eta\]M) vs retention).
#' @param calibration_col If calibration is a column name in the data,
#'   specify it here.
#' @param intrinsic_visc_col Column containing intrinsic viscosity values.
#'   Required for converting between polymer types.
#' @param K_sample Mark-Houwink K parameter for the sample polymer.
#' @param a_sample Mark-Houwink a (alpha) exponent for the sample polymer.
#' @param K_standard Mark-Houwink K for the calibration standard polymer.
#'   Default is 0.000114 (polystyrene in THF).
#' @param a_standard Mark-Houwink a for the calibration standard.
#'   Default is 0.716 (polystyrene in THF).
#' @param output_col Name for the output molecular weight column.
#'   Default is `"mw_universal"`.
#' @param role Role for generated columns.
#' @param trained Logical indicating if the step has been trained.
#' @param skip Logical. Should the step be skipped when baking?
#' @param id Unique step identifier.
#'
#' @return An updated recipe with molecular weight calculated via universal
#'   calibration.
#'
#' @details
#' Universal calibration is based on the principle that polymers with the
#' same hydrodynamic volume elute at the same retention time, regardless of
#' chemical structure. The hydrodynamic volume is proportional to \[eta\]M:
#'
#' \deqn{V_h \propto [\eta] \cdot M}
#'
#' **The Universal Calibration Curve:**
#' \deqn{\log([\eta] \cdot M)_{sample} = \log([\eta] \cdot M)_{standard}}
#'
#' At the same retention volume, using Mark-Houwink equations:
#' \deqn{[\eta] = K \cdot M^a}
#'
#' We can solve for sample MW:
#' \deqn{M_{sample} = \left(\frac{K_{std} \cdot M_{std}^{1+a_{std}}}{K_{sample}}\right)^{\frac{1}{1+a_{sample}}}}
#'
#' **Mark-Houwink Parameters (THF, 25C):**
#' \itemize{
#'   \item Polystyrene: K = 0.000114, a = 0.716
#'   \item PMMA: K = 0.000128, a = 0.690
#'   \item Polyisoprene: K = 0.000251, a = 0.728
#'   \item Polybutadiene: K = 0.000457, a = 0.693
#' }
#'
#' @note
#' Universal calibration requires:
#' \itemize{
#'   \item Known Mark-Houwink parameters for both standard and sample
#'   \item Calibration with narrow standards (typically polystyrene)
#'   \item Same solvent and temperature for all measurements
#' }
#'
#' For absolute MW determination, consider using MALS detection instead.
#'
#' @family sec-calibration
#' @export
#'
#' @examples
#' \dontrun{
#' library(recipes)
#' library(measure)
#'
#' # Apply universal calibration to convert PS calibration to PMMA
#' rec <- recipe(~., data = pmma_data) |>
#'   step_measure_input_long(ri, location = vars(time), col_name = "ri") |>
#'   step_sec_baseline() |>
#'   step_sec_universal_cal(
#'     calibration = ps_calibration,
#'     K_sample = 0.000128,      # PMMA
#'     a_sample = 0.690,
#'     K_standard = 0.000114,    # PS (default)
#'     a_standard = 0.716
#'   ) |>
#'   prep()
#' }
step_sec_universal_cal <- function(
    recipe,
    measures = NULL,
    calibration = NULL,
    calibration_col = NULL,
    intrinsic_visc_col = NULL,
    K_sample,
    a_sample,
    K_standard = 0.000114,
    a_standard = 0.716,
    output_col = "mw_universal",
    role = NA,
    trained = FALSE,
    skip = FALSE,
    id = recipes::rand_id("sec_universal_cal")
) {
  # Validate Mark-Houwink parameters
  if (missing(K_sample) || missing(a_sample)) {
    cli::cli_abort(
      "Mark-Houwink parameters {.arg K_sample} and {.arg a_sample} are required."
    )
  }

  if (!is.numeric(K_sample) || K_sample <= 0) {
    cli::cli_abort("{.arg K_sample} must be a positive number.")
  }
  if (!is.numeric(a_sample) || a_sample <= 0 || a_sample > 1) {
    cli::cli_abort("{.arg a_sample} must be between 0 and 1.")
  }

  if (is.null(calibration) && is.null(calibration_col)) {
    cli::cli_abort(
      "Either {.arg calibration} or {.arg calibration_col} must be specified."
    )
  }

  recipes::add_step(
    recipe,
    step_sec_universal_cal_new(
      measures = measures,
      calibration = calibration,
      calibration_col = calibration_col,
      intrinsic_visc_col = intrinsic_visc_col,
      K_sample = K_sample,
      a_sample = a_sample,
      K_standard = K_standard,
      a_standard = a_standard,
      output_col = output_col,
      role = role,
      trained = trained,
      skip = skip,
      id = id
    )
  )
}

step_sec_universal_cal_new <- function(
    measures,
    calibration,
    calibration_col,
    intrinsic_visc_col,
    K_sample,
    a_sample,
    K_standard,
    a_standard,
    output_col,
    role,
    trained,
    skip,
    id
) {
  recipes::step(
    subclass = "sec_universal_cal",
    measures = measures,
    calibration = calibration,
    calibration_col = calibration_col,
    intrinsic_visc_col = intrinsic_visc_col,
    K_sample = K_sample,
    a_sample = a_sample,
    K_standard = K_standard,
    a_standard = a_standard,
    output_col = output_col,
    role = role,
    trained = trained,
    skip = skip,
    id = id
  )
}

#' @export
prep.step_sec_universal_cal <- function(x, training, info = NULL, ...) {
  check_for_measure(training)

  # Find measure columns if not specified
  if (is.null(x$measures)) {
    measures <- find_measure_cols(training)
  } else {
    measures <- x$measures
  }

  # Validate calibration
  calibration <- x$calibration
  if (!is.null(x$calibration_col)) {
    if (!x$calibration_col %in% names(training)) {
      cli::cli_abort(
        "Calibration column {.val {x$calibration_col}} not found."
      )
    }
  }

  step_sec_universal_cal_new(
    measures = measures,
    calibration = calibration,
    calibration_col = x$calibration_col,
    intrinsic_visc_col = x$intrinsic_visc_col,
    K_sample = x$K_sample,
    a_sample = x$a_sample,
    K_standard = x$K_standard,
    a_standard = x$a_standard,
    output_col = x$output_col,
    role = x$role,
    trained = TRUE,
    skip = x$skip,
    id = x$id
  )
}

#' Convert MW from standard to sample using Mark-Houwink
#' @noRd
.convert_mw_universal <- function(mw_std, K_std, a_std, K_sample, a_sample) {
  # At same retention volume: [eta]_std * M_std = [eta]_sample * M_sample
  # [eta] = K * M^a
  # K_std * M_std^(1+a_std) = K_sample * M_sample^(1+a_sample)
  # M_sample = (K_std * M_std^(1+a_std) / K_sample)^(1/(1+a_sample))

  log_hydrodynamic_vol <- log10(K_std) + (1 + a_std) * log10(mw_std)
  log_mw_sample <- (log_hydrodynamic_vol - log10(K_sample)) / (1 + a_sample)

  10^log_mw_sample
}

#' @export
bake.step_sec_universal_cal <- function(object, new_data, ...) {
  calibration <- object$calibration
  K_sample <- object$K_sample
  a_sample <- object$a_sample
  K_standard <- object$K_standard
  a_standard <- object$a_standard
  output_col <- object$output_col
  measures <- object$measures

  # Get calibration data (assumes it's a data frame with retention and mw columns)
  # or a fitted model that can predict log(M) from retention
  if (is.data.frame(calibration)) {
    # Assume columns: retention (or time/volume), mw (or log_mw)
    if ("log_mw" %in% names(calibration)) {
      cal_fit <- stats::lm(log_mw ~ poly(retention, 3), data = calibration)
    } else if ("mw" %in% names(calibration)) {
      calibration$log_mw <- log10(calibration$mw)
      cal_fit <- stats::lm(log_mw ~ poly(retention, 3), data = calibration)
    } else {
      cli::cli_abort(
        "Calibration data must contain 'mw' or 'log_mw' column."
      )
    }
  } else if (inherits(calibration, "lm")) {
    cal_fit <- calibration
  } else {
    cli::cli_abort(
      "Calibration must be a data frame or fitted model."
    )
  }

  # Calculate MW for each sample
  col <- measures[1]  # Use first measure column for location

  mw_list <- purrr::map(new_data[[col]], function(m) {
    location <- m$location

    # Predict log(MW) from calibration (for standard polymer)
    pred_data <- data.frame(retention = location)
    log_mw_std <- stats::predict(cal_fit, newdata = pred_data)
    mw_std <- 10^log_mw_std

    # Convert to sample polymer using universal calibration
    mw_sample <- .convert_mw_universal(
      mw_std, K_standard, a_standard, K_sample, a_sample
    )

    # Handle edge cases
    mw_sample[mw_sample < 100] <- NA_real_
    mw_sample[mw_sample > 1e10] <- NA_real_

    new_measure_tbl(location = location, value = mw_sample)
  })

  new_data[[output_col]] <- new_measure_list(mw_list)

  tibble::as_tibble(new_data)
}

#' @export
print.step_sec_universal_cal <- function(
    x,
    width = max(20, options()$width - 30),
    ...
) {
  title <- sprintf(
    "SEC universal calibration (K=%.2e, a=%.3f)",
    x$K_sample, x$a_sample
  )
  if (x$trained) {
    cat(title, " -> ", x$output_col, sep = "")
  } else {
    cat(title)
  }
  cat("\n")
  invisible(x)
}

#' @rdname tidy.step_sec
#' @export
#' @keywords internal
tidy.step_sec_universal_cal <- function(x, ...) {
  tibble::tibble(
    K_sample = x$K_sample,
    a_sample = x$a_sample,
    K_standard = x$K_standard,
    a_standard = x$a_standard,
    output_col = x$output_col,
    id = x$id
  )
}

#' @rdname required_pkgs.step_sec
#' @export
#' @keywords internal
required_pkgs.step_sec_universal_cal <- function(x, ...) {
  c("measure.sec", "measure")
}
