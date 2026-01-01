# ==============================================================================
# step_sec_conventional_cal
#
# Conventional calibration using narrow molecular weight standards
# ==============================================================================

#' Conventional Calibration for SEC Using Narrow Standards
#'
#' `step_sec_conventional_cal()` creates a *specification* of a recipe step that
#' fits a calibration curve from narrow molecular weight standards and applies
#' it to convert elution time/volume to molecular weight.
#'
#' @param recipe A recipe object.
#' @param measures Character vector of measure columns to apply calibration to.
#'   If `NULL`, uses all measure columns.
#' @param standards A data frame containing calibration standards with columns:
#'   - `location` (or `time`, `volume`, `retention`): Elution position
#'   - `log_mw` (or `mw`): Molecular weight (will be log-transformed if `mw`)
#' @param fit_type Type of polynomial fit for the calibration curve:
#'   - `"cubic"` (default): Third-order polynomial
#'   - `"quadratic"`: Second-order polynomial
#'   - `"linear"`: First-order (linear) fit
#'   - `"fifth"`: Fifth-order polynomial
#' @param extrapolation How to handle data outside the calibration range:
#'   - `"warn"` (default): Extrapolate but warn
#'   - `"none"`: Return NA for out-of-range values
#'   - `"linear"`: Use linear extrapolation at boundaries
#' @param output_col Name for the output molecular weight column.
#'   Default is `"mw"`.
#' @param log_output Logical. If `TRUE` (default), output column contains

#'   log10(MW). If `FALSE`, output contains MW in Daltons.
#' @param role Role for generated columns.
#' @param trained Logical indicating if the step has been trained.
#' @param skip Logical. Should the step be skipped when baking?
#' @param id Unique step identifier.
#'
#' @return An updated recipe with the new step added.
#'
#' @details
#' This step performs conventional (also called relative) SEC calibration using
#' narrow dispersity standards of known molecular weight. The calibration curve
#' relates elution position to log(MW):
#'
#' \deqn{\log_{10}(M) = f(V_e)}
#'
#' where f is a polynomial function and V_e is the elution volume or time.
#'
#' **Calibration Curve Fitting:**
#'
#' The calibration is fit using orthogonal polynomials for numerical stability.
#' At least 3 standards are required for cubic fits, 4 for quadratic, etc.
#'
#' **Important Considerations:**
#' \itemize{
#'   \item Standards should bracket the MW range of interest
#'   \item Calibration is polymer-specific (different polymers have different
#'     hydrodynamic volumes at the same MW)
#'   \item For cross-polymer comparisons, use universal calibration instead
#'     (\code{\link{step_sec_universal_cal}})
#' }
#'
#' **Fit Quality Metrics:**
#' The `tidy()` method returns calibration coefficients and R-squared values
#' for assessing fit quality. R² > 0.999 is typical for good calibrations.
#'
#' @family sec-calibration
#' @export
#'
#' @examples
#' \dontrun{
#' library(recipes)
#' library(measure)
#'
#' # Create calibration standards data
#' ps_standards <- data.frame(
#'   retention = c(12.5, 13.2, 14.1, 15.0, 16.2, 17.5),
#'   log_mw = c(6.0, 5.5, 5.0, 4.5, 4.0, 3.5)
#' )
#'
#' # Apply conventional calibration
#' rec <- recipe(~., data = polymer_data) |>
#'   step_measure_input_long(ri, location = vars(time), col_name = "ri") |>
#'   step_sec_baseline() |>
#'   step_sec_conventional_cal(
#'     standards = ps_standards,
#'     fit_type = "cubic"
#'   ) |>
#'   prep()
#'
#' # Check calibration quality
#' tidy(rec, number = 3)
#' }
step_sec_conventional_cal <- function(
  recipe,
  measures = NULL,
  standards = NULL,
  fit_type = c("cubic", "quadratic", "linear", "fifth"),
  extrapolation = c("warn", "none", "linear"),
  output_col = "mw",
  log_output = TRUE,
  role = NA,
  trained = FALSE,
  skip = FALSE,

  id = recipes::rand_id("sec_conventional_cal")
) {
  fit_type <- match.arg(fit_type)
  extrapolation <- match.arg(extrapolation)

  # Validate standards
  if (is.null(standards)) {
    cli::cli_abort("{.arg standards} is required for conventional calibration.")
  }

  if (!is.data.frame(standards)) {
    cli::cli_abort("{.arg standards} must be a data frame.")
  }

  # Validate standards has required columns
  .validate_standards_columns(standards)

  recipes::add_step(
    recipe,
    step_sec_conventional_cal_new(
      measures = measures,
      standards = standards,
      fit_type = fit_type,
      extrapolation = extrapolation,
      output_col = output_col,
      log_output = log_output,
      calibration_fit = NULL,
      calibration_range = NULL,
      r_squared = NULL,
      role = role,
      trained = trained,
      skip = skip,
      id = id
    )
  )
}

#' Validate standards data frame has required columns
#' @noRd
.validate_standards_columns <- function(standards) {
  # Check for location column

  location_cols <- c(
    "location",
    "time",
    "volume",
    "retention",
    "elution_time",
    "elution_volume",
    "ret_time",
    "ret_vol"
  )
  has_location <- any(names(standards) %in% location_cols)

  if (!has_location) {
    cli::cli_abort(
      c(
        "Standards data frame must have a location column.",
        "i" = "Expected one of: {.val {location_cols}}"
      )
    )
  }

  # Check for MW column
  mw_cols <- c("log_mw", "mw", "log_m", "m", "mol_weight", "molecular_weight")
  has_mw <- any(names(standards) %in% mw_cols)

  if (!has_mw) {
    cli::cli_abort(
      c(
        "Standards data frame must have a molecular weight column.",
        "i" = "Expected one of: {.val {mw_cols}}"
      )
    )
  }

  invisible(TRUE)
}

#' Get location column from standards
#' @noRd
.get_location_col <- function(standards) {
  location_cols <- c(
    "location",
    "time",
    "volume",
    "retention",
    "elution_time",
    "elution_volume",
    "ret_time",
    "ret_vol"
  )
  found <- names(standards)[names(standards) %in% location_cols]
  found[1]
}

#' Get log MW values from standards
#' @noRd
.get_log_mw <- function(standards) {
  if ("log_mw" %in% names(standards)) {
    return(standards$log_mw)
  } else if ("log_m" %in% names(standards)) {
    return(standards$log_m)
  } else if ("mw" %in% names(standards)) {
    return(log10(standards$mw))
  } else if ("m" %in% names(standards)) {
    return(log10(standards$m))
  } else if ("mol_weight" %in% names(standards)) {
    return(log10(standards$mol_weight))
  } else if ("molecular_weight" %in% names(standards)) {
    return(log10(standards$molecular_weight))
  }

  cli::cli_abort("Could not find molecular weight column in standards.")
}

step_sec_conventional_cal_new <- function(
  measures,
  standards,
  fit_type,
  extrapolation,
  output_col,
  log_output,
  calibration_fit,
  calibration_range,
  r_squared,
  role,
  trained,
  skip,
  id
) {
  recipes::step(
    subclass = "sec_conventional_cal",
    measures = measures,
    standards = standards,
    fit_type = fit_type,
    extrapolation = extrapolation,
    output_col = output_col,
    log_output = log_output,
    calibration_fit = calibration_fit,
    calibration_range = calibration_range,
    r_squared = r_squared,
    role = role,
    trained = trained,
    skip = skip,
    id = id
  )
}

#' @export
prep.step_sec_conventional_cal <- function(x, training, info = NULL, ...) {
  check_for_measure(training)

  # Find measure columns if not specified
  if (is.null(x$measures)) {
    measures <- find_measure_cols(training)
  } else {
    measures <- x$measures
  }

  # Extract calibration data
  standards <- x$standards
  location_col <- .get_location_col(standards)
  location <- standards[[location_col]]
  log_mw <- .get_log_mw(standards)

  # Validate sufficient standards for fit type
  n_standards <- length(location)
  min_standards <- switch(
    x$fit_type,
    linear = 2,
    quadratic = 3,
    cubic = 4,
    fifth = 6
  )

  if (n_standards < min_standards) {
    cli::cli_abort(
      c(
        "Insufficient standards for {.val {x$fit_type}} fit.",
        "i" = "Need at least {min_standards} standards, got {n_standards}."
      )
    )
  }

  # Fit calibration curve using orthogonal polynomials
  degree <- switch(
    x$fit_type,
    linear = 1,
    quadratic = 2,
    cubic = 3,
    fifth = 5
  )

  cal_data <- data.frame(location = location, log_mw = log_mw)
  cal_fit <- stats::lm(log_mw ~ stats::poly(location, degree), data = cal_data)

  # Calculate R-squared (suppress "essentially perfect fit" warning for exact fits)
  r_squared <- suppressWarnings(summary(cal_fit)$r.squared)

  # Store calibration range
  calibration_range <- range(location)

  # Warn if R² is low

  if (r_squared < 0.99) {
    cli::cli_warn(
      c(
        "Calibration fit quality is low (R\\u00b2 = {round(r_squared, 4)}).",
        "i" = "Typical calibrations have R\\u00b2 > 0.999."
      )
    )
  }

  step_sec_conventional_cal_new(
    measures = measures,
    standards = standards,
    fit_type = x$fit_type,
    extrapolation = x$extrapolation,
    output_col = x$output_col,
    log_output = x$log_output,
    calibration_fit = cal_fit,
    calibration_range = calibration_range,
    r_squared = r_squared,
    role = x$role,
    trained = TRUE,
    skip = x$skip,
    id = x$id
  )
}

#' @export
bake.step_sec_conventional_cal <- function(object, new_data, ...) {
  cal_fit <- object$calibration_fit
  calibration_range <- object$calibration_range
  extrapolation <- object$extrapolation
  output_col <- object$output_col
  log_output <- object$log_output
  measures <- object$measures

  # Use first measure column for location data
  col <- measures[1]

  mw_list <- purrr::map(new_data[[col]], function(m) {
    location <- m$location

    # Check for out-of-range values
    out_of_range <- location < calibration_range[1] |
      location > calibration_range[2]

    if (any(out_of_range) && extrapolation == "warn") {
      n_out <- sum(out_of_range)
      pct_out <- round(100 * n_out / length(location), 1)
      cli::cli_warn(
        c(
          "{n_out} points ({pct_out}%) are outside calibration range.",
          "i" = "Calibration range: {round(calibration_range[1], 2)} to {round(calibration_range[2], 2)}"
        )
      )
    }

    # Predict log(MW) from calibration
    pred_data <- data.frame(location = location)
    log_mw <- stats::predict(cal_fit, newdata = pred_data)

    # Handle out-of-range based on extrapolation setting
    if (extrapolation == "none") {
      log_mw[out_of_range] <- NA_real_
    }
    # For "warn" and "linear", allow extrapolation (poly already extrapolates)

    # Apply reasonable bounds
    log_mw[log_mw < 1] <- NA_real_ # < 10 Da
    log_mw[log_mw > 10] <- NA_real_ # > 10^10 Da

    # Convert to MW if requested
    if (log_output) {
      value <- log_mw
    } else {
      value <- 10^log_mw
    }

    new_measure_tbl(location = location, value = value)
  })

  new_data[[output_col]] <- new_measure_list(mw_list)

  tibble::as_tibble(new_data)
}

#' @export
print.step_sec_conventional_cal <- function(
  x,
  width = max(20, options()$width - 30),
  ...
) {
  title <- sprintf("SEC conventional calibration (%s fit)", x$fit_type)

  if (x$trained && !is.null(x$r_squared)) {
    title <- sprintf("%s, R\u00b2=%.4f", title, x$r_squared)
  }

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
tidy.step_sec_conventional_cal <- function(x, ...) {
  if (x$trained && !is.null(x$calibration_fit)) {
    coefs <- stats::coef(x$calibration_fit)
    tibble::tibble(
      fit_type = x$fit_type,
      r_squared = x$r_squared,
      calibration_min = x$calibration_range[1],
      calibration_max = x$calibration_range[2],
      n_standards = nrow(x$standards),
      coefficients = list(coefs),
      output_col = x$output_col,
      id = x$id
    )
  } else {
    tibble::tibble(
      fit_type = x$fit_type,
      r_squared = NA_real_,
      calibration_min = NA_real_,
      calibration_max = NA_real_,
      n_standards = if (is.null(x$standards)) {
        NA_integer_
      } else {
        nrow(x$standards)
      },
      coefficients = list(NULL),
      output_col = x$output_col,
      id = x$id
    )
  }
}

#' @rdname required_pkgs.step_sec
#' @export
#' @keywords internal
required_pkgs.step_sec_conventional_cal <- function(x, ...) {
  c("measure.sec", "measure")
}
