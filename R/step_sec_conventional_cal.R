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
#'   Required unless `calibration` is provided.
#' @param calibration A pre-loaded calibration object from
#'   [load_sec_calibration()]. When provided, skips fitting and uses the saved
#'   calibration directly. Takes precedence over `standards`.
#' @param fit_type Type of fit for the calibration curve:
#'   - `"cubic"` (default): Third-order polynomial
#'   - `"quadratic"`: Second-order polynomial
#'   - `"linear"`: First-order (linear) fit
#'   - `"fifth"`: Fifth-order polynomial
#'   - `"gam"`: Generalized Additive Model with cubic splines (requires mgcv)
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
  calibration = NULL,
  fit_type = c("cubic", "quadratic", "linear", "fifth", "gam"),
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

  # Check if using pre-loaded calibration
  if (!is.null(calibration)) {
    if (!inherits(calibration, "sec_calibration")) {
      cli::cli_abort(
        c(
          "{.arg calibration} must be a {.cls sec_calibration} object.",
          "i" = "Use {.fn load_sec_calibration} to load a saved calibration."
        )
      )
    }
    # Use settings from calibration object
    fit_type <- calibration$fit_type
    extrapolation <- calibration$extrapolation
    output_col <- calibration$output_col
    log_output <- calibration$log_output
    standards <- calibration$standards
  } else if (is.null(standards)) {
    cli::cli_abort(
      c(
        "Either {.arg standards} or {.arg calibration} is required.",
        "i" = "Provide calibration standards or use {.fn load_sec_calibration}."
      )
    )
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
      calibration = calibration,
      fit_type = fit_type,
      extrapolation = extrapolation,
      output_col = output_col,
      log_output = log_output,
      calibration_fit = NULL,
      calibration_range = NULL,
      calibration_diagnostics = NULL,
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

#' Fit GAM calibration curve using mgcv
#'
#' Fits a Generalized Additive Model with cubic regression splines for
#' flexible calibration curves. Falls back to linear model if too few points.
#'
#' @param cal_data Data frame with `location` and `log_mw` columns.
#' @return A fitted GAM model (mgcv::gam object) or lm object if fallback.
#' @noRd
.fit_gam_calibration <- function(cal_data) {
  rlang::check_installed("mgcv", reason = "for GAM calibration")

  n_unique <- length(unique(cal_data$location))

  if (n_unique < 4) {
    # Fall back to linear model for very few points
    cli::cli_warn(
      c(
        "Too few unique points ({n_unique}) for GAM, using linear model.",
        "i" = "GAM requires at least 4 unique calibration points."
      )
    )
    return(stats::lm(log_mw ~ location, data = cal_data))
  }

  # Set k (number of basis functions) based on available data

  # k must be less than number of unique values, and >= 3 for a spline
  # Use n_unique - 1 as upper bound, minimum of 4 for reasonable flexibility
  k <- min(n_unique - 1, 10)
  k <- max(k, 4)

  # Fit GAM with cubic regression splines (bs = 'cs')
  # Using REML for smoothing parameter estimation (most robust)
  mgcv::gam(
    stats::as.formula(paste0("log_mw ~ s(location, bs = 'cs', k = ", k, ")")),
    data = cal_data,
    method = "REML"
  )
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
  calibration,
  fit_type,
  extrapolation,
  output_col,
  log_output,
  calibration_fit,
  calibration_range,
  calibration_diagnostics,
  role,
  trained,
  skip,
  id
) {
  recipes::step(
    subclass = "sec_conventional_cal",
    measures = measures,
    standards = standards,
    calibration = calibration,
    fit_type = fit_type,
    extrapolation = extrapolation,
    output_col = output_col,
    log_output = log_output,
    calibration_fit = calibration_fit,
    calibration_range = calibration_range,
    calibration_diagnostics = calibration_diagnostics,
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

  # Check if using pre-loaded calibration
  if (!is.null(x$calibration)) {
    cal <- x$calibration

    # Use the pre-fitted model directly
    cal_fit <- cal$fit_object
    calibration_range <- cal$calibration_range
    calibration_diagnostics <- cal$diagnostics

    cli::cli_alert_info(
      "Using pre-loaded {.field {cal$fit_type}} calibration (R\\u00b2 = {round(cal$diagnostics$r_squared, 4)})"
    )

    return(
      step_sec_conventional_cal_new(
        measures = measures,
        standards = x$standards,
        calibration = x$calibration,
        fit_type = x$fit_type,
        extrapolation = x$extrapolation,
        output_col = x$output_col,
        log_output = x$log_output,
        calibration_fit = cal_fit,
        calibration_range = calibration_range,
        calibration_diagnostics = calibration_diagnostics,
        role = x$role,
        trained = TRUE,
        skip = x$skip,
        id = x$id
      )
    )
  }

  # Otherwise fit calibration from standards
  standards <- x$standards
  location_col <- .get_location_col(standards)
  location <- standards[[location_col]]
  log_mw <- .get_log_mw(standards)

  # Get MW in Daltons for % deviation calculation
  mw_daltons <- 10^log_mw

  # Validate sufficient standards for fit type
  n_standards <- length(location)
  min_standards <- switch(
    x$fit_type,
    linear = 2,
    quadratic = 3,
    cubic = 4,
    fifth = 6,
    gam = 4
  )

  if (n_standards < min_standards) {
    cli::cli_abort(
      c(
        "Insufficient standards for {.val {x$fit_type}} fit.",
        "i" = "Need at least {min_standards} standards, got {n_standards}."
      )
    )
  }

  # Fit calibration curve
  cal_data <- data.frame(location = location, log_mw = log_mw)

  if (x$fit_type == "gam") {
    # Fit GAM using mgcv with cubic splines
    cal_fit <- .fit_gam_calibration(cal_data)
  } else {
    # Fit polynomial using orthogonal polynomials
    degree <- switch(
      x$fit_type,
      linear = 1,
      quadratic = 2,
      cubic = 3,
      fifth = 5
    )
    cal_fit <- stats::lm(
      log_mw ~ stats::poly(location, degree),
      data = cal_data
    )
  }

  # Calculate comprehensive diagnostics
  cal_summary <- suppressWarnings(summary(cal_fit))

  if (x$fit_type == "gam") {
    # GAM summary has different structure
    r_squared <- cal_summary$r.sq
    adj_r_squared <- cal_summary$r.sq # GAM doesn't have adj R² in same way
    sigma <- sqrt(cal_summary$scale) # Residual standard error
  } else {
    r_squared <- cal_summary$r.squared
    adj_r_squared <- cal_summary$adj.r.squared
    sigma <- cal_summary$sigma # Residual standard error
  }

  # Predictions and residuals for each standard
  predicted_log_mw <- stats::predict(cal_fit)
  residuals_log_mw <- log_mw - predicted_log_mw

  # Calculate MW from predictions
  predicted_mw <- 10^predicted_log_mw

  # Calculate % deviation in MW space (more intuitive for scientists)
  pct_deviation <- 100 * (predicted_mw - mw_daltons) / mw_daltons

  # RMSE in log(MW) space
  rmse_log <- sqrt(mean(residuals_log_mw^2))

  # Calculate prediction intervals at standard locations
  if (x$fit_type == "gam") {
    # GAM uses se.fit differently
    pred_with_se <- stats::predict(cal_fit, newdata = cal_data, se.fit = TRUE)
    prediction_se <- pred_with_se$se.fit
    # Approximate prediction intervals using 1.96 * SE
    ci_lower_log_mw <- as.numeric(predicted_log_mw) - 1.96 * prediction_se
    ci_upper_log_mw <- as.numeric(predicted_log_mw) + 1.96 * prediction_se
  } else {
    pred_with_se <- stats::predict(
      cal_fit,
      newdata = cal_data,
      se.fit = TRUE,
      interval = "prediction",
      level = 0.95
    )
    prediction_se <- pred_with_se$se.fit
    ci_lower_log_mw <- as.numeric(pred_with_se$fit[, "lwr"])
    ci_upper_log_mw <- as.numeric(pred_with_se$fit[, "upr"])
  }

  # Build per-standard diagnostics table
  standard_diagnostics <- tibble::tibble(
    location = location,
    actual_log_mw = log_mw,
    predicted_log_mw = as.numeric(predicted_log_mw),
    residual_log_mw = as.numeric(residuals_log_mw),
    actual_mw = mw_daltons,
    predicted_mw = as.numeric(predicted_mw),
    pct_deviation = as.numeric(pct_deviation),
    prediction_se = as.numeric(prediction_se),
    ci_lower_log_mw = ci_lower_log_mw,
    ci_upper_log_mw = ci_upper_log_mw
  )

  # Store calibration range
  calibration_range <- range(location)

  # Build comprehensive diagnostics object
  # GAM uses different df accessor

  if (x$fit_type == "gam") {
    df_residual <- cal_fit$df.residual
  } else {
    df_residual <- cal_fit$df.residual
  }

  calibration_diagnostics <- list(
    r_squared = r_squared,
    adj_r_squared = adj_r_squared,
    rmse_log_mw = rmse_log,
    residual_std_error = sigma,
    degrees_of_freedom = df_residual,
    max_abs_residual = max(abs(residuals_log_mw)),
    max_abs_pct_deviation = max(abs(pct_deviation)),
    mean_abs_pct_deviation = mean(abs(pct_deviation)),
    standard_results = standard_diagnostics
  )

  # Warn if R² is low
  if (r_squared < 0.99) {
    cli::cli_warn(
      c(
        "Calibration fit quality is low (R\\u00b2 = {round(r_squared, 4)}).",
        "i" = "Typical calibrations have R\\u00b2 > 0.999."
      )
    )
  }

  # Warn if any standard has high deviation
  if (max(abs(pct_deviation)) > 5) {
    worst_idx <- which.max(abs(pct_deviation))
    cli::cli_warn(
      c(
        "Standard at {round(location[worst_idx], 2)} has {round(pct_deviation[worst_idx], 1)}% MW deviation.",
        "i" = "Consider removing outlier standards or using a different fit type."
      )
    )
  }

  step_sec_conventional_cal_new(
    measures = measures,
    standards = standards,
    calibration = NULL,
    fit_type = x$fit_type,
    extrapolation = x$extrapolation,
    output_col = x$output_col,
    log_output = x$log_output,
    calibration_fit = cal_fit,
    calibration_range = calibration_range,
    calibration_diagnostics = calibration_diagnostics,
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

  if (x$trained && !is.null(x$calibration_diagnostics)) {
    diag <- x$calibration_diagnostics
    title <- sprintf(
      "%s, R\u00b2=%.4f, RMSE=%.4f",
      title,
      diag$r_squared,
      diag$rmse_log_mw
    )
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
  if (x$trained && !is.null(x$calibration_diagnostics)) {
    diag <- x$calibration_diagnostics

    # Get coefficients - GAM stores these differently
    if (x$fit_type == "gam") {
      coefs <- x$calibration_fit$coefficients
    } else {
      coefs <- stats::coef(x$calibration_fit)
    }

    tibble::tibble(
      fit_type = x$fit_type,
      n_standards = nrow(x$standards),
      r_squared = diag$r_squared,
      adj_r_squared = diag$adj_r_squared,
      rmse_log_mw = diag$rmse_log_mw,
      residual_std_error = diag$residual_std_error,
      max_abs_pct_deviation = diag$max_abs_pct_deviation,
      mean_abs_pct_deviation = diag$mean_abs_pct_deviation,
      calibration_min = x$calibration_range[1],
      calibration_max = x$calibration_range[2],
      degrees_of_freedom = diag$degrees_of_freedom,
      coefficients = list(coefs),
      standard_results = list(diag$standard_results),
      output_col = x$output_col,
      id = x$id
    )
  } else {
    tibble::tibble(
      fit_type = x$fit_type,
      n_standards = if (is.null(x$standards)) {
        NA_integer_
      } else {
        nrow(x$standards)
      },
      r_squared = NA_real_,
      adj_r_squared = NA_real_,
      rmse_log_mw = NA_real_,
      residual_std_error = NA_real_,
      max_abs_pct_deviation = NA_real_,
      mean_abs_pct_deviation = NA_real_,
      calibration_min = NA_real_,
      calibration_max = NA_real_,
      degrees_of_freedom = NA_integer_,
      coefficients = list(NULL),
      standard_results = list(NULL),
      output_col = x$output_col,
      id = x$id
    )
  }
}

#' @rdname required_pkgs.step_sec
#' @export
#' @keywords internal
required_pkgs.step_sec_conventional_cal <- function(x, ...) {
  pkgs <- c("measure.sec", "measure")
  if (!is.null(x$fit_type) && x$fit_type == "gam") {
    pkgs <- c(pkgs, "mgcv")
  }
  pkgs
}
