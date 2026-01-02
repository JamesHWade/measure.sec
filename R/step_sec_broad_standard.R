# ==============================================================================
# step_sec_broad_standard
#
# Calibration using broad/polydisperse molecular weight standards
# ==============================================================================

#' Broad Standard Calibration for SEC/GPC
#'
#' `step_sec_broad_standard()` creates a *specification* of a recipe step that
#' fits a calibration curve using a polydisperse (broad) molecular weight
#' standard with known Mn and Mw values.
#'
#' @param recipe A recipe object.
#' @param measures Character vector of measure columns to apply calibration to.
#'   If `NULL`, uses all measure columns.
#' @param broad_standard A data frame containing the broad standard chromatogram:
#'   - `location` (or `time`, `volume`, `retention`, `elution_time`,
#'     `elution_volume`): Elution position
#'   - `value` (or `signal`, `response`, `intensity`, `ri`, `uv`): Detector response
#' @param known_mn Known number-average molecular weight (Mn) in Daltons.
#' @param known_mw Known weight-average molecular weight (Mw) in Daltons.
#' @param fit_type Type of calibration curve:
#'   - `"linear"` (default): log10(M) = C1 + C2*V (classic Hamielec)
#'   - `"quadratic"`: log10(M) = C1 + C2*V + C3*V^2
#' @param method Calibration method:
#'   - `"hamielec"` (default): Optimize to match Mn and Mw
#'   - `"integral"`: Use integral MWD matching (not yet implemented)
#' @param integration_range Optional numeric vector `c(min, max)` specifying
#'   the elution range to use for the broad standard. If `NULL`, auto-detects
#'   peak region.
#' @param extrapolation How to handle data outside the calibration range:
#'   - `"warn"` (default): Extrapolate but warn
#'   - `"none"`: Return NA for out-of-range values
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
#' This step implements broad standard calibration methods for SEC/GPC, which
#' use a single polydisperse standard with known molecular weight averages
#' instead of multiple narrow standards.
#'
#' **Hamielec Method:**
#'
#' Based on Balke, Hamielec, LeClair, and Pearce (1969). The original method
#' assumes a linear calibration relationship, though this implementation also
#' supports a quadratic extension for better fit over wide MW ranges:
#'
#' \deqn{\log_{10}(M) = C_1 + C_2 \cdot V}
#'
#' The algorithm simultaneously optimizes coefficients C1 and C2 (and C3 for
#' quadratic fits) using Nelder-Mead optimization to minimize the squared
#' relative errors between calculated and known Mn and Mw values.
#'
#' **When to Use Broad Standard Calibration:**
#'
#' - QC labs running the same polymer type repeatedly
#' - When narrow standards aren't available for your polymer
#' - When you have a well-characterized in-house reference material
#' - Provides "absolute" molecular weights for the same polymer type
#'
#' **Limitations:**
#'
#' - Results are only valid for polymers with similar hydrodynamic behavior
#' - Linear calibration may not fit well over very wide MW ranges
#' - Requires well-characterized broad standard (accurate Mn and Mw)
#'
#' @references
#' Balke, S.T., Hamielec, A.E., LeClair, B.P., and Pearce, S.L. (1969).
#' Gel permeation chromatography.
#' *Industrial & Engineering Chemistry Product Research and Development*,
#' 8(1), 54-57.
#'
#' @family sec-calibration
#' @export
#'
#' @examples
#' \dontrun{
#' library(recipes)
#' library(measure)
#'
#' data(sec_triple_detect)
#'
#' # Create broad standard chromatogram data
#' broad_std <- data.frame(
#'   time = seq(10, 20, by = 0.1),
#'   signal = dnorm(seq(10, 20, by = 0.1), mean = 15, sd = 1.5)
#' )
#'
#' # Apply broad standard calibration
#' rec <- recipe(~., data = sec_triple_detect) |>
#'   step_measure_input_long(
#'     ri_signal,
#'     location = vars(elution_time),
#'     col_name = "ri"
#'   ) |>
#'   step_sec_baseline() |>
#'   step_sec_broad_standard(
#'     broad_standard = broad_std,
#'     known_mn = 50000,
#'     known_mw = 150000
#'   ) |>
#'   prep()
#'
#' # Check calibration results
#' tidy(rec, number = 3)
#' }
step_sec_broad_standard <- function(
  recipe,
  measures = NULL,
  broad_standard = NULL,
  known_mn = NULL,
  known_mw = NULL,
  fit_type = c("linear", "quadratic"),
  method = c("hamielec", "integral"),
  integration_range = NULL,
  extrapolation = c("warn", "none"),
  output_col = "mw",
  log_output = TRUE,
  role = NA,

  trained = FALSE,
  skip = FALSE,
  id = recipes::rand_id("sec_broad_standard")
) {
  fit_type <- match.arg(fit_type)
  method <- match.arg(method)
  extrapolation <- match.arg(extrapolation)

  # Validate required parameters
  if (is.null(broad_standard)) {
    cli::cli_abort(
      "{.arg broad_standard} is required for broad standard calibration."
    )
  }

  if (!is.data.frame(broad_standard)) {
    cli::cli_abort("{.arg broad_standard} must be a data frame.")
  }

  if (is.null(known_mn) || is.null(known_mw)) {
    cli::cli_abort("Both {.arg known_mn} and {.arg known_mw} are required.")
  }

  if (!is.numeric(known_mn) || known_mn <= 0) {
    cli::cli_abort("{.arg known_mn} must be a positive number.")
  }

  if (!is.numeric(known_mw) || known_mw <= 0) {
    cli::cli_abort("{.arg known_mw} must be a positive number.")
  }

  if (known_mw < known_mn) {
    cli::cli_abort(
      "{.arg known_mw} must be greater than or equal to {.arg known_mn}."
    )
  }

  # Validate broad standard has required columns
  .validate_broad_standard_columns(broad_standard)

  # Warn about unimplemented methods
  if (method == "integral") {
    cli::cli_warn(
      c(
        "The {.val integral} method is not yet implemented.",
        "i" = "Using {.val hamielec} method instead."
      )
    )
    method <- "hamielec"
  }

  recipes::add_step(
    recipe,
    step_sec_broad_standard_new(
      measures = measures,
      broad_standard = broad_standard,
      known_mn = known_mn,
      known_mw = known_mw,
      fit_type = fit_type,
      method = method,
      integration_range = integration_range,
      extrapolation = extrapolation,
      output_col = output_col,
      log_output = log_output,
      calibration_coefficients = NULL,
      calibration_range = NULL,
      calibration_diagnostics = NULL,
      role = role,
      trained = trained,
      skip = skip,
      id = id
    )
  )
}

#' Validate broad standard data frame has required columns
#' @noRd
.validate_broad_standard_columns <- function(broad_standard) {
  # Check for location column
  location_cols <- c(
    "location",
    "time",
    "volume",
    "retention",
    "elution_time",
    "elution_volume"
  )
  has_location <- any(names(broad_standard) %in% location_cols)

  if (!has_location) {
    cli::cli_abort(
      c(
        "Broad standard data frame must have a location column.",
        "i" = "Expected one of: {.val {location_cols}}"
      )
    )
  }

  # Check for signal column
  signal_cols <- c("value", "signal", "response", "intensity", "ri", "uv")
  has_signal <- any(names(broad_standard) %in% signal_cols)

  if (!has_signal) {
    cli::cli_abort(
      c(
        "Broad standard data frame must have a signal column.",
        "i" = "Expected one of: {.val {signal_cols}}"
      )
    )
  }

  invisible(TRUE)
}

#' Get location column from broad standard
#' @noRd
.get_broad_location_col <- function(broad_standard) {
  location_cols <- c(
    "location",
    "time",
    "volume",
    "retention",
    "elution_time",
    "elution_volume"
  )
  found <- names(broad_standard)[names(broad_standard) %in% location_cols]
  found[1]
}

#' Get signal column from broad standard
#' @noRd
.get_broad_signal_col <- function(broad_standard) {
  signal_cols <- c("value", "signal", "response", "intensity", "ri", "uv")
  found <- names(broad_standard)[names(broad_standard) %in% signal_cols]
  found[1]
}

step_sec_broad_standard_new <- function(
  measures,
  broad_standard,
  known_mn,
  known_mw,
  fit_type,
  method,
  integration_range,
  extrapolation,
  output_col,
  log_output,
  calibration_coefficients,
  calibration_range,
  calibration_diagnostics,
  role,
  trained,
  skip,
  id
) {
  recipes::step(
    subclass = "sec_broad_standard",
    measures = measures,
    broad_standard = broad_standard,
    known_mn = known_mn,
    known_mw = known_mw,
    fit_type = fit_type,
    method = method,
    integration_range = integration_range,
    extrapolation = extrapolation,
    output_col = output_col,
    log_output = log_output,
    calibration_coefficients = calibration_coefficients,
    calibration_range = calibration_range,
    calibration_diagnostics = calibration_diagnostics,
    role = role,
    trained = trained,
    skip = skip,
    id = id
  )
}

#' Calculate Mn and Mw for given calibration coefficients
#' @noRd
.calc_mw_from_coefficients <- function(
  location,
  signal,
  coefficients,
  fit_type
) {
  # Apply calibration to get log10(MW)
  if (fit_type == "linear") {
    log_mw <- coefficients[1] + coefficients[2] * location
  } else {
    # quadratic
    log_mw <- coefficients[1] +
      coefficients[2] * location +
      coefficients[3] * location^2
  }

  # Convert to MW
  mw <- 10^log_mw

  # Weight is proportional to signal (ensure non-negative)
  w <- pmax(signal, 0)

  # Only use positive weights
  valid <- w > 0 & is.finite(mw) & mw > 0
  if (sum(valid) < 2) {
    return(list(mn = NA_real_, mw = NA_real_, dispersity = NA_real_))
  }

  mw <- mw[valid]
  w <- w[valid]
  sum_w <- sum(w)

  # Calculate Mn and Mw
  mn <- sum_w / sum(w / mw)
  mw_calc <- sum(w * mw) / sum_w
  dispersity <- mw_calc / mn

  list(mn = mn, mw = mw_calc, dispersity = dispersity)
}

#' Objective function for Hamielec optimization
#' @noRd
.hamielec_objective <- function(
  params,
  location,
  signal,
  known_mn,
  known_mw,
  fit_type
) {
  result <- .calc_mw_from_coefficients(location, signal, params, fit_type)

  if (is.na(result$mn) || is.na(result$mw)) {
    return(1e10)
  }

  # Objective: minimize squared relative errors in Mn and Mw
  err_mn <- ((result$mn - known_mn) / known_mn)^2
  err_mw <- ((result$mw - known_mw) / known_mw)^2

  err_mn + err_mw
}

#' Fit Hamielec calibration
#' @noRd
.fit_hamielec <- function(
  location,
  signal,
  known_mn,
  known_mw,
  fit_type,
  integration_range
) {
  # Apply integration range if specified
  if (!is.null(integration_range)) {
    idx <- location >= integration_range[1] & location <= integration_range[2]
    location <- location[idx]
    signal <- signal[idx]
  } else {
    # Auto-detect peak region (above 5% of max)
    threshold <- 0.05 * max(signal, na.rm = TRUE)
    idx <- signal > threshold
    if (sum(idx) > 10) {
      location <- location[idx]
      signal <- signal[idx]
    }
  }

  if (length(location) < 10) {
    cli::cli_abort("Insufficient data points in broad standard chromatogram.")
  }

  known_dispersity <- known_mw / known_mn
  known_log_mw <- log10(known_mw)
  known_log_mn <- log10(known_mn)

  # Initial guess for coefficients
  # Typical SEC: log(MW) decreases with elution time
  # Rough estimate: MW range spans calibration range
  v_range <- range(location)
  v_mid <- mean(v_range)

  if (fit_type == "linear") {
    # log10(M) = C1 + C2*V
    # At v_mid, estimate log10(MW) as average of log10(Mn) and log10(Mw)
    # This equals log10(geometric_mean(Mn, Mw)), not log10(arithmetic_mean)
    avg_log_mw <- (known_log_mn + known_log_mw) / 2

    # Slope is typically negative (higher MW elutes first)
    # Estimate slope from MW range and elution range
    c2_init <- -(known_log_mw - known_log_mn) / diff(v_range) * 2
    c1_init <- avg_log_mw - c2_init * v_mid

    initial_params <- c(c1_init, c2_init)

    # Optimize
    result <- stats::optim(
      par = initial_params,
      fn = .hamielec_objective,
      location = location,
      signal = signal,
      known_mn = known_mn,
      known_mw = known_mw,
      fit_type = fit_type,
      method = "Nelder-Mead",
      control = list(maxit = 1000)
    )

    coefficients <- result$par
    names(coefficients) <- c("intercept", "slope")
  } else {
    # Quadratic fit
    avg_log_mw <- (known_log_mn + known_log_mw) / 2
    c2_init <- -(known_log_mw - known_log_mn) / diff(v_range) * 2
    c1_init <- avg_log_mw - c2_init * v_mid
    c3_init <- 0

    initial_params <- c(c1_init, c2_init, c3_init)

    result <- stats::optim(
      par = initial_params,
      fn = .hamielec_objective,
      location = location,
      signal = signal,
      known_mn = known_mn,
      known_mw = known_mw,
      fit_type = fit_type,
      method = "Nelder-Mead",
      control = list(maxit = 2000)
    )

    coefficients <- result$par
    names(coefficients) <- c("intercept", "slope", "quadratic")
  }

  # Calculate final MW values with fitted coefficients
  final_result <- .calc_mw_from_coefficients(
    location,
    signal,
    coefficients,
    fit_type
  )

  # Calculate diagnostics
  mn_error_pct <- 100 * (final_result$mn - known_mn) / known_mn
  mw_error_pct <- 100 * (final_result$mw - known_mw) / known_mw
  dispersity_error_pct <- 100 *
    (final_result$dispersity - known_dispersity) /
    known_dispersity

  diagnostics <- list(
    known_mn = known_mn,
    known_mw = known_mw,
    known_dispersity = known_dispersity,
    calculated_mn = final_result$mn,
    calculated_mw = final_result$mw,
    calculated_dispersity = final_result$dispersity,
    mn_error_pct = mn_error_pct,
    mw_error_pct = mw_error_pct,
    dispersity_error_pct = dispersity_error_pct,
    convergence = result$convergence,
    iterations = result$counts["function"],
    final_objective = result$value
  )

  # Warn if fit is poor
  if (abs(mn_error_pct) > 5 || abs(mw_error_pct) > 5) {
    cli::cli_warn(
      c(
        "Broad standard calibration fit has >5% error.",
        "i" = "Mn error: {round(mn_error_pct, 2)}%, Mw error: {round(mw_error_pct, 2)}%",
        "i" = "Consider adjusting integration range or using different fit type."
      )
    )
  }

  list(
    coefficients = coefficients,
    calibration_range = range(location),
    diagnostics = diagnostics
  )
}

#' @export
prep.step_sec_broad_standard <- function(x, training, info = NULL, ...) {
  check_for_measure(training)

  # Find measure columns if not specified
  if (is.null(x$measures)) {
    measures <- find_measure_cols(training)
  } else {
    measures <- x$measures
  }

  # Extract broad standard data
  broad_standard <- x$broad_standard
  location_col <- .get_broad_location_col(broad_standard)
  signal_col <- .get_broad_signal_col(broad_standard)

  location <- broad_standard[[location_col]]
  signal <- broad_standard[[signal_col]]

  # Fit calibration
  fit_result <- .fit_hamielec(
    location = location,
    signal = signal,
    known_mn = x$known_mn,
    known_mw = x$known_mw,
    fit_type = x$fit_type,
    integration_range = x$integration_range
  )

  step_sec_broad_standard_new(
    measures = measures,
    broad_standard = x$broad_standard,
    known_mn = x$known_mn,
    known_mw = x$known_mw,
    fit_type = x$fit_type,
    method = x$method,
    integration_range = x$integration_range,
    extrapolation = x$extrapolation,
    output_col = x$output_col,
    log_output = x$log_output,
    calibration_coefficients = fit_result$coefficients,
    calibration_range = fit_result$calibration_range,
    calibration_diagnostics = fit_result$diagnostics,
    role = x$role,
    trained = TRUE,
    skip = x$skip,
    id = x$id
  )
}

#' @export
bake.step_sec_broad_standard <- function(object, new_data, ...) {
  coefficients <- object$calibration_coefficients
  calibration_range <- object$calibration_range
  extrapolation <- object$extrapolation
  output_col <- object$output_col
  log_output <- object$log_output
  fit_type <- object$fit_type
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

    # Apply calibration to get log10(MW)
    if (fit_type == "linear") {
      log_mw <- coefficients["intercept"] + coefficients["slope"] * location
    } else {
      log_mw <- coefficients["intercept"] +
        coefficients["slope"] * location +
        coefficients["quadratic"] * location^2
    }

    # Handle out-of-range based on extrapolation setting
    if (extrapolation == "none") {
      log_mw[out_of_range] <- NA_real_
    }

    # Apply reasonable bounds
    log_mw[log_mw < 1] <- NA_real_ # < 10 Da
    log_mw[log_mw > 10] <- NA_real_ # > 10^10 Da

    # Convert to MW if requested
    if (log_output) {
      value <- as.numeric(log_mw)
    } else {
      value <- 10^log_mw
    }

    new_measure_tbl(location = location, value = value)
  })

  new_data[[output_col]] <- new_measure_list(mw_list)

  tibble::as_tibble(new_data)
}

#' @export
print.step_sec_broad_standard <- function(
  x,
  width = max(20, options()$width - 30),
  ...
) {
  title <- sprintf(
    "SEC broad standard calibration (%s, %s)",
    x$method,
    x$fit_type
  )

  if (x$trained && !is.null(x$calibration_diagnostics)) {
    diag <- x$calibration_diagnostics
    title <- sprintf(
      "%s, Mn err=%.1f%%, Mw err=%.1f%%",
      title,
      diag$mn_error_pct,
      diag$mw_error_pct
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
tidy.step_sec_broad_standard <- function(x, ...) {
  if (x$trained && !is.null(x$calibration_diagnostics)) {
    diag <- x$calibration_diagnostics
    tibble::tibble(
      method = x$method,
      fit_type = x$fit_type,
      known_mn = diag$known_mn,
      known_mw = diag$known_mw,
      known_dispersity = diag$known_dispersity,
      calculated_mn = diag$calculated_mn,
      calculated_mw = diag$calculated_mw,
      calculated_dispersity = diag$calculated_dispersity,
      mn_error_pct = diag$mn_error_pct,
      mw_error_pct = diag$mw_error_pct,
      dispersity_error_pct = diag$dispersity_error_pct,
      calibration_min = x$calibration_range[1],
      calibration_max = x$calibration_range[2],
      coefficients = list(x$calibration_coefficients),
      convergence = diag$convergence,
      output_col = x$output_col,
      id = x$id
    )
  } else {
    tibble::tibble(
      method = x$method,
      fit_type = x$fit_type,
      known_mn = x$known_mn,
      known_mw = x$known_mw,
      known_dispersity = if (!is.null(x$known_mw) && !is.null(x$known_mn)) {
        x$known_mw / x$known_mn
      } else {
        NA_real_
      },
      calculated_mn = NA_real_,
      calculated_mw = NA_real_,
      calculated_dispersity = NA_real_,
      mn_error_pct = NA_real_,
      mw_error_pct = NA_real_,
      dispersity_error_pct = NA_real_,
      calibration_min = NA_real_,
      calibration_max = NA_real_,
      coefficients = list(NULL),
      convergence = NA_integer_,
      output_col = x$output_col,
      id = x$id
    )
  }
}

#' @rdname required_pkgs.step_sec
#' @export
#' @keywords internal
required_pkgs.step_sec_broad_standard <- function(x, ...) {
  c("measure.sec", "measure")
}
