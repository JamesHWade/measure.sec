# ==============================================================================
# Calibration Save/Load Helper Functions
#
# Functions for saving and loading SEC calibration data for reuse
# ==============================================================================

#' Save SEC Calibration Parameters
#'
#' Extracts calibration parameters from a prepped recipe and saves them to a
#' file for later reuse. This allows calibrations to be established once and
#' applied to future analyses without re-fitting.
#'
#' @param prepped_recipe A prepped recipe containing a trained calibration step
#'   (e.g., `step_sec_conventional_cal`).
#' @param file Path to save the calibration file. File extension determines
#'   format: `.rds` for RDS format (default), `.yaml` or `.yml` for YAML format.
#' @param step_number Integer. Which step to extract calibration from. If NULL

#'   (default), finds the first calibration step automatically.
#' @param metadata Optional named list of additional metadata to include (e.g.,
#'   `instrument`, `column`, `analyst`, `date`).
#' @param overwrite Logical. If TRUE, overwrite existing file. Default FALSE.
#'
#' @return Invisibly returns the calibration object that was saved.
#'
#' @details
#' The saved calibration includes:
#' \itemize{
#'   \item Calibration fit coefficients and polynomial degree
#'   \item Calibration range (valid elution range)
#'   \item Fit diagnostics (R², RMSE, per-standard results)
#'   \item Original standards data
#'   \item Settings (fit_type, extrapolation, log_output)
#'   \item Timestamp and version information
#'   \item Optional user-provided metadata
#' }
#'
#' RDS format preserves all R objects exactly and is recommended for most uses.
#' YAML format is human-readable and useful for documentation or version control.
#'
#' @family sec-calibration
#' @seealso [load_sec_calibration()], [step_sec_conventional_cal()]
#' @export
#'
#' @examples
#' \dontrun{
#' # Create and prep a calibration recipe
#' ps_standards <- data.frame(
#'   retention = c(12.5, 13.2, 14.1, 15.0, 16.2, 17.5),
#'   log_mw = c(6.0, 5.5, 5.0, 4.5, 4.0, 3.5)
#' )
#'
#' rec <- recipe(~., data = sample_data) |>
#'   step_sec_conventional_cal(standards = ps_standards, fit_type = "cubic") |>
#'   prep()
#'
#' # Save calibration for later use
#' save_sec_calibration(
#'   rec,
#'   "ps_calibration.rds",
#'   metadata = list(
#'     column = "PLgel 5um Mixed-C",
#'     instrument = "Agilent 1260",
#'     analyst = "JW"
#'   )
#' )
#' }
save_sec_calibration <- function(
  prepped_recipe,
  file,
  step_number = NULL,
  metadata = NULL,
  overwrite = FALSE
) {
  # Validate inputs
  if (!inherits(prepped_recipe, "recipe")) {
    cli::cli_abort("{.arg prepped_recipe} must be a prepped recipe object.")
  }

  # Find calibration step
  cal_step <- .find_calibration_step(prepped_recipe, step_number)

  if (is.null(cal_step)) {
    cli::cli_abort(
      c(
        "No trained calibration step found in recipe.",
        "i" = "Recipe must contain a trained {.fn step_sec_conventional_cal}."
      )
    )
  }

  # Check if step is trained
  if (!cal_step$trained) {
    cli::cli_abort(
      c(
        "Calibration step is not trained.",
        "i" = "Use {.fn prep} on the recipe before saving calibration."
      )
    )
  }

  # Check file doesn't exist (unless overwrite)
  if (file.exists(file) && !overwrite) {
    cli::cli_abort(
      c(
        "File already exists: {.file {file}}",
        "i" = "Use {.arg overwrite = TRUE} to replace existing file."
      )
    )
  }

  # Determine format from extension
  ext <- tolower(tools::file_ext(file))
  if (!ext %in% c("rds", "yaml", "yml")) {
    cli::cli_warn(
      c(
        "Unrecognized file extension {.val {ext}}, using RDS format.",
        "i" = "Use {.val .rds} or {.val .yaml} extension."
      )
    )
    ext <- "rds"
  }

  # Build calibration object
  cal_object <- .build_calibration_object(cal_step, metadata)

  # Save based on format
  if (ext == "rds") {
    saveRDS(cal_object, file = file)
  } else {
    .save_calibration_yaml(cal_object, file)
  }

  cli::cli_alert_success(
    "Saved {.field {cal_step$fit_type}} calibration to {.file {file}}"
  )

  invisible(cal_object)
}


#' Load SEC Calibration Parameters
#'
#' Loads previously saved SEC calibration parameters for use in
#' `step_sec_conventional_cal()`. This enables reusing established calibrations
#' without refitting.
#'
#' @param file Path to the calibration file (`.rds` or `.yaml`/`.yml`).
#'
#' @return A `sec_calibration` object that can be passed to the `calibration`
#'   argument of [step_sec_conventional_cal()].
#'
#' @details
#' The loaded calibration can be used in two ways:
#'
#' 1. **Direct use**: Pass to `step_sec_conventional_cal(calibration = cal)`
#'    to skip fitting and use the pre-established calibration.
#'
#' 2. **Inspection**: Use `print()` or `summary()` to view calibration details.
#'
#' @family sec-calibration
#' @seealso [save_sec_calibration()], [step_sec_conventional_cal()]
#' @export
#'
#' @examples
#' \dontrun{
#' # Load a saved calibration
#' cal <- load_sec_calibration("ps_calibration.rds")
#' print(cal)
#'
#' # Use in a recipe (no fitting needed)
#' rec <- recipe(~., data = new_samples) |>
#'   step_sec_conventional_cal(calibration = cal) |>
#'   prep()
#' }
load_sec_calibration <- function(file) {
  if (!file.exists(file)) {
    cli::cli_abort("Calibration file not found: {.file {file}}")
  }

  # Determine format from extension
  ext <- tolower(tools::file_ext(file))

  if (ext == "rds") {
    cal_object <- readRDS(file)
  } else if (ext %in% c("yaml", "yml")) {
    cal_object <- .load_calibration_yaml(file)
  } else {
    # Try RDS first, fall back to YAML
    cal_object <- tryCatch(
      readRDS(file),
      error = function(e) .load_calibration_yaml(file)
    )
  }

  # Validate loaded object
  if (!inherits(cal_object, "sec_calibration")) {
    cli::cli_abort(
      c(
        "File does not contain a valid SEC calibration.",
        "i" = "Use {.fn save_sec_calibration} to create calibration files."
      )
    )
  }

  # Validate required components
  .validate_calibration_object(cal_object)

  cli::cli_alert_info(
    "Loaded {.field {cal_object$fit_type}} calibration (R\\u00b2 = {round(cal_object$diagnostics$r_squared, 4)})"
  )

  cal_object
}


# ==============================================================================
# Internal Helper Functions
# ==============================================================================

#' Find calibration step in recipe
#' @noRd
.find_calibration_step <- function(recipe, step_number = NULL) {
  steps <- recipe$steps

  if (!is.null(step_number)) {
    if (step_number < 1 || step_number > length(steps)) {
      cli::cli_abort(
        "Step number {step_number} is out of range (1-{length(steps)})."
      )
    }
    return(steps[[step_number]])
  }

  # Find first calibration step
  cal_classes <- c("step_sec_conventional_cal", "step_sec_universal_cal")

  for (step in steps) {
    if (any(class(step) %in% cal_classes)) {
      return(step)
    }
  }

  NULL
}


#' Build calibration object from step
#' @noRd
.build_calibration_object <- function(step, metadata = NULL) {
  # Extract coefficients from lm object
  fit <- step$calibration_fit
  coefs <- stats::coef(fit)

  # Get polynomial degree
  degree <- switch(
    step$fit_type,
    linear = 1,
    quadratic = 2,
    cubic = 3,
    fifth = 5
  )

  # Build calibration object
  cal <- structure(
    list(
      # Core calibration data
      fit_type = step$fit_type,
      degree = degree,
      coefficients = as.numeric(coefs),
      coef_names = names(coefs),
      calibration_range = step$calibration_range,

      # Model object (for RDS only - YAML will reconstruct)
      fit_object = fit,

      # Standards data
      standards = step$standards,

      # Diagnostics
      diagnostics = step$calibration_diagnostics,

      # Settings
      extrapolation = step$extrapolation,
      output_col = step$output_col,
      log_output = step$log_output,

      # Metadata
      created = Sys.time(),
      measure_sec_version = as.character(utils::packageVersion("measure.sec")),
      r_version = R.version.string,
      user_metadata = metadata
    ),
    class = c("sec_calibration", "list")
  )

  cal
}


#' Validate calibration object has required components
#' @noRd
.validate_calibration_object <- function(cal) {
  required <- c(
    "fit_type",
    "degree",
    "coefficients",
    "calibration_range",
    "standards"
  )

  missing <- setdiff(required, names(cal))

  if (length(missing) > 0) {
    cli::cli_abort(
      c(
        "Calibration object is missing required components.",
        "i" = "Missing: {.val {missing}}"
      )
    )
  }

  invisible(TRUE)
}


#' Save calibration to YAML format
#' @noRd
.save_calibration_yaml <- function(cal, file) {
  # Check if yaml package is available
  rlang::check_installed("yaml", reason = "to save calibration in YAML format")

  # Convert to YAML-friendly format
  yaml_obj <- list(
    sec_calibration = list(
      fit_type = cal$fit_type,
      degree = cal$degree,
      coefficients = as.list(cal$coefficients),
      coef_names = cal$coef_names,
      calibration_range = as.list(cal$calibration_range),
      extrapolation = cal$extrapolation,
      output_col = cal$output_col,
      log_output = cal$log_output,
      created = format(cal$created, "%Y-%m-%d %H:%M:%S"),
      measure_sec_version = cal$measure_sec_version,
      r_version = cal$r_version
    ),
    standards = as.list(cal$standards),
    diagnostics = list(
      r_squared = cal$diagnostics$r_squared,
      adj_r_squared = cal$diagnostics$adj_r_squared,
      rmse_log_mw = cal$diagnostics$rmse_log_mw,
      residual_std_error = cal$diagnostics$residual_std_error,
      max_abs_pct_deviation = cal$diagnostics$max_abs_pct_deviation,
      mean_abs_pct_deviation = cal$diagnostics$mean_abs_pct_deviation
    )
  )

  if (!is.null(cal$user_metadata)) {
    yaml_obj$user_metadata <- as.list(cal$user_metadata)
  }

  yaml::write_yaml(yaml_obj, file)
}


#' Load calibration from YAML format
#' @noRd
.load_calibration_yaml <- function(file) {
  rlang::check_installed(
    "yaml",
    reason = "to load calibration from YAML format"
  )

  yaml_obj <- yaml::read_yaml(file)

  if (is.null(yaml_obj$sec_calibration)) {
    cli::cli_abort("YAML file does not contain SEC calibration data.")
  }

  sec_cal <- yaml_obj$sec_calibration

  # Reconstruct standards data frame
  standards <- as.data.frame(yaml_obj$standards)

  # Reconstruct diagnostics
  diag <- yaml_obj$diagnostics

  # Reconstruct the lm fit from coefficients and standards
  fit_object <- .reconstruct_lm_fit(
    standards = standards,
    degree = sec_cal$degree,
    coefficients = unlist(sec_cal$coefficients)
  )

  # Build calibration object
  structure(
    list(
      fit_type = sec_cal$fit_type,
      degree = sec_cal$degree,
      coefficients = unlist(sec_cal$coefficients),
      coef_names = sec_cal$coef_names,
      calibration_range = unlist(sec_cal$calibration_range),
      fit_object = fit_object,
      standards = standards,
      diagnostics = list(
        r_squared = diag$r_squared,
        adj_r_squared = diag$adj_r_squared,
        rmse_log_mw = diag$rmse_log_mw,
        residual_std_error = diag$residual_std_error,
        max_abs_pct_deviation = diag$max_abs_pct_deviation,
        mean_abs_pct_deviation = diag$mean_abs_pct_deviation
      ),
      extrapolation = sec_cal$extrapolation,
      output_col = sec_cal$output_col,
      log_output = sec_cal$log_output,
      created = as.POSIXct(sec_cal$created),
      measure_sec_version = sec_cal$measure_sec_version,
      r_version = sec_cal$r_version,
      user_metadata = yaml_obj$user_metadata
    ),
    class = c("sec_calibration", "list")
  )
}


#' Reconstruct lm fit from coefficients (for YAML loading)
#' @noRd
.reconstruct_lm_fit <- function(standards, degree, coefficients) {
  # Get location column
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
  location_col <- names(standards)[names(standards) %in% location_cols][1]
  location <- standards[[location_col]]

  # Get log_mw
  if ("log_mw" %in% names(standards)) {
    log_mw <- standards$log_mw
  } else if ("mw" %in% names(standards)) {
    log_mw <- log10(standards$mw)
  } else {
    log_mw <- standards[[2]] # Fall back to second column
  }

  # Refit the model (this ensures predict() works correctly)
  cal_data <- data.frame(location = location, log_mw = log_mw)
  fit <- stats::lm(log_mw ~ stats::poly(location, degree), data = cal_data)

  # Replace coefficients with saved ones (should be identical, but ensures consistency)
  fit$coefficients <- coefficients

  fit
}


# ==============================================================================
# Print Methods
# ==============================================================================

#' @export
print.sec_calibration <- function(x, ...) {
  cli::cli_h1("SEC Calibration")

  cli::cli_text("{.field Fit type}: {x$fit_type} (degree {x$degree})")
  cli::cli_text(
    "{.field Calibration range}: {round(x$calibration_range[1], 2)} to {round(x$calibration_range[2], 2)}"
  )
  cli::cli_text("{.field Standards}: {nrow(x$standards)}")

  if (!is.null(x$diagnostics)) {
    cli::cli_h2("Fit Quality")
    cli::cli_text("{.field R\\u00b2}: {round(x$diagnostics$r_squared, 6)}")
    cli::cli_text(
      "{.field RMSE (log MW)}: {round(x$diagnostics$rmse_log_mw, 4)}"
    )
    cli::cli_text(
      "{.field Max % deviation}: {round(x$diagnostics$max_abs_pct_deviation, 2)}%"
    )
  }

  cli::cli_h2("Settings")
  cli::cli_text("{.field Extrapolation}: {x$extrapolation}")
  cli::cli_text("{.field Output column}: {x$output_col}")
  cli::cli_text("{.field Log output}: {x$log_output}")

  cli::cli_h2("Metadata")
  cli::cli_text("{.field Created}: {format(x$created, '%Y-%m-%d %H:%M:%S')}")
  cli::cli_text("{.field measure.sec version}: {x$measure_sec_version}")

  if (!is.null(x$user_metadata)) {
    cli::cli_h3("User Metadata")
    for (nm in names(x$user_metadata)) {
      cli::cli_text("{.field {nm}}: {x$user_metadata[[nm]]}")
    }
  }

  invisible(x)
}


#' @export
summary.sec_calibration <- function(object, ...) {
  print(object, ...)

  if (!is.null(object$diagnostics$standard_results)) {
    cli::cli_h2("Per-Standard Results")
    print(object$diagnostics$standard_results)
  }

  invisible(object)
}
