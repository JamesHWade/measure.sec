# ==============================================================================
# step_sec_mw_fractions
#
# Calculate molecular weight fractions above/below cutoffs
# ==============================================================================

#' Calculate Molecular Weight Fractions for SEC/GPC
#'
#' `step_sec_mw_fractions()` creates a *specification* of a recipe step that
#' calculates weight fractions above and below specified molecular weight cutoffs.
#'
#' @param recipe A recipe object.
#' @param measures An optional character vector of measure column names.
#' @param cutoffs Numeric vector of MW cutoff values. For each cutoff, the step
#'   calculates the weight fraction below and above that value.
#' @param calibration Calibration method for converting x-axis to log(MW).
#'   See [step_sec_mw_averages()] for details.
#' @param integration_range Optional numeric vector `c(min, max)` specifying
#'   the x-axis range for integration. If `NULL`, uses full range.
#' @param prefix Prefix for output column names. Default is `"frac_"`.
#' @param role Role for generated columns. Default is `"predictor"`.
#' @param trained Logical indicating if the step has been trained.
#' @param skip Logical. Should the step be skipped when baking?
#' @param id Unique step identifier.
#'
#' @return An updated recipe with the new step added.
#'
#' @details
#' For each cutoff value `C`, this step calculates:
#' - `frac_below_C`: Weight fraction with MW < C
#' - `frac_above_C`: Weight fraction with MW >= C
#'
#' These fractions sum to 1.0 and are useful for characterizing polymer
#' distributions. Common cutoffs include:
#' - 1000 Da for oligomer content
#' - 10000 Da for low MW fraction
#' - 100000 Da for high MW fraction
#'
#' @family sec-chromatography
#' @export
#'
#' @examples
#' \dontrun{
#' library(recipes)
#' library(measure)
#'
#' # Calculate fractions at multiple cutoffs
#' rec <- recipe(~., data = sec_data) |>
#'   step_measure_input_wide(starts_with("signal_")) |>
#'   step_sec_baseline() |>
#'   step_sec_mw_fractions(cutoffs = c(1000, 10000, 100000)) |>
#'   prep()
#' }
step_sec_mw_fractions <- function(
  recipe,
  measures = NULL,
  cutoffs = c(1000, 10000, 100000),
  calibration = NULL,
  integration_range = NULL,
  prefix = "frac_",
  role = "predictor",
  trained = FALSE,
  skip = FALSE,
  id = recipes::rand_id("sec_mw_fractions")
) {
  if (!is.numeric(cutoffs) || length(cutoffs) < 1) {
    cli::cli_abort(
      "{.arg cutoffs} must be a numeric vector with at least one value."
    )
  }

  if (any(cutoffs <= 0)) {
    cli::cli_abort("All {.arg cutoffs} must be positive values.")
  }

  if (!is.null(integration_range)) {
    if (!is.numeric(integration_range) || length(integration_range) != 2) {
      cli::cli_abort(
        "{.arg integration_range} must be a numeric vector of length 2."
      )
    }
  }

  recipes::add_step(
    recipe,
    step_sec_mw_fractions_new(
      measures = measures,
      cutoffs = sort(cutoffs),
      calibration = calibration,
      integration_range = integration_range,
      prefix = prefix,
      role = role,
      trained = trained,
      skip = skip,
      id = id
    )
  )
}

step_sec_mw_fractions_new <- function(
  measures,
  cutoffs,
  calibration,
  integration_range,
  prefix,
  role,
  trained,
  skip,
  id
) {
  recipes::step(
    subclass = "sec_mw_fractions",
    measures = measures,
    cutoffs = cutoffs,
    calibration = calibration,
    integration_range = integration_range,
    prefix = prefix,
    role = role,
    trained = trained,
    skip = skip,
    id = id
  )
}

#' @export
prep.step_sec_mw_fractions <- function(x, training, info = NULL, ...) {
  check_for_measure(training)

  if (is.null(x$measures)) {
    measure_cols <- find_measure_cols(training)
  } else {
    measure_cols <- x$measures
  }

  step_sec_mw_fractions_new(
    measures = measure_cols,
    cutoffs = x$cutoffs,
    calibration = x$calibration,
    integration_range = x$integration_range,
    prefix = x$prefix,
    role = x$role,
    trained = TRUE,
    skip = x$skip,
    id = x$id
  )
}

#' Calculate MW fractions from a single chromatogram
#' @noRd
.calc_mw_fractions <- function(
  location,
  value,
  cutoffs,
  calibration,
  integration_range
) {
  # Apply integration range
  if (!is.null(integration_range)) {
    idx <- location >= integration_range[1] & location <= integration_range[2]
    location <- location[idx]
    value <- value[idx]
  }

  n_cutoffs <- length(cutoffs)
  result <- numeric(n_cutoffs * 2)
  names(result) <- c(
    paste0("below_", cutoffs),
    paste0("above_", cutoffs)
  )

  if (length(location) < 2) {
    return(rep(NA_real_, length(result)))
  }

  # Convert x-axis to log10(MW) if calibration provided
  if (is.null(calibration)) {
    log_mw <- location
  } else if (is.numeric(calibration) && length(calibration) == 2) {
    log_mw <- calibration[1] * location + calibration[2]
  } else if (identical(calibration, "auto")) {
    log_mw <- seq(7, 2, length.out = length(location))
  } else {
    log_mw <- location
  }

  # Convert to MW
  mw <- 10^log_mw

  # Weight is proportional to signal
  w <- pmax(value, 0)
  total_w <- sum(w)

  if (total_w <= 0) {
    return(rep(NA_real_, length(result)))
  }

  # Calculate fractions for each cutoff
  for (i in seq_along(cutoffs)) {
    cutoff <- cutoffs[i]
    below_idx <- mw < cutoff
    frac_below <- sum(w[below_idx]) / total_w
    result[paste0("below_", cutoff)] <- frac_below
    result[paste0("above_", cutoff)] <- 1 - frac_below
  }

  result
}

#' @export
bake.step_sec_mw_fractions <- function(object, new_data, ...) {
  calibration <- object$calibration
  integration_range <- object$integration_range
  cutoffs <- object$cutoffs
  prefix <- object$prefix

  # Calculate MW fractions for each sample
  all_results <- purrr::map(new_data[[object$measures[1]]], function(m) {
    .calc_mw_fractions(
      m$location,
      m$value,
      cutoffs,
      calibration,
      integration_range
    )
  })

  # Convert to data frame
  result_df <- do.call(rbind, all_results)
  result_df <- tibble::as_tibble(result_df)

  # Add prefix to column names
  names(result_df) <- paste0(prefix, names(result_df))

  # Bind to original data
  new_data <- dplyr::bind_cols(new_data, result_df)

  tibble::as_tibble(new_data)
}

#' @export
print.step_sec_mw_fractions <- function(
  x,
  width = max(20, options()$width - 30),
  ...
) {
  cutoffs_str <- paste(format(x$cutoffs, scientific = FALSE), collapse = ", ")
  title <- paste0("SEC MW fractions at cutoffs: ", cutoffs_str)
  if (x$trained) {
    cat(title, " on <internal measurements>", sep = "")
  } else {
    cat(title)
  }

  cat("\n")

  invisible(x)
}

#' @rdname tidy.step_sec
#' @export
#' @keywords internal
tidy.step_sec_mw_fractions <- function(x, ...) {
  tibble::tibble(
    cutoffs = list(x$cutoffs),
    prefix = x$prefix,
    id = x$id
  )
}

#' @rdname required_pkgs.step_sec
#' @export
#' @keywords internal
required_pkgs.step_sec_mw_fractions <- function(x, ...) {
  c("measure.sec", "measure")
}
