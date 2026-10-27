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
#' @param mw_col Optional per-slice molecular-weight measure column. A matching
#'   curve produced by [step_sec_conventional_cal()] is detected automatically.
#' @param mw_scale Scale of `mw_col`: `"log10"` or `"absolute"`. Required for
#'   untagged curves; see [step_sec_mw_averages()].
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
#' Slice weights use trapezoidal integration on the elution axis. The cutoff is
#' applied to the sampled molecular-weight values, so accuracy near a cutoff
#' depends on the sampling resolution. These fractions sum to 1.0 and are useful for characterizing polymer
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
#' rec <- recipe(~., data = sec_triple_detect) |>
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
  id = recipes::rand_id("sec_mw_fractions"),
  mw_col = NULL,
  mw_scale = NULL
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
      mw_col = mw_col,
      mw_scale = mw_scale,
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
  mw_col,
  mw_scale,
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
    mw_col = mw_col,
    mw_scale = mw_scale,
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

  resolved <- .sec_resolve_mw_input(
    training,
    x$measures,
    x$mw_col,
    x$mw_scale,
    x$calibration
  )

  step_sec_mw_fractions_new(
    measures = resolved$measures,
    cutoffs = x$cutoffs,
    calibration = x$calibration,
    mw_col = resolved$mw_col,
    mw_scale = resolved$mw_scale,
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
  integration_range,
  mw_measure = NULL,
  mw_scale = NULL
) {
  slices <- .sec_mw_slices(
    list(location = location, value = value),
    mw_measure,
    mw_scale,
    calibration,
    integration_range
  )

  n_cutoffs <- length(cutoffs)
  result <- numeric(n_cutoffs * 2)
  names(result) <- c(
    paste0("below_", cutoffs),
    paste0("above_", cutoffs)
  )

  if (length(slices$location) < 2) {
    result[] <- NA_real_
    return(result)
  }
  mw <- slices$mw
  w <- slices$weight
  valid <- is.finite(mw) & mw > 0 & w > 0
  mw <- mw[valid]
  w <- w[valid]
  total_w <- sum(w)

  if (total_w <= 0) {
    result[] <- NA_real_
    return(result)
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
  mw_curves <- if (is.null(object$mw_col)) {
    rep(list(NULL), nrow(new_data))
  } else {
    new_data[[object$mw_col]]
  }
  all_results <- purrr::map2(
    new_data[[object$measures[1]]],
    mw_curves,
    function(m, mw_curve) {
      .calc_mw_fractions(
        m$location,
        m$value,
        cutoffs,
        calibration,
        integration_range,
        mw_measure = mw_curve,
        mw_scale = object$mw_scale
      )
    }
  )

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
    mw_col = x$mw_col %||% NA_character_,
    mw_scale = x$mw_scale %||% NA_character_,
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
