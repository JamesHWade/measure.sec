# ==============================================================================
# step_sec_mw_distribution
#
# Generate molecular weight distribution curves
# ==============================================================================

#' Generate Molecular Weight Distribution Curve
#'
#' `step_sec_mw_distribution()` creates a *specification* of a recipe step
#' that generates molecular weight distribution curves from SEC/GPC data.
#'
#' @param recipe A recipe object.
#' @param measures An optional character vector of measure column names.
#' @param type Type of distribution to generate:
#'   - `"differential"` (default): dW/d(log M) differential distribution
#'   - `"cumulative"`: Cumulative weight fraction distribution
#'   - `"both"`: Generate both distributions
#' @param calibration Calibration method for converting x-axis to log(MW).
#'   See [step_sec_mw_averages()] for details.
#' @param n_points Number of points in the output distribution. Default is 100.
#'   If `NULL`, uses the original data resolution.
#' @param mw_range Optional numeric vector `c(min, max)` specifying the MW range
#'   for the output distribution. If `NULL`, uses the range from data.
#' @param normalize Logical. Should the differential distribution be normalized
#'   to integrate to 1? Default is `TRUE`.
#' @param role Role for generated columns. Default is `"predictor"`.
#' @param trained Logical indicating if the step has been trained.
#' @param skip Logical. Should the step be skipped when baking?
#' @param id Unique step identifier.
#'
#' @return An updated recipe with the new step added.
#'
#' @details
#' This step transforms the raw chromatogram into standard MW distribution
#' representations:
#'
#' **Differential Distribution (dW/d(log M)):**
#' The weight fraction per unit log(MW). This representation is preferred
#' because the area under the curve represents the weight fraction in that
#' MW range.
#'
#' **Cumulative Distribution:**
#' The cumulative weight fraction from low to high MW. Values range from 0 to 1.
#'
#' The output replaces the `.measures` column with the distribution data,
#' where `location` contains log10(MW) values and `value` contains the
#' distribution values.
#'
#' @family sec-chromatography
#' @export
#'
#' @examples
#' \dontrun{
#' library(recipes)
#' library(measure)
#'
#' # Generate differential MW distribution
#' rec <- recipe(~., data = sec_triple_detect) |>
#'   step_measure_input_wide(starts_with("signal_")) |>
#'   step_sec_baseline() |>
#'   step_sec_mw_distribution(type = "differential") |>
#'   prep()
#' }
step_sec_mw_distribution <- function(
  recipe,
  measures = NULL,
  type = c("differential", "cumulative", "both"),
  calibration = NULL,
  n_points = 100L,
  mw_range = NULL,
  normalize = TRUE,
  role = NA,
  trained = FALSE,
  skip = FALSE,
  id = recipes::rand_id("sec_mw_distribution")
) {
  type <- match.arg(type)

  if (!is.null(n_points)) {
    if (!is.numeric(n_points) || length(n_points) != 1 || n_points < 10) {
      cli::cli_abort(
        "{.arg n_points} must be a positive integer >= 10 or NULL."
      )
    }
    n_points <- as.integer(n_points)
  }

  if (!is.null(mw_range)) {
    if (!is.numeric(mw_range) || length(mw_range) != 2) {
      cli::cli_abort("{.arg mw_range} must be a numeric vector of length 2.")
    }
    if (mw_range[1] >= mw_range[2]) {
      cli::cli_abort("{.arg mw_range} must have min < max.")
    }
  }

  recipes::add_step(
    recipe,
    step_sec_mw_distribution_new(
      measures = measures,
      type = type,
      calibration = calibration,
      n_points = n_points,
      mw_range = mw_range,
      normalize = normalize,
      role = role,
      trained = trained,
      skip = skip,
      id = id
    )
  )
}

step_sec_mw_distribution_new <- function(
  measures,
  type,
  calibration,
  n_points,
  mw_range,
  normalize,
  role,
  trained,
  skip,
  id
) {
  recipes::step(
    subclass = "sec_mw_distribution",
    measures = measures,
    type = type,
    calibration = calibration,
    n_points = n_points,
    mw_range = mw_range,
    normalize = normalize,
    role = role,
    trained = trained,
    skip = skip,
    id = id
  )
}

#' @export
prep.step_sec_mw_distribution <- function(x, training, info = NULL, ...) {
  check_for_measure(training)

  if (is.null(x$measures)) {
    measure_cols <- find_measure_cols(training)
  } else {
    measure_cols <- x$measures
  }

  step_sec_mw_distribution_new(
    measures = measure_cols,
    type = x$type,
    calibration = x$calibration,
    n_points = x$n_points,
    mw_range = x$mw_range,
    normalize = x$normalize,
    role = x$role,
    trained = TRUE,
    skip = x$skip,
    id = x$id
  )
}

#' Generate MW distribution from a single chromatogram
#' @noRd
.generate_mw_distribution <- function(
  location,
  value,
  type,
  calibration,
  n_points,
  mw_range,
  normalize
) {
  if (length(location) < 2) {
    return(list(
      differential = new_measure_tbl(numeric(0), numeric(0)),
      cumulative = new_measure_tbl(numeric(0), numeric(0))
    ))
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

  # Weight is proportional to signal
  w <- pmax(value, 0)

  # Sort by log_mw (ascending order for cumulative)
  ord <- order(log_mw)
  log_mw <- log_mw[ord]
  w <- w[ord]

  # Determine MW range for output
  if (is.null(mw_range)) {
    out_range <- range(log_mw)
  } else {
    out_range <- log10(mw_range)
  }

  # Create output grid
  if (is.null(n_points)) {
    out_log_mw <- log_mw
  } else {
    out_log_mw <- seq(out_range[1], out_range[2], length.out = n_points)
  }

  result <- list()

  # Calculate differential distribution: dW/d(log M)
  if (type %in% c("differential", "both")) {
    # Interpolate weights to output grid
    if (length(log_mw) >= 2) {
      interp_w <- stats::approx(log_mw, w, xout = out_log_mw, rule = 2)$y
    } else {
      interp_w <- rep(w[1], length(out_log_mw))
    }

    # For differential distribution, we want dW/d(log M)
    # The detector signal is already proportional to dW/dV (volume)
    # We need to convert to dW/d(log M)

    # If calibration is linear: log M = a*V + b
    # Then d(log M)/dV = a (constant)
    # So dW/d(log M) = (1/a) * dW/dV

    # For now, we use the interpolated weights directly as the differential
    # This assumes uniform spacing in log M space after interpolation
    diff_dist <- interp_w

    if (normalize && sum(diff_dist) > 0) {
      # Normalize so integral = 1
      d_log_mw <- if (length(out_log_mw) > 1) {
        mean(diff(out_log_mw))
      } else {
        1
      }
      diff_dist <- diff_dist / (sum(diff_dist) * d_log_mw)
    }

    result$differential <- new_measure_tbl(out_log_mw, diff_dist)
  }

  # Calculate cumulative distribution
  if (type %in% c("cumulative", "both")) {
    # Interpolate weights to output grid
    if (length(log_mw) >= 2) {
      interp_w <- stats::approx(log_mw, w, xout = out_log_mw, rule = 2)$y
    } else {
      interp_w <- rep(w[1], length(out_log_mw))
    }

    # Cumulative sum normalized to \[0, 1\]
    cum_dist <- cumsum(interp_w)
    if (max(cum_dist) > 0) {
      cum_dist <- cum_dist / max(cum_dist)
    }

    result$cumulative <- new_measure_tbl(out_log_mw, cum_dist)
  }

  result
}

#' @export
bake.step_sec_mw_distribution <- function(object, new_data, ...) {
  type <- object$type
  calibration <- object$calibration
  n_points <- object$n_points
  mw_range <- object$mw_range
  normalize <- object$normalize
  measure_col <- object$measures[1]

  # Generate distributions for each sample
  all_results <- purrr::map(new_data[[measure_col]], function(m) {
    .generate_mw_distribution(
      m$location,
      m$value,
      type,
      calibration,
      n_points,
      mw_range,
      normalize
    )
  })

  # Replace measure column with distribution data
  if (type == "differential") {
    new_data[[measure_col]] <- new_measure_list(
      purrr::map(all_results, ~ .x$differential)
    )
  } else if (type == "cumulative") {
    new_data[[measure_col]] <- new_measure_list(
      purrr::map(all_results, ~ .x$cumulative)
    )
  } else {
    # Both - create two columns
    new_data[[paste0(measure_col, "_differential")]] <- new_measure_list(
      purrr::map(all_results, ~ .x$differential)
    )
    new_data[[paste0(measure_col, "_cumulative")]] <- new_measure_list(
      purrr::map(all_results, ~ .x$cumulative)
    )
    # Keep original column as differential
    new_data[[measure_col]] <- new_measure_list(
      purrr::map(all_results, ~ .x$differential)
    )
  }

  tibble::as_tibble(new_data)
}

#' @export
print.step_sec_mw_distribution <- function(
  x,
  width = max(20, options()$width - 30),
  ...
) {
  title <- paste0("SEC MW distribution (", x$type, ")")
  if (!is.null(x$n_points)) {
    title <- paste0(title, ", ", x$n_points, " points")
  }
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
tidy.step_sec_mw_distribution <- function(x, ...) {
  if (is.null(x$n_points)) {
    n_pts <- NA_integer_
  } else {
    n_pts <- x$n_points
  }
  tibble::tibble(
    type = x$type,
    n_points = n_pts,
    normalize = x$normalize,
    id = x$id
  )
}

#' @rdname required_pkgs.step_sec
#' @export
#' @keywords internal
required_pkgs.step_sec_mw_distribution <- function(x, ...) {
  c("measure.sec", "measure")
}
