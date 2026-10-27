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
#' @param mw_col Optional per-slice molecular-weight measure column. A matching
#'   curve produced by [step_sec_conventional_cal()] is detected automatically.
#' @param mw_scale Scale of `mw_col`: `"log10"` or `"absolute"`. Required for
#'   untagged curves; see [step_sec_mw_averages()].
#' @param integration_range Optional limits on the original elution axis.
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
#' Weight is integrated on the original elution axis before transforming to
#' log(MW), so a nonlinear calibration is handled without treating the signal
#' as a density on the molecular-weight axis. The selected MW curve must be
#' strictly monotonic and valid over the integration range. When `mw_range`
#' restricts the plotted range, it displays a portion of the full distribution;
#' it does not renormalize that portion to 1.
#'
#' The selected measure column is replaced by the requested distribution.
#' With `type = "both"`, two additional columns, `<measure>_differential` and
#' `<measure>_cumulative`, are added and the selected column is replaced by the
#' differential curve. Their `location` values are log10(MW), and their
#' `value` entries are the corresponding distribution values.
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
  id = recipes::rand_id("sec_mw_distribution"),
  mw_col = NULL,
  mw_scale = NULL,
  integration_range = NULL
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
    if (
      !is.numeric(mw_range) ||
        length(mw_range) != 2 ||
        !all(is.finite(mw_range)) ||
        any(mw_range <= 0)
    ) {
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
      mw_col = mw_col,
      mw_scale = mw_scale,
      integration_range = integration_range,
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
  mw_col,
  mw_scale,
  integration_range,
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
    mw_col = mw_col,
    mw_scale = mw_scale,
    integration_range = integration_range,
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

  resolved <- .sec_resolve_mw_input(
    training,
    x$measures,
    x$mw_col,
    x$mw_scale,
    x$calibration
  )

  step_sec_mw_distribution_new(
    measures = resolved$measures,
    type = x$type,
    calibration = x$calibration,
    mw_col = resolved$mw_col,
    mw_scale = resolved$mw_scale,
    integration_range = x$integration_range,
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
  normalize,
  mw_measure = NULL,
  mw_scale = NULL,
  integration_range = NULL
) {
  slices <- .sec_mw_slices(
    list(location = location, value = value),
    mw_measure,
    mw_scale,
    calibration,
    integration_range
  )
  if (length(slices$location) < 2) {
    return(
      list(
        differential = new_measure_tbl(numeric(0), numeric(0)),
        cumulative = new_measure_tbl(numeric(0), numeric(0))
      )
    )
  }

  if (!all(is.finite(slices$mw)) || any(slices$mw <= 0)) {
    cli::cli_abort(
      "The selected distribution requires a valid molecular weight at every slice."
    )
  }
  log_mw <- log10(slices$mw)
  if (!(all(diff(log_mw) > 0) || all(diff(log_mw) < 0))) {
    cli::cli_abort(
      "A molecular-weight distribution requires a strictly monotonic calibration over the selected range."
    )
  }
  ord <- order(log_mw)
  log_mw <- log_mw[ord]
  elution <- slices$location[ord]
  signal <- slices$signal[ord]
  segment_mass <- abs(diff(elution)) *
    (utils::head(signal, -1) + utils::tail(signal, -1)) /
    2
  cumulative_mass <- c(0, cumsum(segment_mass))
  total_mass <- utils::tail(cumulative_mass, 1)

  # Determine MW range for output
  if (is.null(mw_range)) {
    out_range <- range(log_mw)
  } else {
    out_range <- log10(mw_range)
  }

  # Create output grid
  if (is.null(n_points)) {
    out_log_mw <- if (is.null(mw_range)) {
      log_mw
    } else {
      sort(unique(c(
        out_range,
        log_mw[log_mw > out_range[1] & log_mw < out_range[2]]
      )))
    }
  } else {
    out_log_mw <- seq(out_range[1], out_range[2], length.out = n_points)
  }

  result <- list()
  grid_mass <- stats::approx(
    log_mw,
    cumulative_mass,
    xout = out_log_mw,
    rule = 2
  )$y
  if (type %in% c("differential", "both")) {
    slope <- diff(grid_mass) / diff(out_log_mw)
    density <- c(
      slope[1],
      if (length(slope) > 1) {
        (utils::head(slope, -1) + utils::tail(slope, -1)) / 2
      },
      utils::tail(slope, 1)
    )
    area <- sum(
      diff(out_log_mw) *
        (utils::head(density, -1) + utils::tail(density, -1)) /
        2
    )
    covered_mass <- utils::tail(grid_mass, 1) - grid_mass[1]
    if (area > 0) {
      density <- density * covered_mass / area
    }
    if (normalize && total_mass > 0) {
      density <- density / total_mass
    }
    result$differential <- new_measure_tbl(out_log_mw, density)
  }
  if (type %in% c("cumulative", "both")) {
    cum_dist <- if (total_mass > 0) {
      grid_mass / total_mass
    } else {
      rep(0, length(grid_mass))
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
  mw_curves <- if (is.null(object$mw_col)) {
    rep(list(NULL), nrow(new_data))
  } else {
    new_data[[object$mw_col]]
  }
  all_results <- purrr::map2(
    new_data[[measure_col]],
    mw_curves,
    function(m, mw_curve) {
      .generate_mw_distribution(
        m$location,
        m$value,
        type,
        calibration,
        n_points,
        mw_range,
        normalize,
        mw_measure = mw_curve,
        mw_scale = object$mw_scale,
        integration_range = object$integration_range
      )
    }
  )

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
    mw_col = x$mw_col %||% NA_character_,
    mw_scale = x$mw_scale %||% NA_character_,
    id = x$id
  )
}

#' @rdname required_pkgs.step_sec
#' @export
#' @keywords internal
required_pkgs.step_sec_mw_distribution <- function(x, ...) {
  c("measure.sec", "measure")
}
