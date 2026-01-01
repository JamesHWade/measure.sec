# ==============================================================================
# step_sec_viscometer
#
# Differential viscometer processing for SEC
# ==============================================================================

#' Differential Viscometer Processing for SEC
#'
#' `step_sec_viscometer()` creates a *specification* of a recipe step that
#' processes differential viscometer signals to calculate specific viscosity
#' at each elution point.
#'
#' @param recipe A recipe object.
#' @param dp_col Name of the differential pressure (DP) measure column.
#' @param ip_col Name of the inlet pressure (IP) measure column. If `NULL`,
#'   assumes a single-capillary viscometer where DP is directly proportional
#'   to specific viscosity.
#' @param output_col Name for the specific viscosity output column.
#'   Default is `"specific_visc"`.
#' @param viscometer_constant Instrument calibration constant. Default is 1.0.
#'   Obtain from viscosity standard calibration.
#' @param min_signal Minimum signal threshold (as fraction of max) below which
#'   viscosity is set to NA. Default is 0.01.
#' @param role Role for generated columns.
#' @param trained Logical indicating if the step has been trained.
#' @param skip Logical. Should the step be skipped when baking?
#' @param id Unique step identifier.
#'
#' @return An updated recipe with the new step added, containing specific
#'   viscosity at each elution point.
#'
#' @details
#' Differential viscometers measure the pressure difference between a sample
#' and reference capillary to determine solution viscosity. The specific
#' viscosity is calculated from the differential pressure (DP) and inlet
#' pressure (IP):
#'
#' \deqn{\eta_{sp} = \frac{4 \cdot DP}{IP - 2 \cdot DP}}
#'
#' For single-capillary viscometers or when IP is not available:
#' \deqn{\eta_{sp} = K \cdot DP}
#'
#' where K is a calibration constant.
#'
#' **Viscometry in SEC:**
#' \itemize{
#'   \item Provides specific viscosity at each MW slice
#'   \item Combined with concentration gives intrinsic viscosity \[eta\]
#'   \item Used for universal calibration: log(\[eta\] * M) vs retention
#'   \item Essential for branching analysis (g' = \[eta\]_branched / \[eta\]_linear)
#' }
#'
#' @note
#' For intrinsic viscosity calculation, use `step_sec_intrinsic_visc()` after
#' this step.
#'
#' @family sec-detectors
#' @seealso [step_sec_intrinsic_visc()] for intrinsic viscosity calculation
#' @export
#'
#' @examples
#' \dontrun{
#' library(recipes)
#' library(measure)
#'
#' # Process differential viscometer data
#' rec <- recipe(~., data = sec_visc_data) |>
#'   step_measure_input_long(dp_signal, location = vars(elution_time), col_name = "dp") |>
#'   step_measure_input_long(ip_signal, location = vars(elution_time), col_name = "ip") |>
#'   step_sec_baseline() |>
#'   step_sec_viscometer(dp_col = "dp", ip_col = "ip") |>
#'   prep()
#' }
step_sec_viscometer <- function(
    recipe,
    dp_col = NULL,
    ip_col = NULL,
    output_col = "specific_visc",
    viscometer_constant = 1.0,
    min_signal = 0.01,
    role = NA,
    trained = FALSE,
    skip = FALSE,
    id = recipes::rand_id("sec_viscometer")
) {

  if (!is.numeric(viscometer_constant) || viscometer_constant <= 0) {
    cli::cli_abort("{.arg viscometer_constant} must be a positive number.")
  }

  recipes::add_step(
    recipe,
    step_sec_viscometer_new(
      dp_col = dp_col,
      ip_col = ip_col,
      output_col = output_col,
      viscometer_constant = viscometer_constant,
      min_signal = min_signal,
      role = role,
      trained = trained,
      skip = skip,
      id = id
    )
  )
}

step_sec_viscometer_new <- function(
    dp_col,
    ip_col,
    output_col,
    viscometer_constant,
    min_signal,
    role,
    trained,
    skip,
    id
) {
  recipes::step(
    subclass = "sec_viscometer",
    dp_col = dp_col,
    ip_col = ip_col,
    output_col = output_col,
    viscometer_constant = viscometer_constant,
    min_signal = min_signal,
    role = role,
    trained = trained,
    skip = skip,
    id = id
  )
}

#' @export
prep.step_sec_viscometer <- function(x, training, info = NULL, ...) {
  check_for_measure(training)

  measure_cols <- find_measure_cols(training)

  # Find DP column if not specified
  if (is.null(x$dp_col)) {
    dp_cols <- measure_cols[grepl("dp|visc|pressure", measure_cols, ignore.case = TRUE)]
    if (length(dp_cols) == 0) {
      cli::cli_abort(
        "No viscometer column found. Specify {.arg dp_col} explicitly."
      )
    }
    dp_col <- dp_cols[1]
  } else {
    dp_col <- x$dp_col
    if (!dp_col %in% measure_cols) {
      cli::cli_abort("DP column {.val {dp_col}} not found in measure columns.")
    }
  }

  # Validate IP column if specified
  ip_col <- x$ip_col
  if (!is.null(ip_col) && !ip_col %in% measure_cols) {
    cli::cli_abort("IP column {.val {ip_col}} not found in measure columns.")
  }

  step_sec_viscometer_new(
    dp_col = dp_col,
    ip_col = ip_col,
    output_col = x$output_col,
    viscometer_constant = x$viscometer_constant,
    min_signal = x$min_signal,
    role = x$role,
    trained = TRUE,
    skip = x$skip,
    id = x$id
  )
}

#' @export
bake.step_sec_viscometer <- function(object, new_data, ...) {
  dp_col <- object$dp_col
  ip_col <- object$ip_col
  output_col <- object$output_col
  viscometer_constant <- object$viscometer_constant
  min_signal <- object$min_signal

  # Calculate specific viscosity for each sample
  if (is.null(ip_col)) {
    # Single capillary: eta_sp = K * DP
    visc_list <- purrr::map(new_data[[dp_col]], function(dp_m) {
      dp_val <- dp_m$value
      location <- dp_m$location

      # Threshold based on DP signal
      max_dp <- max(abs(dp_val), na.rm = TRUE)
      threshold <- min_signal * max_dp

      eta_sp <- rep(NA_real_, length(dp_val))
      valid <- abs(dp_val) > threshold & !is.na(dp_val)

      if (any(valid)) {
        eta_sp[valid] <- viscometer_constant * dp_val[valid]
      }

      new_measure_tbl(location = location, value = eta_sp)
    })
  } else {
    # Differential viscometer: eta_sp = 4*DP / (IP - 2*DP)
    visc_list <- purrr::map2(
      new_data[[dp_col]],
      new_data[[ip_col]],
      function(dp_m, ip_m) {
        dp_val <- dp_m$value
        ip_val <- ip_m$value
        location <- dp_m$location

        # Threshold based on DP signal
        max_dp <- max(abs(dp_val), na.rm = TRUE)
        threshold <- min_signal * max_dp

        eta_sp <- rep(NA_real_, length(dp_val))
        denominator <- ip_val - 2 * dp_val

        # Valid where signal is above threshold and denominator is positive
        valid <- abs(dp_val) > threshold & !is.na(dp_val) & !is.na(ip_val) &
          denominator > 0

        if (any(valid)) {
          eta_sp[valid] <- viscometer_constant * 4 * dp_val[valid] / denominator[valid]
        }

        new_measure_tbl(location = location, value = eta_sp)
      }
    )
  }

  new_data[[output_col]] <- new_measure_list(visc_list)

  tibble::as_tibble(new_data)
}

#' @export
print.step_sec_viscometer <- function(
    x,
    width = max(20, options()$width - 30),
    ...
) {
  title <- "SEC viscometer processing"
  if (x$trained) {
    if (is.null(x$ip_col)) {
      cat(title, " (", x$dp_col, " -> ", x$output_col, ")", sep = "")
    } else {
      cat(title, " (", x$dp_col, ", ", x$ip_col, " -> ", x$output_col, ")", sep = "")
    }
  } else {
    cat(title)
  }
  cat("\n")
  invisible(x)
}

#' @rdname tidy.step_sec
#' @export
#' @keywords internal
tidy.step_sec_viscometer <- function(x, ...) {
  tibble::tibble(
    dp_col = x$dp_col %||% NA_character_,
    ip_col = x$ip_col %||% NA_character_,
    output_col = x$output_col,
    viscometer_constant = x$viscometer_constant,
    id = x$id
  )
}

#' @rdname required_pkgs.step_sec
#' @export
#' @keywords internal
required_pkgs.step_sec_viscometer <- function(x, ...) {
  c("measure.sec", "measure")
}
