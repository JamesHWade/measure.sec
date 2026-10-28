# ==============================================================================
# step_sec_uv_ri_ratio
#
# Calculate UV/RI ratio for composition analysis
# ==============================================================================

#' Calculate UV/RI Ratio for Composition Analysis
#'
#' `step_sec_uv_ri_ratio()` creates a *specification* of a recipe step that
#' calculates the ratio of UV to RI detector signals at each elution point.
#' This ratio is useful for detecting compositional heterogeneity in copolymers
#' and conjugates.
#'
#' @param recipe A recipe object.
#' @param uv_col Name of the UV detector measure column.
#' @param ri_col Name of the RI detector measure column.
#' @param output_col Name for the output ratio column. Default is `"uv_ri_ratio"`.
#' @param min_signal Minimum signal threshold (as fraction of max) below which
#'   ratio is set to NA. Default is 0.01 (1%). Prevents noisy ratios in baseline.
#' @param smooth Logical. Apply smoothing to the ratio? Default is `TRUE`.
#' @param smooth_window Window size for smoothing (number of points). Default is 5.
#' @param role Role for generated columns.
#' @param trained Logical indicating if the step has been trained.
#' @param skip Logical. Should the step be skipped when baking?
#' @param id Unique step identifier.
#'
#' @return An updated recipe with the new step added.
#'
#' @details
#' The UV/RI ratio provides information about chemical composition across the
#' molecular weight distribution:
#'
#' \deqn{Ratio(t) = \frac{UV(t)}{RI(t)} = \frac{\varepsilon \cdot c(t)}{(dn/dc) \cdot c(t) \cdot K}}
#'
#' Since concentration cancels out, the ratio reflects the relative detector
#' response factors, which depend on chemical composition.
#'
#' **Applications:**
#' \itemize{
#'   \item Copolymer composition drift with molecular weight
#'   \item Block copolymer characterization
#'   \item PEGylation analysis (protein-PEG conjugates)
#'   \item Detection of chemical heterogeneity
#'   \item End-group analysis with UV labels
#' }
#'
#' **Interpretation:**
#' \itemize{
#'   \item Constant ratio: Uniform composition across MW
#'   \item Increasing ratio with MW: More chromophore in higher MW species
#'   \item Decreasing ratio with MW: Less chromophore in higher MW species
#' }
#'
#' @note
#' Both UV and RI signals should be baseline-corrected and properly aligned
#' (using `step_sec_detector_delay()`) before calculating the ratio.
#'
#' @family sec-composition
#' @export
#'
#' @examples
#' \dontrun{
#' library(recipes)
#' library(measure)
#'
#' # Calculate UV/RI ratio for copolymer
#' rec <- recipe(~., data = sec_copolymer) |>
#'   step_measure_input_long(ri_signal, location = vars(elution_time), col_name = "ri") |>
#'   step_measure_input_long(uv_signal, location = vars(elution_time), col_name = "uv") |>
#'   step_sec_detector_delay(reference = "ri", delay_volumes = c(uv = 0.05)) |>
#'   step_sec_baseline() |>
#'   step_sec_uv_ri_ratio(uv_col = "uv", ri_col = "ri") |>
#'   prep()
#' }
step_sec_uv_ri_ratio <- function(
  recipe,
  uv_col = NULL,
  ri_col = NULL,
  output_col = "uv_ri_ratio",
  min_signal = 0.01,
  smooth = TRUE,
  smooth_window = 5,
  role = NA,
  trained = FALSE,
  skip = FALSE,
  id = recipes::rand_id("sec_uv_ri_ratio")
) {
  # Validate inputs

  if (!is.numeric(min_signal) || min_signal < 0 || min_signal > 1) {
    cli::cli_abort("{.arg min_signal} must be between 0 and 1.")
  }

  if (!is.numeric(smooth_window) || smooth_window < 1) {
    cli::cli_abort("{.arg smooth_window} must be a positive integer.")
  }

  recipes::add_step(
    recipe,
    step_sec_uv_ri_ratio_new(
      uv_col = uv_col,
      ri_col = ri_col,
      output_col = output_col,
      min_signal = min_signal,
      smooth = smooth,
      smooth_window = as.integer(smooth_window),
      role = role,
      trained = trained,
      skip = skip,
      id = id
    )
  )
}

step_sec_uv_ri_ratio_new <- function(
  uv_col,
  ri_col,
  output_col,
  min_signal,
  smooth,
  smooth_window,
  role,
  trained,
  skip,
  id
) {
  recipes::step(
    subclass = "sec_uv_ri_ratio",
    uv_col = uv_col,
    ri_col = ri_col,
    output_col = output_col,
    min_signal = min_signal,
    smooth = smooth,
    smooth_window = smooth_window,
    role = role,
    trained = trained,
    skip = skip,
    id = id
  )
}

#' @export
prep.step_sec_uv_ri_ratio <- function(x, training, info = NULL, ...) {
  check_for_measure(training)

  measure_cols <- find_measure_cols(training)

  # Find UV column if not specified
  if (is.null(x$uv_col)) {
    uv_cols <- measure_cols[grepl("uv", measure_cols, ignore.case = TRUE)]
    if (length(uv_cols) == 0) {
      cli::cli_abort("No UV column found. Specify {.arg uv_col} explicitly.")
    }
    uv_col <- uv_cols[1]
  } else {
    uv_col <- x$uv_col
    if (!uv_col %in% measure_cols) {
      cli::cli_abort("UV column {.val {uv_col}} not found in measure columns.")
    }
  }

  # Find RI column if not specified
  if (is.null(x$ri_col)) {
    ri_cols <- measure_cols[grepl("ri", measure_cols, ignore.case = TRUE)]
    if (length(ri_cols) == 0) {
      cli::cli_abort("No RI column found. Specify {.arg ri_col} explicitly.")
    }
    ri_col <- ri_cols[1]
  } else {
    ri_col <- x$ri_col
    if (!ri_col %in% measure_cols) {
      cli::cli_abort("RI column {.val {ri_col}} not found in measure columns.")
    }
  }

  step_sec_uv_ri_ratio_new(
    uv_col = uv_col,
    ri_col = ri_col,
    output_col = x$output_col,
    min_signal = x$min_signal,
    smooth = x$smooth,
    smooth_window = x$smooth_window,
    role = x$role,
    trained = TRUE,
    skip = x$skip,
    id = x$id
  )
}

#' Simple moving average smoother
#' @noRd
.moving_average <- function(x, window) {
  n <- length(x)
  if (n < window) {
    return(x)
  }

  result <- x
  half_window <- floor(window / 2)

  for (i in seq_len(n)) {
    start_idx <- max(1, i - half_window)
    end_idx <- min(n, i + half_window)
    # Only average non-NA values
    vals <- x[start_idx:end_idx]
    vals <- vals[!is.na(vals)]
    if (length(vals) > 0) {
      result[i] <- mean(vals)
    }
  }

  result
}

#' @export
bake.step_sec_uv_ri_ratio <- function(object, new_data, ...) {
  uv_col <- object$uv_col
  ri_col <- object$ri_col
  output_col <- object$output_col
  min_signal <- object$min_signal
  smooth <- object$smooth
  smooth_window <- object$smooth_window

  # Calculate ratio for each sample
  ratio_list <- purrr::map2(
    new_data[[uv_col]],
    new_data[[ri_col]],
    function(uv_m, ri_m) {
      uv_val <- uv_m$value
      ri_val <- ri_m$value
      location <- uv_m$location

      # Determine signal threshold
      max_ri <- max(abs(ri_val), na.rm = TRUE)
      threshold <- min_signal * max_ri

      # Calculate ratio where signal is above threshold
      ratio <- rep(NA_real_, length(uv_val))
      valid <- abs(ri_val) > threshold & !is.na(uv_val) & !is.na(ri_val)
      ratio[valid] <- uv_val[valid] / ri_val[valid]

      # Apply smoothing if requested
      if (smooth && sum(valid) > smooth_window) {
        ratio <- .moving_average(ratio, smooth_window)
      }

      # Return as measure object
      new_measure_tbl(location = location, value = ratio)
    }
  )

  new_data[[output_col]] <- new_measure_list(ratio_list)

  tibble::as_tibble(new_data)
}

#' @export
print.step_sec_uv_ri_ratio <- function(
  x,
  width = max(20, options()$width - 30),
  ...
) {
  title <- "SEC UV/RI ratio"
  if (x$trained) {
    cat(
      title,
      " (",
      x$uv_col,
      "/",
      x$ri_col,
      " -> ",
      x$output_col,
      ")",
      sep = ""
    )
  } else {
    cat(title)
  }
  cat("\n")
  invisible(x)
}

#' @rdname tidy.step_sec
#' @export
#' @keywords internal
tidy.step_sec_uv_ri_ratio <- function(x, ...) {
  tibble::tibble(
    uv_col = x$uv_col %||% NA_character_,
    ri_col = x$ri_col %||% NA_character_,
    output_col = x$output_col,
    min_signal = x$min_signal,
    smooth = x$smooth,
    id = x$id
  )
}

#' @rdname required_pkgs.step_sec
#' @export
#' @keywords internal
required_pkgs.step_sec_uv_ri_ratio <- function(x, ...) {
  c("measure.sec", "measure")
}
