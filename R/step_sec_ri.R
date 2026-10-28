# ==============================================================================
# step_sec_ri
#
# RI detector processing with dn/dc handling
# ==============================================================================

#' RI Detector Processing for SEC
#'
#' `step_sec_ri()` creates a *specification* of a recipe step that processes
#' refractive index (RI) detector signals for SEC analysis, including
#' application of dn/dc (refractive index increment) values.
#'
#' @param recipe A recipe object.
#' @param measures Character vector of RI detector column names to process.
#'   If `NULL`, will look for columns containing "ri" in the name.
#' @param dn_dc Refractive index increment (mL/g). Can be:
#'   - A single numeric value applied to all samples
#'   - `NULL` to skip dn/dc normalization (signal remains in detector units)
#' @param concentration_units Units for the signal after dn/dc normalization:
#'   g/mL (default, the natural unit from RI divided by dn/dc) or mg/mL.
#'   This is ignored when dn/dc normalization is skipped.
#' @param dn_dc_column Character name of a column containing sample-specific
#'   positive finite dn/dc values. Overrides `dn_dc` if provided.
#' @param instrument_constant RI detector instrument constant. Default is 1.0
#'   (no adjustment). This converts raw detector response to refractive index
#'   units.
#' @param output_col Name for the output column. Default is to modify in place.
#' @param role Role for generated columns.
#' @param trained Logical indicating if the step has been trained.
#' @param skip Logical. Should the step be skipped when baking?
#' @param id Unique step identifier.
#'
#' @return An updated recipe with the new step added.
#'
#' @details
#' The refractive index (RI) detector is the most common concentration detector
#' in SEC/GPC. The detector response is proportional to the product of
#' concentration and dn/dc:
#'
#' \deqn{RI_{signal} = K \times c \times (dn/dc)}
#'
#' where:
#' - K is the instrument constant
#' - c is mass concentration in g/mL when dn/dc is in mL/g
#' - dn/dc is the refractive index increment (mL/g)
#'
#' Therefore dividing a calibrated refractive-index response by dn/dc gives
#' concentration in g/mL. Set concentration_units to mg/mL to scale those
#' values by 1000. The default preserves the current numeric behavior and can
#' be read directly by scattering steps that use a tagged concentration curve.
#'
#' **Common dn/dc values (in water at 25°C, 633 nm):**
#' - Polystyrene in THF: 0.185 mL/g
#' - PMMA in THF: 0.084 mL/g
#' - PEG in water: 0.135 mL/g
#' - Proteins: ~0.185 mL/g
#' - DNA: ~0.170 mL/g
#'
#' @family sec-detectors
#' @export
#'
#' @examples
#' \dontrun{
#' library(recipes)
#' library(measure)
#'
#' # Apply fixed dn/dc value
#' rec <- recipe(~., data = sec_triple_detect) |>
#'   step_measure_input_long(ri_signal, location = vars(elution_time), col_name = "ri") |>
#'   step_sec_ri(dn_dc = 0.185) |>
#'   prep()
#'
#' # Use sample-specific dn/dc from a column
#' rec <- recipe(~., data = sec_triple_detect) |>
#'   step_measure_input_long(ri_signal, location = vars(elution_time), col_name = "ri") |>
#'   step_sec_ri(dn_dc_column = "dn_dc") |>
#'   prep()
#' }
step_sec_ri <- function(
  recipe,
  measures = NULL,
  dn_dc = NULL,
  dn_dc_column = NULL,
  instrument_constant = 1.0,
  output_col = NULL,
  role = NA,
  trained = FALSE,
  skip = FALSE,
  id = recipes::rand_id("sec_ri"),
  concentration_units = "g/mL"
) {
  .sec_validate_concentration_unit(
    concentration_units,
    "concentration_units",
    allow_null = FALSE
  )

  # Validate inputs
  if (!is.null(dn_dc) && !is.null(dn_dc_column)) {
    cli::cli_warn(
      "{.arg dn_dc_column} takes precedence over {.arg dn_dc}."
    )
  }

  if (!is.null(dn_dc)) {
    .sec_validate_positive_scalar(dn_dc, "dn_dc")
  }

  .sec_validate_positive_scalar(instrument_constant, "instrument_constant")

  recipes::add_step(
    recipe,
    step_sec_ri_new(
      measures = measures,
      dn_dc = dn_dc,
      dn_dc_column = dn_dc_column,
      instrument_constant = instrument_constant,
      concentration_units = concentration_units,
      output_col = output_col,
      role = role,
      trained = trained,
      skip = skip,
      id = id
    )
  )
}

step_sec_ri_new <- function(
  measures,
  dn_dc,
  dn_dc_column,
  instrument_constant,
  concentration_units,
  output_col,
  role,
  trained,
  skip,
  id
) {
  recipes::step(
    subclass = "sec_ri",
    measures = measures,
    dn_dc = dn_dc,
    dn_dc_column = dn_dc_column,
    instrument_constant = instrument_constant,
    concentration_units = concentration_units,
    output_col = output_col,
    role = role,
    trained = trained,
    skip = skip,
    id = id
  )
}

#' @export
prep.step_sec_ri <- function(x, training, info = NULL, ...) {
  check_for_measure(training)

  # Find RI columns if not specified
  if (is.null(x$measures)) {
    measure_cols <- find_measure_cols(training)
    # Look for columns containing "ri" (case insensitive)
    ri_cols <- measure_cols[grepl("ri", measure_cols, ignore.case = TRUE)]
    if (length(ri_cols) == 0) {
      cli::cli_abort(
        "No RI detector columns found. Specify {.arg measures} explicitly."
      )
    }
    measures <- ri_cols
  } else {
    measures <- x$measures
  }

  # Validate dn_dc_column exists if specified
  if (!is.null(x$dn_dc_column)) {
    if (!x$dn_dc_column %in% names(training)) {
      cli::cli_abort(
        "Column {.val {x$dn_dc_column}} not found in training data."
      )
    }
  }

  step_sec_ri_new(
    measures = measures,
    dn_dc = x$dn_dc,
    dn_dc_column = x$dn_dc_column,
    instrument_constant = x$instrument_constant,
    concentration_units = x$concentration_units,
    output_col = x$output_col,
    role = x$role,
    trained = TRUE,
    skip = x$skip,
    id = x$id
  )
}

#' @export
bake.step_sec_ri <- function(object, new_data, ...) {
  measures <- object$measures
  dn_dc <- object$dn_dc
  dn_dc_column <- object$dn_dc_column
  instrument_constant <- object$instrument_constant
  concentration_units <- object$concentration_units
  normalize_to_concentration <- !is.null(dn_dc) || !is.null(dn_dc_column)

  for (col in measures) {
    # Get dn/dc values for each sample
    if (!is.null(dn_dc_column)) {
      dn_dc_values <- new_data[[dn_dc_column]]
    } else if (!is.null(dn_dc)) {
      dn_dc_values <- rep(dn_dc, nrow(new_data))
    } else {
      dn_dc_values <- rep(1.0, nrow(new_data))
    }
    if (
      normalize_to_concentration &&
        (!is.numeric(dn_dc_values) ||
          length(dn_dc_values) != nrow(new_data) ||
          !all(is.finite(dn_dc_values)) ||
          any(dn_dc_values <= 0))
    ) {
      cli::cli_abort("Each dn/dc value must be positive and finite.")
    }

    # Process each sample
    output <- new_measure_list(
      purrr::map2(new_data[[col]], dn_dc_values, function(m, dndc) {
        m$value <- m$value * instrument_constant
        if (normalize_to_concentration) {
          m$value <- m$value / dndc
          if (identical(concentration_units, "mg/mL")) {
            m$value <- m$value * 1000
          }
        }
        m
      })
    )
    if (normalize_to_concentration) {
      attr(output, "sec_concentration_unit") <- concentration_units
    }
    new_data[[col]] <- output
  }

  tibble::as_tibble(new_data)
}

#' @export
print.step_sec_ri <- function(
  x,
  width = max(20, options()$width - 30),
  ...
) {
  title <- "SEC RI detector processing"
  if (!is.null(x$dn_dc)) {
    title <- paste0(title, " (dn/dc = ", x$dn_dc, ")")
  } else if (!is.null(x$dn_dc_column)) {
    title <- paste0(title, " (dn/dc from ", x$dn_dc_column, ")")
  }

  if (x$trained) {
    cols_str <- paste(x$measures, collapse = ", ")
    cat(title, " on ", cols_str, sep = "")
  } else {
    cat(title)
  }
  cat("\n")
  invisible(x)
}

#' @rdname tidy.step_sec
#' @export
#' @keywords internal
tidy.step_sec_ri <- function(x, ...) {
  tibble::tibble(
    measures = list(x$measures),
    dn_dc = x$dn_dc %||% NA_real_,
    dn_dc_column = x$dn_dc_column %||% NA_character_,
    instrument_constant = x$instrument_constant,
    concentration_units = x$concentration_units,
    id = x$id
  )
}

#' @rdname required_pkgs.step_sec
#' @export
#' @keywords internal
required_pkgs.step_sec_ri <- function(x, ...) {
  c("measure.sec", "measure")
}
