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
#' @param measures Character vector of RI detector columns to process. If
#'   `NULL`, the step looks for columns containing "ri" in the name. Multiple
#'   columns can be processed in place; a named output requires exactly one.
#' @param dn_dc Refractive index increment (mL/g). Can be:
#'   - A single finite, nonzero numeric value applied to all samples
#'   - `NULL` to skip dn/dc normalization (signal remains in detector units)
#' @param concentration_units Units for the signal after dn/dc normalization:
#'   g/mL (default, the natural unit from RI divided by dn/dc) or mg/mL.
#'   This is ignored when dn/dc normalization is skipped.
#' @param dn_dc_column Character name of a column containing sample-specific
#'   finite, nonzero dn/dc values. Values may be positive or negative.
#'   Overrides `dn_dc` if provided.
#' @param instrument_constant RI detector instrument constant. Default is 1.0
#'   (no adjustment). This positive conversion gain multiplies raw detector
#'   response to obtain refractive-index difference.
#' @param output_col Name for the output column. `NULL` modifies the selected
#'   RI measure in place. A different name creates a converted column and keeps
#'   the source measure unchanged.
#' @param role Role for generated columns.
#' @param trained Logical indicating if the step has been trained.
#' @param skip Logical. Should the step be skipped when baking?
#' @param id Unique step identifier.
#'
#' @return An updated recipe with the new step added.
#'
#' @details
#' The RI detector response is proportional to concentration and dn/dc. Let
#' `S_raw` be the stored detector response and `K` be the positive conversion
#' gain supplied as `instrument_constant`. The calibrated refractive-index
#' difference and concentration are related by:
#'
#' \deqn{\Delta n = K \times S_{raw} = c \times (dn/dc)}
#'
#' Thus the step calculates concentration in g/mL as
#' `S_raw * instrument_constant / dn_dc`, when dn/dc is in mL/g. A negative
#' dn/dc is physically valid for some solute/solvent combinations; the
#' [Brookhaven BI-DNDC manual](https://www.brookhaveninstruments.com/wp-content/uploads/2023/08/BI-DNDC-Manual-Static-Mode.pdf)
#' describes negative values. Wyatt notes that dn/dc depends on solute,
#' solvent, wavelength, and temperature, so supply a value measured or
#' qualified for the relevant conditions
#' ([application note](https://www.wyatt.com/library/application-notes/determination-of-dn-dc.html)).
#'
#' Set `concentration_units` to `"mg/mL"` to scale the result by 1000. The
#' default preserves the existing numeric behavior and can be read directly by
#' scattering steps that use a tagged concentration curve. With `dn_dc = NULL`
#' and no `dn_dc_column`, no concentration normalization is performed. Using
#' the default `instrument_constant = 1.0` leaves signal values unchanged and
#' adds no concentration-unit tag.
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
  .sec_validate_ri_output_col(output_col)

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
    .sec_validate_nonzero_scalar(dn_dc, "dn_dc")
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

.sec_validate_ri_output_col <- function(output_col) {
  if (is.null(output_col)) {
    return(invisible(NULL))
  }

  if (
    !is.character(output_col) ||
      length(output_col) != 1L ||
      is.na(output_col) ||
      !nzchar(trimws(output_col))
  ) {
    cli::cli_abort(
      "{.arg output_col} must be NULL or a single non-missing, non-blank character name."
    )
  }

  invisible(output_col)
}

.sec_validate_ri_measures <- function(measures, require_single = FALSE) {
  if (
    !is.character(measures) ||
      anyNA(measures) ||
      !all(nzchar(measures))
  ) {
    cli::cli_abort("RI measures must be non-missing, non-blank column names.")
  }

  if (require_single && length(measures) != 1L) {
    cli::cli_abort("Select exactly one RI measure column for a named output.")
  }
  if (length(measures) == 0L) {
    cli::cli_abort("Select at least one RI measure column.")
  }

  invisible(measures)
}

.sec_check_ri_output_collision <- function(
  output_col,
  source_col,
  data,
  stage
) {
  if (
    !is.null(output_col) &&
      !identical(unname(output_col), unname(source_col)) &&
      output_col %in% names(data)
  ) {
    cli::cli_abort(
      "Output column {.val {output_col}} already exists in {stage} data."
    )
  }
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

  measure_cols <- find_measure_cols(training)

  # Find RI columns if not specified
  if (is.null(x$measures)) {
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

  .sec_validate_ri_measures(
    measures,
    require_single = !is.null(x$output_col)
  )
  .sec_validate_ri_output_col(x$output_col)

  missing_measures <- setdiff(measures, names(training))
  if (length(missing_measures) > 0L) {
    cli::cli_abort(
      "RI measure column{?s} {.and {.val {missing_measures}}} not found in training data."
    )
  }
  invalid_measures <- setdiff(measures, measure_cols)
  if (length(invalid_measures) > 0L) {
    cli::cli_abort(
      "Selected column{?s} {.and {.val {invalid_measures}}} {?is/are} not measure column{?s}."
    )
  }
  if (!is.null(x$output_col)) {
    .sec_check_ri_output_collision(
      x$output_col,
      measures[[1]],
      training,
      "training"
    )
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
  .sec_validate_ri_measures(
    measures,
    require_single = !is.null(object$output_col)
  )
  .sec_validate_ri_output_col(object$output_col)

  missing_measures <- setdiff(measures, names(new_data))
  if (length(missing_measures) > 0L) {
    cli::cli_abort(
      "RI measure column{?s} {.and {.val {missing_measures}}} not found in new data."
    )
  }
  if (!is.null(object$output_col)) {
    .sec_check_ri_output_collision(
      object$output_col,
      measures[[1]],
      new_data,
      "baked"
    )
  }

  dn_dc <- object$dn_dc
  dn_dc_column <- object$dn_dc_column
  instrument_constant <- object$instrument_constant
  concentration_units <- object$concentration_units
  normalize_to_concentration <- !is.null(dn_dc) || !is.null(dn_dc_column)

  for (source_col in measures) {
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
          any(dn_dc_values == 0))
    ) {
      cli::cli_abort("Each dn/dc value must be nonzero and finite.")
    }

    # Process each sample
    output <- new_measure_list(
      purrr::map2(new_data[[source_col]], dn_dc_values, function(m, dndc) {
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
    output_col <- object$output_col %||% source_col
    new_data[[output_col]] <- output
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
    output_col = x$output_col %||% NA_character_,
    id = x$id
  )
}

#' @rdname required_pkgs.step_sec
#' @export
#' @keywords internal
required_pkgs.step_sec_ri <- function(x, ...) {
  c("measure.sec", "measure")
}
