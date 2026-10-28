# ==============================================================================
# step_sec_uv
#
# UV detector processing with extinction coefficient handling
# ==============================================================================

#' UV Detector Processing for SEC
#'
#' `step_sec_uv()` creates a *specification* of a recipe step that processes
#' UV detector signals for SEC analysis, including application of extinction
#' coefficients for concentration determination.
#'
#' @param recipe A recipe object.
#' @param measures Character vector of UV detector column names to process.
#'   If `NULL`, will look for columns containing "uv" in the name.
#' @param extinction_coef Extinction coefficient in mL/(mg*cm) or L/(g*cm).
#'   Can be:
#'   - A single numeric value applied to all samples
#'   - `NULL` to skip normalization (signal remains in AU)
#' @param concentration_units Units for the signal after extinction
#'   normalization: mg/mL (default) or g/mL. This is ignored when extinction
#'   normalization is skipped.
#' @param extinction_column Character name of a column containing sample-specific
#'   extinction coefficients. Overrides `extinction_coef` if provided. Values
#'   must be positive and finite because a zero coefficient cannot be inverted
#'   to determine concentration from absorbance.
#' @param wavelength UV detection wavelength in nm. For documentation only.
#' @param path_length Path length of the flow cell in cm. Default is 1.0.
#' @param output_col Name for the output column. Default is to modify in place.
#' @param role Role for generated columns.
#' @param trained Logical indicating if the step has been trained.
#' @param skip Logical. Should the step be skipped when baking?
#' @param id Unique step identifier.
#'
#' @return An updated recipe with the new step added.
#'
#' @details
#' The UV detector measures absorbance according to the Beer-Lambert law:
#'
#' \deqn{A = \varepsilon \times c \times l}
#'
#' where:
#' - A is absorbance (AU)
#' - epsilon is a mass-specific extinction coefficient in mL/(mg*cm)
#' - c is the concentration in mg/mL
#' - l is the path length (cm)
#'
#' This step can convert UV absorbance to concentration-proportional signals
#' by dividing by the extinction coefficient and path length. With the
#' documented coefficient units, the result is mg/mL; set concentration_units
#' to g/mL to divide by 1000.
#'
#' **Common UV applications in SEC:**
#' - Proteins at 280 nm (aromatic amino acids)
#' - Nucleic acids at 260 nm
#' - Conjugated polymers
#' - UV-active end groups or labels
#'
#' **UV vs RI for concentration:**
#' - UV is more sensitive for chromophore-containing analytes
#' - UV response depends on chemical composition (may vary with MW)
#' - RI is more universal but less sensitive
#' - For accurate MW, combine both detectors
#'
#' @family sec-detectors
#' @export
#'
#' @examples
#' \dontrun{
#' library(recipes)
#' library(measure)
#'
#' # Apply fixed extinction coefficient
#' rec <- recipe(~., data = sec_triple_detect) |>
#'   step_measure_input_long(uv_signal, location = vars(elution_time), col_name = "uv") |>
#'   step_sec_uv(extinction_coef = 1.0, wavelength = 280) |>
#'   prep()
#'
#' # Use sample-specific extinction coefficients
#' rec <- recipe(~., data = sec_triple_detect) |>
#'   step_measure_input_long(uv_signal, location = vars(elution_time), col_name = "uv") |>
#'   step_sec_uv(extinction_column = "ext_coef") |>
#'   prep()
#' }
step_sec_uv <- function(
  recipe,
  measures = NULL,
  extinction_coef = NULL,
  extinction_column = NULL,
  wavelength = 280,
  path_length = 1.0,
  output_col = NULL,
  role = NA,
  trained = FALSE,
  skip = FALSE,
  id = recipes::rand_id("sec_uv"),
  concentration_units = "mg/mL"
) {
  .sec_validate_concentration_unit(
    concentration_units,
    "concentration_units",
    allow_null = FALSE
  )

  # Validate inputs
  if (!is.null(extinction_coef) && !is.null(extinction_column)) {
    cli::cli_warn(
      "{.arg extinction_column} takes precedence over {.arg extinction_coef}."
    )
  }

  if (!is.null(extinction_coef)) {
    .sec_validate_positive_scalar(extinction_coef, "extinction_coef")
  }

  .sec_validate_positive_scalar(wavelength, "wavelength")
  .sec_validate_positive_scalar(path_length, "path_length")

  recipes::add_step(
    recipe,
    step_sec_uv_new(
      measures = measures,
      extinction_coef = extinction_coef,
      extinction_column = extinction_column,
      concentration_units = concentration_units,
      wavelength = wavelength,
      path_length = path_length,
      output_col = output_col,
      role = role,
      trained = trained,
      skip = skip,
      id = id
    )
  )
}

step_sec_uv_new <- function(
  measures,
  extinction_coef,
  extinction_column,
  concentration_units,
  wavelength,
  path_length,
  output_col,
  role,
  trained,
  skip,
  id
) {
  recipes::step(
    subclass = "sec_uv",
    measures = measures,
    extinction_coef = extinction_coef,
    extinction_column = extinction_column,
    concentration_units = concentration_units,
    wavelength = wavelength,
    path_length = path_length,
    output_col = output_col,
    role = role,
    trained = trained,
    skip = skip,
    id = id
  )
}

#' @export
prep.step_sec_uv <- function(x, training, info = NULL, ...) {
  check_for_measure(training)

  # Find UV columns if not specified
  if (is.null(x$measures)) {
    measure_cols <- find_measure_cols(training)
    # Look for columns containing "uv" (case insensitive)
    uv_cols <- measure_cols[grepl("uv", measure_cols, ignore.case = TRUE)]
    if (length(uv_cols) == 0) {
      cli::cli_abort(
        "No UV detector columns found. Specify {.arg measures} explicitly."
      )
    }
    measures <- uv_cols
  } else {
    measures <- x$measures
  }

  # Validate extinction_column exists if specified
  if (!is.null(x$extinction_column)) {
    if (!x$extinction_column %in% names(training)) {
      cli::cli_abort(
        "Column {.val {x$extinction_column}} not found in training data."
      )
    }
  }

  step_sec_uv_new(
    measures = measures,
    extinction_coef = x$extinction_coef,
    extinction_column = x$extinction_column,
    concentration_units = x$concentration_units,
    wavelength = x$wavelength,
    path_length = x$path_length,
    output_col = x$output_col,
    role = x$role,
    trained = TRUE,
    skip = x$skip,
    id = x$id
  )
}

#' @export
bake.step_sec_uv <- function(object, new_data, ...) {
  measures <- object$measures
  extinction_coef <- object$extinction_coef
  extinction_column <- object$extinction_column
  concentration_units <- object$concentration_units
  normalize_to_concentration <- !is.null(extinction_coef) ||
    !is.null(extinction_column)
  path_length <- object$path_length

  for (col in measures) {
    # Get extinction coefficient values for each sample
    if (!is.null(extinction_column)) {
      ext_values <- new_data[[extinction_column]]
    } else if (!is.null(extinction_coef)) {
      ext_values <- rep(extinction_coef, nrow(new_data))
    } else {
      ext_values <- rep(1.0, nrow(new_data))
    }
    if (
      normalize_to_concentration &&
        (!is.numeric(ext_values) ||
          length(ext_values) != nrow(new_data) ||
          !all(is.finite(ext_values)) ||
          any(ext_values <= 0))
    ) {
      cli::cli_abort(
        "Each extinction coefficient must be positive and finite."
      )
    }

    # Process each sample
    output <- new_measure_list(
      purrr::map2(new_data[[col]], ext_values, function(m, ext) {
        # Convert absorbance to concentration-proportional signal
        # A = epsilon * c * l, so c ~ A / (epsilon * l)
        if (normalize_to_concentration) {
          m$value <- m$value / (ext * path_length)
          if (identical(concentration_units, "g/mL")) {
            m$value <- m$value / 1000
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
print.step_sec_uv <- function(
  x,
  width = max(20, options()$width - 30),
  ...
) {
  title <- paste0("SEC UV detector processing (", x$wavelength, " nm)")
  if (!is.null(x$extinction_coef)) {
    title <- paste0(title, ", e = ", x$extinction_coef)
  } else if (!is.null(x$extinction_column)) {
    title <- paste0(title, ", e from ", x$extinction_column)
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
tidy.step_sec_uv <- function(x, ...) {
  tibble::tibble(
    measures = list(x$measures),
    extinction_coef = x$extinction_coef %||% NA_real_,
    extinction_column = x$extinction_column %||% NA_character_,
    concentration_units = x$concentration_units,
    wavelength = x$wavelength,
    path_length = x$path_length,
    id = x$id
  )
}

#' @rdname required_pkgs.step_sec
#' @export
#' @keywords internal
required_pkgs.step_sec_uv <- function(x, ...) {
  c("measure.sec", "measure")
}
