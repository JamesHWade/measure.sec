# ==============================================================================
# step_sec_dad
#
# Diode array detector processing for multi-wavelength UV
# ==============================================================================

#' Diode Array Detector Processing for SEC
#'
#' `step_sec_dad()` creates a *specification* of a recipe step that processes
#' diode array detector (DAD) signals across multiple wavelengths. It can apply
#' wavelength-specific extinction coefficients and optionally compute ratios to
#' a reference wavelength.
#'
#' @param recipe A recipe object.
#' @param measures Character vector of DAD/UV measure columns. If `NULL`, the
#'   step searches for measure columns containing "dad" or "uv".
#' @param wavelengths Numeric vector of wavelengths (nm) aligned with
#'   `measures`. If `NULL`, attempts to parse wavelengths from `measures` names.
#' @param extinction_coefs Extinction coefficients by wavelength. Accepts a
#'   named numeric vector (names are wavelengths), an unnamed vector aligned
#'   with `wavelengths`, or a data frame with columns `wavelength` and
#'   `extinction_coef`.
#' @param reference_wavelength Optional wavelength used as denominator for
#'   ratio calculations.
#' @param output_prefix Prefix used to name output columns (e.g., `uv_254`).
#' @param path_length Path length of the flow cell in cm. Default is 1.0.
#' @param role Role for generated columns.
#' @param trained Logical indicating if the step has been trained.
#' @param skip Logical. Should the step be skipped when baking?
#' @param id Unique step identifier.
#'
#' @return An updated recipe with the new step added.
#'
#' @details
#' For each wavelength, the signal can be normalized using the Beer-Lambert law:
#'
#' \deqn{A = \varepsilon \times c \times l}
#'
#' where A is absorbance (AU), epsilon is the extinction coefficient, and
#' l is the path length. When `reference_wavelength` is provided, the step
#' additionally creates ratio columns for each wavelength vs the reference.
#'
#' @family sec-detectors
#' @export
#'
#' @examples
#' \dontrun{
#' library(recipes)
#' library(measure)
#'
#' rec <- recipe(~., data = sec_triple_detect) |>
#'   step_measure_input_long(uv_254, location = vars(elution_time), col_name = "uv_254") |>
#'   step_measure_input_long(uv_280, location = vars(elution_time), col_name = "uv_280") |>
#'   step_sec_dad(
#'     measures = c("uv_254", "uv_280"),
#'     wavelengths = c(254, 280),
#'     extinction_coefs = c(`254` = 1.2, `280` = 1.0),
#'     reference_wavelength = 280
#'   ) |>
#'   prep()
#' }
step_sec_dad <- function(
  recipe,
  measures = NULL,
  wavelengths = c(254, 280, 220),
  extinction_coefs = NULL,
  reference_wavelength = NULL,
  output_prefix = "uv",
  path_length = 1.0,
  role = NA,
  trained = FALSE,
  skip = FALSE,
  id = recipes::rand_id("sec_dad")
) {
  if (!is.null(reference_wavelength) && !is.numeric(reference_wavelength)) {
    cli::cli_abort("{.arg reference_wavelength} must be numeric (nm).")
  }

  if (!is.numeric(path_length) || path_length <= 0) {
    cli::cli_abort("{.arg path_length} must be a positive numeric value.")
  }

  recipes::add_step(
    recipe,
    step_sec_dad_new(
      measures = measures,
      wavelengths = wavelengths,
      extinction_coefs = extinction_coefs,
      reference_wavelength = reference_wavelength,
      output_prefix = output_prefix,
      path_length = path_length,
      role = role,
      trained = trained,
      skip = skip,
      id = id
    )
  )
}

step_sec_dad_new <- function(
  measures,
  wavelengths,
  extinction_coefs,
  reference_wavelength,
  output_prefix,
  path_length,
  role,
  trained,
  skip,
  id
) {
  recipes::step(
    subclass = "sec_dad",
    measures = measures,
    wavelengths = wavelengths,
    extinction_coefs = extinction_coefs,
    reference_wavelength = reference_wavelength,
    output_prefix = output_prefix,
    path_length = path_length,
    role = role,
    trained = trained,
    skip = skip,
    id = id
  )
}

#' @export
prep.step_sec_dad <- function(x, training, info = NULL, ...) {
  check_for_measure(training)

  measure_cols <- find_measure_cols(training)

  if (is.null(x$measures)) {
    dad_cols <- measure_cols[grepl("dad|uv", measure_cols, ignore.case = TRUE)]
    if (length(dad_cols) == 0) {
      cli::cli_abort(
        "No DAD/UV measure columns found. Specify {.arg measures} explicitly."
      )
    }
    measures <- dad_cols
  } else {
    measures <- x$measures
  }

  missing_cols <- setdiff(measures, measure_cols)
  if (length(missing_cols) > 0) {
    cli::cli_abort(
      "DAD columns not found in measure columns: {.val {missing_cols}}."
    )
  }

  if (is.null(x$wavelengths)) {
    wavelengths <- .dad_parse_wavelengths(measures)
  } else {
    wavelengths <- x$wavelengths
  }

  if (length(wavelengths) != length(measures)) {
    cli::cli_abort(
      "{.arg wavelengths} must have the same length as {.arg measures}."
    )
  }

  if (!is.numeric(wavelengths) || !all(is.finite(wavelengths))) {
    cli::cli_abort("{.arg wavelengths} must be numeric values (nm).")
  }

  if (
    !is.null(x$reference_wavelength) &&
      !x$reference_wavelength %in% wavelengths
  ) {
    cli::cli_abort(
      "{.arg reference_wavelength} must be one of the provided wavelengths."
    )
  }

  step_sec_dad_new(
    measures = measures,
    wavelengths = wavelengths,
    extinction_coefs = x$extinction_coefs,
    reference_wavelength = x$reference_wavelength,
    output_prefix = x$output_prefix,
    path_length = x$path_length,
    role = x$role,
    trained = TRUE,
    skip = x$skip,
    id = x$id
  )
}

#' @export
bake.step_sec_dad <- function(object, new_data, ...) {
  measures <- object$measures
  wavelengths <- object$wavelengths
  output_prefix <- object$output_prefix
  path_length <- object$path_length

  ext_coefs <- .dad_resolve_extinction_coefs(
    object$extinction_coefs,
    wavelengths
  )

  output_cols <- character(length(measures))
  for (i in seq_along(measures)) {
    in_col <- measures[[i]]
    out_col <- if (is.null(output_prefix) || output_prefix == "") {
      in_col
    } else {
      paste0(output_prefix, "_", wavelengths[[i]])
    }

    ext_value <- ext_coefs[[i]]

    new_data[[out_col]] <- new_measure_list(
      purrr::map(new_data[[in_col]], function(m) {
        values <- m$value
        if (!is.na(ext_value) && ext_value > 0) {
          values <- values / (ext_value * path_length)
        }
        new_measure_tbl(location = m$location, value = values)
      })
    )

    output_cols[[i]] <- out_col
  }

  if (!is.null(object$reference_wavelength)) {
    ref_idx <- which(wavelengths == object$reference_wavelength)
    ref_col <- output_cols[[ref_idx]]

    for (i in seq_along(output_cols)) {
      if (i == ref_idx) {
        next
      }
      ratio_col <- paste0(
        output_cols[[i]],
        "_to_",
        object$reference_wavelength
      )

      new_data[[ratio_col]] <- new_measure_list(
        purrr::map2(
          new_data[[output_cols[[i]]]],
          new_data[[ref_col]],
          function(m, ref_m) {
            ratio_values <- m$value / ref_m$value
            new_measure_tbl(location = m$location, value = ratio_values)
          }
        )
      )
    }
  }

  tibble::as_tibble(new_data)
}

#' @export
print.step_sec_dad <- function(
  x,
  width = max(20, options()$width - 30),
  ...
) {
  title <- "SEC DAD processing"
  if (x$trained) {
    cols <- paste(x$measures, collapse = ", ")
    cat(title, " on ", cols, sep = "")
  } else {
    cat(title)
  }
  cat("\n")
  invisible(x)
}

#' @rdname tidy.step_sec
#' @export
#' @keywords internal
tidy.step_sec_dad <- function(x, ...) {
  tibble::tibble(
    measures = list(x$measures %||% NA_character_),
    wavelengths = list(x$wavelengths %||% NA_real_),
    reference_wavelength = x$reference_wavelength %||% NA_real_,
    output_prefix = x$output_prefix %||% NA_character_,
    id = x$id
  )
}

#' @rdname required_pkgs.step_sec
#' @export
#' @keywords internal
required_pkgs.step_sec_dad <- function(x, ...) {
  c("measure.sec", "measure")
}

.dad_parse_wavelengths <- function(measures) {
  parsed <- vapply(
    measures,
    function(name) {
      match <- regmatches(name, regexpr("[0-9]+", name))
      if (length(match) == 0) {
        NA_real_
      } else {
        as.numeric(match)
      }
    },
    numeric(1)
  )

  if (anyNA(parsed)) {
    cli::cli_abort(
      "Unable to parse wavelengths from {.arg measures}; supply {.arg wavelengths}."
    )
  }

  parsed
}

.dad_resolve_extinction_coefs <- function(extinction_coefs, wavelengths) {
  if (is.null(extinction_coefs)) {
    return(rep(1.0, length(wavelengths)))
  }

  if (is.data.frame(extinction_coefs)) {
    if (!all(c("wavelength", "extinction_coef") %in% names(extinction_coefs))) {
      cli::cli_abort(
        "{.arg extinction_coefs} data frame must include columns `wavelength` and `extinction_coef`."
      )
    }
    values <- extinction_coefs$extinction_coef
    names(values) <- as.character(extinction_coefs$wavelength)
  } else {
    if (!is.numeric(extinction_coefs)) {
      cli::cli_abort("{.arg extinction_coefs} must be numeric.")
    }
    values <- extinction_coefs
  }

  if (length(values) == 1) {
    return(rep(values, length(wavelengths)))
  }

  if (length(values) == length(wavelengths) && is.null(names(values))) {
    return(values)
  }

  if (!is.null(names(values))) {
    out <- rep(NA_real_, length(wavelengths))
    names(out) <- as.character(wavelengths)
    for (i in seq_along(wavelengths)) {
      key <- as.character(wavelengths[[i]])
      if (!key %in% names(values)) {
        cli::cli_abort(
          "Missing extinction coefficient for wavelength {.val {key}}."
        )
      }
      out[[i]] <- values[[key]]
    }
    return(out)
  }

  cli::cli_abort(
    "{.arg extinction_coefs} must be length 1, match {.arg wavelengths}, or be named by wavelength."
  )
}
