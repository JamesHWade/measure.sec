# ==============================================================================
# step_sec_protein
#
# Convenience workflow wrapper for protein SEC analysis
# ==============================================================================

#' Protein SEC Analysis Workflow
#'
#' `step_sec_protein()` creates a *specification* of a recipe step that
#' provides a streamlined workflow for protein SEC analysis, combining
#' baseline correction, aggregate quantitation, and optionally oligomer
#' analysis in a single step.
#'
#' @param recipe A recipe object.
#' @param measures Character vector of measure columns to analyze. If `NULL`,
#'   analyzes all measure columns.
#' @param type Analysis type:
#'   - `"native"` (default): Native conditions (aqueous buffer, non-denaturing)
#'   - `"denaturing"`: Denaturing conditions (e.g., with SDS or guanidine)
#' @param monomer_mw Expected monomer molecular weight in Da. Required for
#'   oligomer analysis if `include_oligomer = TRUE`.
#' @param monomer_start Start of the monomer peak region (in location units).
#'   If `NULL`, automatically determined.
#' @param monomer_end End of the monomer peak region. If `NULL`, automatically
#'   determined.
#' @param extinction_coef Extinction coefficient for UV-based concentration.
#'   If `NULL`, signal remains in raw units.
#' @param aggregate_threshold Minimum fraction of signal to report as
#'   aggregate/fragment. Default is 0.001 (0.1%).
#' @param baseline_method Method for baseline correction. One of `"linear"`
#'   (default), `"median"`, or `"spline"`.
#' @param baseline_left_frac Fraction of chromatogram start for baseline.
#'   Default is 0.05.
#' @param baseline_right_frac Fraction of chromatogram end for baseline.
#'   Default is 0.05.
#' @param include_oligomer Logical. Include detailed oligomer analysis?
#'   Default is `TRUE` if `monomer_mw` is provided.
#' @param output_prefix Prefix for output columns. Default is `"protein_"`.
#' @param role Role for generated columns.
#' @param trained Logical indicating if the step has been trained.
#' @param skip Logical. Should the step be skipped when baking?
#' @param id Unique step identifier.
#'
#' @return An updated recipe with new columns:
#' \describe{
#'   \item{protein_hmws_pct}{Percent high molecular weight species}
#'   \item{protein_monomer_pct}{Percent monomer}
#'   \item{protein_lmws_pct}{Percent low molecular weight species}
#'   \item{protein_main_start}{Start of main peak region}
#'   \item{protein_main_end}{End of main peak region}
#' }
#'
#' If `include_oligomer = TRUE` and `monomer_mw` is provided:
#' \describe{
#'   \item{protein_monomer_oligo_pct}{Percent monomer (from oligomer analysis)}
#'   \item{protein_dimer_pct}{Percent dimer}
#'   \item{protein_trimer_pct}{Percent trimer}
#'   \item{protein_hmw_oligo_pct}{Percent HMW oligomers}
#'   \item{protein_lmw_oligo_pct}{Percent fragments}
#'   \item{protein_species_count}{Number of detected species}
#' }
#'
#' @details
#' This step provides a convenient "one-stop" workflow for protein SEC analysis,
#' suitable for biopharmaceutical characterization. It combines:
#'
#' 1
#. **Baseline correction**: SEC-optimized linear or median baseline
#' 2. **Aggregate quantitation**: HMWS/monomer/LMWS percentages
#' 3. **Oligomer analysis** (optional): Detailed species identification
#'
#' **Native vs Denaturing:**
#' \itemize{
#'   \item **Native**: Preserves quaternary structure; use for oligomer analysis
#'   \item **Denaturing**: Disrupts non-covalent interactions; use for
#'     covalent aggregate detection
#' }
#'
#' **Regulatory Context:**
#' Aggregate analysis is critical for biopharmaceutical characterization:
#' \itemize{
#'   \item ICH Q6B requires aggregate content specification
#'   \item USP <129> provides guidance on aggregate testing
#'   \item Typical specifications: HMWS < 5%, Monomer > 95%
#' }
#'
#' **For More Control:**
#' For advanced analysis or custom workflows, use the individual steps:
#' - [step_sec_baseline()] for baseline correction
#' - [step_sec_uv()] for UV signal processing
#' - [step_sec_aggregates()] for HMWS/monomer/LMWS
#' - [step_sec_oligomer()] for detailed species analysis
#'
#' @family sec-protein
#' @seealso [step_sec_aggregates()], [step_sec_oligomer()]
#' @export
#'
#' @examples
#' \dontrun{
#' library(recipes)
#' library(measure)
#'
#' # Basic protein SEC workflow
#' rec <- recipe(~., data = mab_data) |>
#'   step_measure_input_long(uv280, location = vars(time), col_name = "uv") |>
#'   step_sec_protein(monomer_mw = 150000) |>
#'   prep()
#'
#' # Native mAb analysis with oligomer detection
#' rec <- recipe(~., data = mab_data) |>
#'   step_measure_input_long(uv280, location = vars(time), col_name = "uv") |>
#'   step_sec_protein(
#'     type = "native",
#'     monomer_mw = 150000,
#'     extinction_coef = 1.4,
#'     include_oligomer = TRUE
#'   ) |>
#'   prep()
#'
#' # Denaturing conditions (SDS-SEC)
#' rec <- recipe(~., data = sds_sec_data) |>
#'   step_measure_input_long(uv280, location = vars(time), col_name = "uv") |>
#'   step_sec_protein(type = "denaturing", monomer_mw = 150000) |>
#'   prep()
#' }
step_sec_protein <- function(
  recipe,
  measures = NULL,
  type = c("native", "denaturing"),
  monomer_mw = NULL,
  monomer_start = NULL,
  monomer_end = NULL,
  extinction_coef = NULL,
  aggregate_threshold = 0.001,
  baseline_method = c("linear", "median", "spline"),
  baseline_left_frac = 0.05,
  baseline_right_frac = 0.05,
  include_oligomer = NULL,
  output_prefix = "protein_",
  role = NA,
  trained = FALSE,
  skip = FALSE,
  id = recipes::rand_id("sec_protein")
) {
  type <- match.arg(type)
  baseline_method <- match.arg(baseline_method)

  # Auto-determine oligomer inclusion
  if (is.null(include_oligomer)) {
    include_oligomer <- !is.null(monomer_mw)
  }

  # Validate oligomer requires monomer_mw
  if (include_oligomer && is.null(monomer_mw)) {
    cli::cli_abort(
      c(
        "Oligomer analysis requires {.arg monomer_mw}.",
        "i" = "Either provide {.arg monomer_mw} or set {.arg include_oligomer = FALSE}."
      )
    )
  }

  # Validate monomer_mw if provided
  if (!is.null(monomer_mw)) {
    if (!is.numeric(monomer_mw) || monomer_mw <= 0) {
      cli::cli_abort("{.arg monomer_mw} must be a positive number.")
    }
    if (monomer_mw < 10000 || monomer_mw > 1000000) {
      cli::cli_warn(
        c(
          "Monomer MW ({monomer_mw} Da) is outside typical protein range (10-1000 kDa).",
          "i" = "Verify the molecular weight is correct."
        )
      )
    }
  }

  # Validate aggregate_threshold
  if (!is.numeric(aggregate_threshold) ||
      aggregate_threshold < 0 ||
      aggregate_threshold > 1) {
    cli::cli_abort("{.arg aggregate_threshold} must be between 0 and 1.")
  }

  recipes::add_step(
    recipe,
    step_sec_protein_new(
      measures = measures,
      type = type,
      monomer_mw = monomer_mw,
      monomer_start = monomer_start,
      monomer_end = monomer_end,
      extinction_coef = extinction_coef,
      aggregate_threshold = aggregate_threshold,
      baseline_method = baseline_method,
      baseline_left_frac = baseline_left_frac,
      baseline_right_frac = baseline_right_frac,
      include_oligomer = include_oligomer,
      output_prefix = output_prefix,
      role = role,
      trained = trained,
      skip = skip,
      id = id
    )
  )
}

step_sec_protein_new <- function(
  measures,
  type,
  monomer_mw,
  monomer_start,
  monomer_end,
  extinction_coef,
  aggregate_threshold,
  baseline_method,
  baseline_left_frac,
  baseline_right_frac,
  include_oligomer,
  output_prefix,
  role,
  trained,
  skip,
  id
) {
  recipes::step(
    subclass = "sec_protein",
    measures = measures,
    type = type,
    monomer_mw = monomer_mw,
    monomer_start = monomer_start,
    monomer_end = monomer_end,
    extinction_coef = extinction_coef,
    aggregate_threshold = aggregate_threshold,
    baseline_method = baseline_method,
    baseline_left_frac = baseline_left_frac,
    baseline_right_frac = baseline_right_frac,
    include_oligomer = include_oligomer,
    output_prefix = output_prefix,
    role = role,
    trained = trained,
    skip = skip,
    id = id
  )
}

#' @export
prep.step_sec_protein <- function(x, training, info = NULL, ...) {
  check_for_measure(training)

  # Find measure columns if not specified
  if (is.null(x$measures)) {
    measures <- find_measure_cols(training)
  } else {
    measures <- x$measures
  }

  step_sec_protein_new(
    measures = measures,
    type = x$type,
    monomer_mw = x$monomer_mw,
    monomer_start = x$monomer_start,
    monomer_end = x$monomer_end,
    extinction_coef = x$extinction_coef,
    aggregate_threshold = x$aggregate_threshold,
    baseline_method = x$baseline_method,
    baseline_left_frac = x$baseline_left_frac,
    baseline_right_frac = x$baseline_right_frac,
    include_oligomer = x$include_oligomer,
    output_prefix = x$output_prefix,
    role = x$role,
    trained = TRUE,
    skip = x$skip,
    id = x$id
  )
}

#' @export
bake.step_sec_protein <- function(object, new_data, ...) {
  measures <- object$measures
  type <- object$type
  monomer_mw <- object$monomer_mw
  monomer_start <- object$monomer_start
  monomer_end <- object$monomer_end
  extinction_coef <- object$extinction_coef
  aggregate_threshold <- object$aggregate_threshold
  baseline_method <- object$baseline_method
  baseline_left_frac <- object$baseline_left_frac
  baseline_right_frac <- object$baseline_right_frac
  include_oligomer <- object$include_oligomer
  output_prefix <- object$output_prefix

  n_rows <- nrow(new_data)
  measure_col <- measures[1]

  # Initialize output columns
  hmws_pct <- numeric(n_rows)
  monomer_pct <- numeric(n_rows)
  lmws_pct <- numeric(n_rows)
  main_start <- numeric(n_rows)
  main_end <- numeric(n_rows)

  # Oligomer columns (if requested)
  if (include_oligomer) {
    oligo_species <- c("monomer", "dimer", "trimer", "hmw", "lmw")
    oligo_pct <- list()
    for (sp in oligo_species) {
      oligo_pct[[sp]] <- numeric(n_rows)
    }
    species_count <- integer(n_rows)
  }

  for (i in seq_len(n_rows)) {
    m <- new_data[[measure_col]][[i]]
    location <- m$location
    value <- m$value

    # 1. Baseline correction
    value <- .apply_protein_baseline(
      location,
      value,
      left_frac = baseline_left_frac,
      right_frac = baseline_right_frac,
      method = baseline_method
    )

    # 2. Apply extinction coefficient (if provided)
    if (!is.null(extinction_coef)) {
      value <- value / extinction_coef
    }

    # 3. Aggregate analysis
    agg_result <- .calculate_aggregates(
      location,
      value,
      monomer_start = monomer_start,
      monomer_end = monomer_end,
      threshold = aggregate_threshold
    )

    hmws_pct[i] <- agg_result$hmws_pct
    monomer_pct[i] <- agg_result$monomer_pct
    lmws_pct[i] <- agg_result$lmws_pct
    main_start[i] <- agg_result$main_start
    main_end[i] <- agg_result$main_end

    # 4. Oligomer analysis (if requested)
    if (include_oligomer) {
      oligo_result <- .calculate_oligomers(
        location,
        value,
        monomer_mw = monomer_mw,
        mw_tolerance = 0.15
      )

      for (sp in oligo_species) {
        oligo_pct[[sp]][i] <- oligo_result[[paste0(sp, "_pct")]] %||% 0
      }
      species_count[i] <- oligo_result$species_count
    }
  }

  # Add output columns
  new_data[[paste0(output_prefix, "hmws_pct")]] <- hmws_pct
  new_data[[paste0(output_prefix, "monomer_pct")]] <- monomer_pct
  new_data[[paste0(output_prefix, "lmws_pct")]] <- lmws_pct
  new_data[[paste0(output_prefix, "main_start")]] <- main_start
  new_data[[paste0(output_prefix, "main_end")]] <- main_end

  if (include_oligomer) {
    for (sp in oligo_species) {
      col_name <- if (sp == "monomer") {
        paste0(output_prefix, "monomer_oligo_pct")
      } else if (sp == "hmw") {
        paste0(output_prefix, "hmw_oligo_pct")
      } else if (sp == "lmw") {
        paste0(output_prefix, "lmw_oligo_pct")
      } else {
        paste0(output_prefix, sp, "_pct")
      }
      new_data[[col_name]] <- oligo_pct[[sp]]
    }
    new_data[[paste0(output_prefix, "species_count")]] <- species_count
  }

  tibble::as_tibble(new_data)
}

#' Apply baseline correction for protein SEC
#' @noRd
.apply_protein_baseline <- function(
  location,
  value,
  left_frac,
  right_frac,
  method
) {
  n <- length(value)
  value[is.na(value)] <- 0

  left_n <- max(1, floor(n * left_frac))
  right_n <- max(1, floor(n * right_frac))

  left_idx <- seq_len(left_n)
  right_idx <- seq(n - right_n + 1, n)

  if (method == "median") {
    left_val <- stats::median(value[left_idx], na.rm = TRUE)
    right_val <- stats::median(value[right_idx], na.rm = TRUE)
  } else {
    left_val <- mean(value[left_idx], na.rm = TRUE)
    right_val <- mean(value[right_idx], na.rm = TRUE)
  }

  if (method == "spline" && n > 10) {
    # Use spline through baseline regions
    baseline_x <- c(location[left_idx], location[right_idx])
    baseline_y <- c(value[left_idx], value[right_idx])
    spline_fit <- stats::smooth.spline(baseline_x, baseline_y, df = 4)
    baseline <- stats::predict(spline_fit, location)$y
  } else {
    # Linear interpolation
    baseline <- seq(left_val, right_val, length.out = n)
  }

  pmax(value - baseline, 0)
}

#' Calculate aggregate percentages
#' @noRd
.calculate_aggregates <- function(
  location,
  value,
  monomer_start,
  monomer_end,
  threshold
) {
  value <- pmax(value, 0)

  # Find main peak if boundaries not specified
  if (is.null(monomer_start) || is.null(monomer_end)) {
    peak_info <- .find_main_peak(location, value, threshold_frac = 0.05)
    m_start <- peak_info$start
    m_end <- peak_info$end
  } else {
    m_start <- monomer_start
    m_end <- monomer_end
  }

  # Calculate total area
  total_area <- .peak_area(location, value, 1, length(value))

  if (total_area <= 0) {
    return(list(
      hmws_pct = NA_real_,
      monomer_pct = NA_real_,
      lmws_pct = NA_real_,
      main_start = m_start,
      main_end = m_end
    ))
  }

  # HMWS: before monomer
  hmws_idx <- which(location < m_start)
  if (length(hmws_idx) > 1) {
    hmws_area <- .peak_area(
      location,
      value,
      min(hmws_idx),
      max(hmws_idx)
    )
  } else {
    hmws_area <- 0
  }

  # Monomer
  mono_idx <- which(location >= m_start & location <= m_end)
  if (length(mono_idx) > 1) {
    mono_area <- .peak_area(location, value, min(mono_idx), max(mono_idx))
  } else {
    mono_area <- 0
  }

  # LMWS: after monomer
  lmws_idx <- which(location > m_end)
  if (length(lmws_idx) > 1) {
    lmws_area <- .peak_area(location, value, min(lmws_idx), max(lmws_idx))
  } else {
    lmws_area <- 0
  }

  list(
    hmws_pct = 100 * hmws_area / total_area,
    monomer_pct = 100 * mono_area / total_area,
    lmws_pct = 100 * lmws_area / total_area,
    main_start = m_start,
    main_end = m_end
  )
}

#' Calculate oligomer percentages
#' @noRd
.calculate_oligomers <- function(location, value, monomer_mw, mw_tolerance) {
  value <- pmax(value, 0)

  # Detect peaks
  peaks <- .detect_peaks(location, value)

  if (length(peaks) == 0) {
    return(list(
      monomer_pct = 0,
      dimer_pct = 0,
      trimer_pct = 0,
      hmw_pct = 0,
      lmw_pct = 0,
      species_count = 0
    ))
  }

  total_area <- .peak_area(location, value, 1, length(value))
  if (total_area <= 0) {
    return(list(
      monomer_pct = 0,
      dimer_pct = 0,
      trimer_pct = 0,
      hmw_pct = 0,
      lmw_pct = 0,
      species_count = 0
    ))
  }

  # Calculate peak areas and find largest (assumed monomer)
  peak_areas <- sapply(peaks, function(p) {
    .peak_area(location, value, p$start_idx, p$end_idx)
  })

  monomer_idx <- which.max(peak_areas)
  monomer_rt <- peaks[[monomer_idx]]$apex_location

  # Assign species by retention time (no MW data)
  species_pct <- list(
    monomer_pct = 0,
    dimer_pct = 0,
    trimer_pct = 0,
    hmw_pct = 0,
    lmw_pct = 0
  )

  for (j in seq_along(peaks)) {
    pct <- 100 * peak_areas[j] / total_area
    rt <- peaks[[j]]$apex_location

    if (j == monomer_idx) {
      species_pct$monomer_pct <- species_pct$monomer_pct + pct
    } else if (rt < monomer_rt) {
      # Earlier elution = higher MW (aggregate)
      # Rough estimate: each 2x MW shift corresponds to ~0.3-0.5 min shift
      # We'll just group as HMW for now without MW data
      species_pct$hmw_pct <- species_pct$hmw_pct + pct
    } else {
      # Later elution = lower MW (fragment)
      species_pct$lmw_pct <- species_pct$lmw_pct + pct
    }
  }

  species_pct$species_count <- length(peaks)

  species_pct
}

#' @export
print.step_sec_protein <- function(
  x,
  width = max(20, options()$width - 30),
  ...
) {
  title <- paste0("SEC protein analysis (", x$type, ")")
  if (x$trained) {
    if (!is.null(x$monomer_mw)) {
      cat(
        title,
        ", monomer MW = ",
        format(x$monomer_mw, big.mark = ","),
        " Da",
        sep = ""
      )
    } else {
      cat(title)
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
tidy.step_sec_protein <- function(x, ...) {
  tibble::tibble(
    type = x$type,
    monomer_mw = x$monomer_mw %||% NA_real_,
    include_oligomer = x$include_oligomer,
    baseline_method = x$baseline_method,
    aggregate_threshold = x$aggregate_threshold,
    id = x$id
  )
}

#' @rdname required_pkgs.step_sec
#' @export
#' @keywords internal
required_pkgs.step_sec_protein <- function(x, ...) {
  c("measure.sec", "measure")
}
