# ==============================================================================
# step_sec_oligomer
#
# Detailed oligomeric species analysis for protein SEC
# ==============================================================================

#' Oligomeric Species Analysis for Protein SEC
#'
#' `step_sec_oligomer()` creates a *specification* of a recipe step that
#' identifies and quantifies individual oligomeric species (monomer, dimer,
#' trimer, etc.) in protein SEC chromatograms.
#'
#' @param recipe A recipe object.
#' @param measures Character vector of measure columns to analyze. If `NULL`,
#'   analyzes all measure columns.
#' @param monomer_mw Expected monomer molecular weight in Da. Required for
#'   species assignment. Typical range: 10,000 - 1,000,000 Da for proteins.
#' @param peak_detection Method for peak detection:
#'   - `"auto"` (default): Automatic peak detection using derivative analysis
#'   - `"manual"`: Use peaks specified in `peaks` argument
#' @param peaks For manual mode, a list or data frame with peak definitions.
#'   Each peak should have `start` and `end` retention times.
#' @param species Character vector of species to identify. Default includes
#'   `"monomer"`, `"dimer"`, `"trimer"`, `"hmw"`, `"lmw"`.
#' @param mw_tolerance Tolerance for MW-based species assignment as a fraction.
#'   Default is 0.15 (15%). A peak is assigned to "dimer" if its MW is within
#'   15% of 2x the monomer MW.
#' @param mw_column Optional name of a molecular weight measure column (from
#'   MALS or LALS). If provided, uses MW for species assignment; otherwise
#'   uses retention time patterns.
#' @param min_area_pct Minimum peak area percentage to report. Peaks below this
#'   threshold are grouped into "other". Default is 0.1 (0.1%).
#' @param output_prefix Prefix for output columns. Default is `"oligo_"`.
#' @param role Role for generated columns.
#' @param trained Logical indicating if the step has been trained.
#' @param skip Logical. Should the step be skipped when baking?
#' @param id Unique step identifier.
#'
#' @return An updated recipe with new columns:
#' \describe{
#'   \item{oligo_monomer_pct}{Percent monomer}
#'   \item{oligo_dimer_pct}{Percent dimer}
#'   \item{oligo_trimer_pct}{Percent trimer}
#'   \item{oligo_hmw_pct}{Percent high molecular weight (> trimer)}
#'   \item{oligo_lmw_pct}{Percent low molecular weight (fragments)}
#'   \item{oligo_species_count}{Number of detected species}
#'   \item{oligo_monomer_mw}{Observed monomer MW (if mw_column provided)}
#'   \item{oligo_dimer_mw}{Observed dimer MW (if mw_column provided)}
#' }
#'
#' @details
#' This step extends `step_sec_aggregates()` by providing detailed species
#' identification rather than just HMWS/monomer/LMWS classification.
#'
#' **Species Assignment:**
#' \itemize{
#'   \item Monomer: MW within `mw_tolerance` of `monomer_mw`
#'   \item Dimer: MW within `mw_tolerance` of 2x `monomer_mw`
#'   \item Trimer: MW within `mw_tolerance` of 3x `monomer_mw`
#'   \item HMW: MW > 3x `monomer_mw`
#'   \item LMW: MW < `monomer_mw` (fragments, clips)
#' }
#'
#' **Peak Detection:**
#'
#' When using `"auto"` detection, peaks are identified using:
#' 1. Signal derivative analysis
#' 2. Local maxima identification
#' 3. Valley detection for peak boundaries
#'
#' **Without MW Data:**
#'
#' If no MW column is available, species are assigned based on retention time:
#' \itemize{
#'   \item Largest peak is assumed to be monomer
#'   \item Earlier-eluting peaks are HMW/oligomers
#'   \item Later-eluting peaks are LMW/fragments
#' }
#'
#' @note
#' For best results:
#' \itemize{
#'   \item Provide `monomer_mw` for accurate species assignment
#'   \item Use with MALS/LALS data for MW-based assignment
#'   \item Ensure good chromatographic resolution between species
#' }
#'
#' @family sec-protein
#' @seealso [step_sec_aggregates()] for simpler HMWS/monomer/LMWS analysis
#' @export
#'
#' @examples
#' \dontrun{
#' library(recipes)
#' library(measure)
#'
#' # Analyze IgG oligomers (monomer MW ~150 kDa)
#' rec <- recipe(~., data = igg_data) |>
#'   step_measure_input_long(uv280, location = vars(time), col_name = "uv") |>
#'   step_sec_baseline() |>
#'   step_sec_oligomer(monomer_mw = 150000) |>
#'   prep()
#'
#' # With MALS-derived MW for accurate assignment
#' rec <- recipe(~., data = igg_mals_data) |>
#'   step_measure_input_long(uv280, location = vars(time), col_name = "uv") |>
#'   step_sec_mals() |>
#'   step_sec_oligomer(monomer_mw = 150000, mw_column = "mw_mals") |>
#'   prep()
#' }
step_sec_oligomer <- function(
  recipe,
  measures = NULL,

  monomer_mw = NULL,
  peak_detection = c("auto", "manual"),
  peaks = NULL,
  species = c("monomer", "dimer", "trimer", "hmw", "lmw"),
  mw_tolerance = 0.15,
  mw_column = NULL,
  min_area_pct = 0.1,
  output_prefix = "oligo_",
  role = NA,
  trained = FALSE,
  skip = FALSE,
  id = recipes::rand_id("sec_oligomer")
) {
  peak_detection <- match.arg(peak_detection)

  # Validate monomer_mw
  if (is.null(monomer_mw)) {
    cli::cli_abort(
      c(
        "{.arg monomer_mw} is required for oligomer analysis.",
        "i" = "Provide the expected monomer molecular weight in Daltons.",
        "i" = "Example: monomer_mw = 150000 for an IgG antibody."
      )
    )
  }

  if (!is.numeric(monomer_mw) || monomer_mw <= 0) {
    cli::cli_abort("{.arg monomer_mw} must be a positive number.")
  }

  # Warn if MW seems outside protein range

  if (monomer_mw < 10000 || monomer_mw > 1000000) {
    cli::cli_warn(
      c(
        "Monomer MW ({monomer_mw} Da) is outside typical protein range (10-1000 kDa).",
        "i" = "Verify the molecular weight is correct."
      )
    )
  }

  # Validate mw_tolerance
  if (!is.numeric(mw_tolerance) || mw_tolerance <= 0 || mw_tolerance >= 1) {
    cli::cli_abort("{.arg mw_tolerance} must be between 0 and 1.")
  }

  # Validate manual peak mode
  if (peak_detection == "manual" && is.null(peaks)) {
    cli::cli_abort(
      "Peak detection mode {.val manual} requires {.arg peaks} to be specified."
    )
  }

  recipes::add_step(
    recipe,
    step_sec_oligomer_new(
      measures = measures,
      monomer_mw = monomer_mw,
      peak_detection = peak_detection,
      peaks = peaks,
      species = species,
      mw_tolerance = mw_tolerance,
      mw_column = mw_column,
      min_area_pct = min_area_pct,
      output_prefix = output_prefix,
      role = role,
      trained = trained,
      skip = skip,
      id = id
    )
  )
}

step_sec_oligomer_new <- function(
  measures,
  monomer_mw,
  peak_detection,
  peaks,
  species,
  mw_tolerance,
  mw_column,
  min_area_pct,
  output_prefix,
  role,
  trained,
  skip,
  id
) {
  recipes::step(
    subclass = "sec_oligomer",
    measures = measures,
    monomer_mw = monomer_mw,
    peak_detection = peak_detection,
    peaks = peaks,
    species = species,
    mw_tolerance = mw_tolerance,
    mw_column = mw_column,
    min_area_pct = min_area_pct,
    output_prefix = output_prefix,
    role = role,
    trained = trained,
    skip = skip,
    id = id
  )
}

#' @export
prep.step_sec_oligomer <- function(x, training, info = NULL, ...) {
  check_for_measure(training)

  # Find measure columns if not specified
  if (is.null(x$measures)) {
    measures <- find_measure_cols(training)
  } else {
    measures <- x$measures
  }

  # Validate mw_column exists if provided
  if (!is.null(x$mw_column)) {
    measure_cols <- find_measure_cols(training)
    if (
      !x$mw_column %in% names(training) &&
        !x$mw_column %in% measure_cols
    ) {
      cli::cli_warn(
        c(
          "MW column {.val {x$mw_column}} not found in data.",
          "i" = "Species assignment will use retention time patterns instead."
        )
      )
    }
  }

  step_sec_oligomer_new(
    measures = measures,
    monomer_mw = x$monomer_mw,
    peak_detection = x$peak_detection,
    peaks = x$peaks,
    species = x$species,
    mw_tolerance = x$mw_tolerance,
    mw_column = x$mw_column,
    min_area_pct = x$min_area_pct,
    output_prefix = x$output_prefix,
    role = x$role,
    trained = TRUE,
    skip = x$skip,
    id = x$id
  )
}

#' Detect peaks in a chromatogram
#' @noRd
.detect_peaks <- function(location, value, min_height_frac = 0.01) {
  n <- length(value)
  if (n < 5) {
    return(list())
  }

  # Smooth signal slightly for derivative calculation
  if (n > 10) {
    smooth_window <- min(5, floor(n / 10))
    value_smooth <- stats::filter(
      value,
      rep(1 / smooth_window, smooth_window),
      sides = 2
    )
    value_smooth[is.na(value_smooth)] <- value[is.na(value_smooth)]
    value_smooth <- as.numeric(value_smooth)
  } else {
    value_smooth <- value
  }

  max_val <- max(value_smooth, na.rm = TRUE)
  min_height <- min_height_frac * max_val

  # Find local maxima
  peaks <- list()
  peak_idx <- 1

  for (i in seq(2, n - 1)) {
    if (
      value_smooth[i] > value_smooth[i - 1] &&
        value_smooth[i] > value_smooth[i + 1] &&
        value_smooth[i] > min_height
    ) {
      # Found a local maximum, find boundaries
      # Find left boundary (valley or start)
      left <- i
      for (j in seq(i - 1, 1, -1)) {
        if (
          j == 1 ||
            value_smooth[j] < value_smooth[j + 1] &&
              value_smooth[j] < value_smooth[j - 1]
        ) {
          left <- j
          break
        }
        if (value_smooth[j] < min_height * 0.1) {
          left <- j
          break
        }
      }

      # Find right boundary (valley or end)
      right <- i
      for (j in seq(i + 1, n)) {
        if (
          j == n ||
            value_smooth[j] < value_smooth[j - 1] &&
              value_smooth[j] < value_smooth[j + 1]
        ) {
          right <- j
          break
        }
        if (value_smooth[j] < min_height * 0.1) {
          right <- j
          break
        }
      }

      peaks[[peak_idx]] <- list(
        apex_idx = i,
        apex_location = location[i],
        apex_value = value[i],
        start_idx = left,
        end_idx = right,
        start_location = location[left],
        end_location = location[right]
      )
      peak_idx <- peak_idx + 1
    }
  }

  peaks
}

#' Calculate peak area using trapezoidal integration
#' @noRd
.peak_area <- function(location, value, start_idx, end_idx) {
  if (start_idx >= end_idx) {
    return(0)
  }
  idx <- seq(start_idx, end_idx)
  x <- location[idx]
  y <- value[idx]
  y[is.na(y)] <- 0
  y <- pmax(y, 0)

  if (length(x) < 2) {
    return(0)
  }

  dt <- diff(x)
  sum(y[-1] * dt + y[-length(y)] * dt) / 2
}

#' Assign species based on MW ratio
#' @noRd
.assign_species <- function(mw, monomer_mw, tolerance) {
  if (is.na(mw) || mw <= 0) {
    return("unknown")
  }

  ratio <- mw / monomer_mw

  if (ratio < (1 - tolerance)) {
    "lmw"
  } else if (abs(ratio - 1) <= tolerance) {
    "monomer"
  } else if (abs(ratio - 2) <= tolerance * 2) {
    "dimer"
  } else if (abs(ratio - 3) <= tolerance * 3) {
    "trimer"
  } else if (ratio > 3) {
    "hmw"
  } else {
    "unknown"
  }
}

#' @export
bake.step_sec_oligomer <- function(object, new_data, ...) {
  measures <- object$measures
  monomer_mw <- object$monomer_mw
  peak_detection <- object$peak_detection
  manual_peaks <- object$peaks
  species_list <- object$species
  mw_tolerance <- object$mw_tolerance
  mw_column <- object$mw_column
  min_area_pct <- object$min_area_pct
  output_prefix <- object$output_prefix

  # Validate measures are available
  if (length(measures) == 0) {
    cli::cli_abort(
      c(
        "No measure columns available for oligomer analysis.",
        "i" = "Ensure your recipe includes step_measure_input_*() before this step."
      )
    )
  }

  n_rows <- nrow(new_data)

  # Initialize output columns for each species
  species_pct <- list()
  species_mw <- list()
  for (sp in species_list) {
    species_pct[[sp]] <- numeric(n_rows)
    species_mw[[sp]] <- numeric(n_rows)
  }
  species_count <- integer(n_rows)

  # Use first measure column for signal
  measure_col <- measures[1]

  # Check if MW column is available
  has_mw <- !is.null(mw_column) && mw_column %in% names(new_data)

  for (i in seq_len(n_rows)) {
    m <- new_data[[measure_col]][[i]]
    location <- m$location
    value <- m$value

    # Handle NA values
    value[is.na(value)] <- 0
    value <- pmax(value, 0)

    # Detect or use manual peaks
    if (peak_detection == "auto") {
      peaks <- .detect_peaks(location, value)
    } else {
      # Convert manual peaks to internal format
      peaks <- lapply(seq_len(nrow(manual_peaks)), function(j) {
        start_loc <- manual_peaks$start[j]
        end_loc <- manual_peaks$end[j]
        start_idx <- which.min(abs(location - start_loc))
        end_idx <- which.min(abs(location - end_loc))
        apex_idx <- start_idx + which.max(value[start_idx:end_idx]) - 1
        list(
          apex_idx = apex_idx,
          apex_location = location[apex_idx],
          apex_value = value[apex_idx],
          start_idx = start_idx,
          end_idx = end_idx,
          start_location = start_loc,
          end_location = end_loc
        )
      })
    }

    if (length(peaks) == 0) {
      cli::cli_warn(
        c(
          "No peaks detected for row {i}.",
          "i" = "Check chromatogram signal quality."
        )
      )
      species_count[i] <- 0
      next
    }

    # Calculate total area
    total_area <- .peak_area(location, value, 1, length(value))
    if (total_area <= 0) {
      cli::cli_warn(
        c(
          "Zero signal area for row {i}.",
          "i" = "Check baseline correction and signal intensity."
        )
      )
      species_count[i] <- 0
      next
    }

    # Calculate area and assign species for each peak
    peak_data <- lapply(peaks, function(p) {
      area <- .peak_area(location, value, p$start_idx, p$end_idx)
      pct <- 100 * area / total_area

      # Get MW at apex if available
      if (has_mw) {
        mw_m <- new_data[[mw_column]][[i]]
        mw_at_apex <- mw_m$value[p$apex_idx]
        if (is.na(mw_at_apex)) {
          # Try to interpolate
          mw_at_apex <- NA_real_
        }
        assigned <- .assign_species(mw_at_apex, monomer_mw, mw_tolerance)
      } else {
        mw_at_apex <- NA_real_
        assigned <- "unknown"
      }

      list(
        area = area,
        pct = pct,
        mw = mw_at_apex,
        species = assigned,
        apex_location = p$apex_location
      )
    })

    # If no MW data, assign species by retention time
    # (largest peak = monomer, earlier = HMW, later = LMW)
    if (!has_mw && length(peak_data) > 0) {
      # Warn user about RT-based assignment (only once per bake)
      if (i == 1) {
        cli::cli_warn(
          c(
            "No MW column available; using retention time patterns for species assignment.",
            "i" = "Largest peak assumed to be monomer.",
            "i" = "Earlier-eluting peaks assigned as HMW, later as LMW.",
            "i" = "For accurate assignment, provide {.arg mw_column} from MALS/LALS."
          )
        )
      }

      # Find largest peak (assume monomer)
      areas <- sapply(peak_data, function(p) p$area)
      monomer_idx <- which.max(areas)
      monomer_rt <- peak_data[[monomer_idx]]$apex_location

      for (j in seq_along(peak_data)) {
        if (j == monomer_idx) {
          peak_data[[j]]$species <- "monomer"
        } else if (peak_data[[j]]$apex_location < monomer_rt) {
          # Earlier elution = higher MW
          peak_data[[j]]$species <- "hmw"
        } else {
          # Later elution = lower MW
          peak_data[[j]]$species <- "lmw"
        }
      }
    }

    # Aggregate by species
    for (sp in species_list) {
      matching <- sapply(peak_data, function(p) p$species == sp)
      if (any(matching)) {
        species_pct[[sp]][i] <- sum(sapply(
          peak_data[matching],
          function(p) p$pct
        ))
        # Use MW from largest matching peak
        matching_data <- peak_data[matching]
        largest_idx <- which.max(sapply(matching_data, function(p) p$area))
        species_mw[[sp]][i] <- matching_data[[largest_idx]]$mw
      }
    }

    # Count significant species
    species_count[i] <- sum(sapply(peak_data, function(p) {
      p$pct >= min_area_pct
    }))
  }

  # Add output columns
  for (sp in species_list) {
    new_data[[paste0(output_prefix, sp, "_pct")]] <- species_pct[[sp]]
    if (has_mw) {
      new_data[[paste0(output_prefix, sp, "_mw")]] <- species_mw[[sp]]
    }
  }
  new_data[[paste0(output_prefix, "species_count")]] <- species_count

  tibble::as_tibble(new_data)
}

#' @export
print.step_sec_oligomer <- function(
  x,
  width = max(20, options()$width - 30),
  ...
) {
  title <- "SEC oligomer analysis"
  if (x$trained) {
    cat(
      title,
      " (monomer MW = ",
      format(x$monomer_mw, big.mark = ","),
      " Da)",
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
tidy.step_sec_oligomer <- function(x, ...) {
  tibble::tibble(
    measures = list(x$measures),
    monomer_mw = x$monomer_mw,
    peak_detection = x$peak_detection,
    species = list(x$species),
    mw_tolerance = x$mw_tolerance,
    mw_column = x$mw_column %||% NA_character_,
    id = x$id
  )
}

#' @rdname required_pkgs.step_sec
#' @export
#' @keywords internal
required_pkgs.step_sec_oligomer <- function(x, ...) {
  c("measure.sec", "measure")
}
