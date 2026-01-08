# ==============================================================================
# standards-helpers.R
#
# Helper functions for matching and normalizing calibration standard names
# ==============================================================================

#' Get Recognized SEC Calibration Standards
#'
#' Returns a list of recognized calibration standards with their matching
#' patterns. This is the central registry for standard name recognition.
#'
#' @return A named list where names are canonical standard names and values
#'   are regex patterns for matching normalized sample names.
#'
#' @details
#' Standard types supported:
#' - **PS (Polystyrene)**: PS A, PS B, PS C, PS D, PS 1683
#' - **PMMA (Polymethyl methacrylate)**: PMMA A, PMMA B, PMMA C, PMMA D
#' - **PEG/PEO (Polyethylene glycol/oxide)**: PEG A, PEG B, PEG C, PEG D
#' - **Pullulan**: Pullulan A, Pullulan B, Pullulan C, Pullulan D
#' - **Generic**: Standard 1-4 (site-specific, no default MW values)
#'
#' Patterns assume input has been normalized via [normalize_standard_name()]
#' (lowercase, "std" converted to "standard", separators normalized).
#'
#' @examples
#' standards <- get_recognized_standards()
#' names(standards)
#'
#' @family standards
#' @export
get_recognized_standards <- function() {
  list(
    # Polystyrene standards (most common)
    "PS 1683" = "\\b(ps\\s*1683|1683)\\b",
    "PS A" = "\\bps\\s*a\\b",
    "PS B" = "\\bps\\s*b\\b",
    "PS C" = "\\bps\\s*c\\b",
    "PS D" = "\\bps\\s*d\\b",

    # PMMA standards
    "PMMA A" = "\\bpmma\\s*a\\b",
    "PMMA B" = "\\bpmma\\s*b\\b",
    "PMMA C" = "\\bpmma\\s*c\\b",
    "PMMA D" = "\\bpmma\\s*d\\b",

    # PEG/PEO standards (aqueous GPC)
    "PEG A" = "\\b(peg|peo)\\s*a\\b",
    "PEG B" = "\\b(peg|peo)\\s*b\\b",
    "PEG C" = "\\b(peg|peo)\\s*c\\b",
    "PEG D" = "\\b(peg|peo)\\s*d\\b",

    # Pullulan standards (aqueous GPC)
    "Pullulan A" = "\\bpullulan\\s*a\\b",
    "Pullulan B" = "\\bpullulan\\s*b\\b",
    "Pullulan C" = "\\bpullulan\\s*c\\b",
    "Pullulan D" = "\\bpullulan\\s*d\\b",

    # Generic standards (site-specific MW values)
    "Standard 1" = "\\bstandard\\s*1\\b",
    "Standard 2" = "\\bstandard\\s*2\\b",
    "Standard 3" = "\\bstandard\\s*3\\b",
    "Standard 4" = "\\bstandard\\s*4\\b"
  )
}

#' Get Standard Names
#'
#' Returns a character vector of all recognized standard names.
#'
#' @return Character vector of canonical standard names.
#'
#' @examples
#' get_standard_names()
#'
#' @family standards
#' @export
get_standard_names <- function() {
  names(get_recognized_standards())
}

#' Get Standard Type
#'
#' Classifies a canonical standard name into its polymer type.
#'
#' @param standard_name Character string of the standard name (canonical form).
#'
#' @return Character string of the standard type: "PS", "PMMA", "PEG",
#'   "Pullulan", "Generic", or "Other".
#'
#' @examples
#' get_standard_type("PS A")
#' get_standard_type("Standard 2")
#' get_standard_type("PMMA C")
#'
#' @family standards
#' @export
get_standard_type <- function(standard_name) {
  dplyr::case_when(
    grepl("^PS ", standard_name) ~ "PS",
    grepl("^PMMA ", standard_name) ~ "PMMA",
    grepl("^PEG ", standard_name) ~ "PEG",
    grepl("^Pullulan ", standard_name) ~ "Pullulan",
    grepl("^Standard ", standard_name) ~ "Generic",
    TRUE ~ "Other"
  )
}

#' Normalize Standard Names for Matching
#'
#' Normalizes sample and standard names to improve matching by:
#' - Converting to lowercase
#' - Replacing separators (-, _) with spaces
#' - Removing date prefixes (YYYYMMDD, YYYY-MM-DD, etc.)
#' - Removing common lab/vendor identifiers
#' - Normalizing "std" variations to "standard"
#' - Cleaning up whitespace
#'
#' @param name Character vector of names to normalize.
#'
#' @return Character vector of normalized names.
#'
#' @examples
#' normalize_standard_name("20231215_PS_A")
#' normalize_standard_name("Std B")
#' normalize_standard_name("Phillips PS-C")
#' normalize_standard_name(c("DOW_PMMA_A", "2023-01-15 StdD"))
#'
#' @family standards
#' @export
normalize_standard_name <- function(name) {
  if (length(name) == 0) {
    return(character(0))
  }

  # Process each name individually for robustness

  vapply(name, function(n) {
    if (is.na(n)) {
      return(NA_character_)
    }

    # Convert to lowercase and replace separators with spaces
    normalized <- tolower(n)
    normalized <- gsub("[-_]", " ", normalized)

    # Normalize all "std" variations to "standard"
    normalized <- gsub("\\bstd([0-9]+)", "standard \\1", normalized)
    normalized <- gsub("\\bstd([a-z]+)", "standard \\1", normalized)
    normalized <- gsub("\\bstd\\s+", "standard ", normalized)
    normalized <- gsub("\\bstd\\b", "standard", normalized)

    # Remove date patterns at the start
    normalized <- gsub("^\\d{8}\\s*", "", normalized)
    normalized <- gsub("^\\d{4} \\d{2} \\d{2}\\s*", "", normalized)
    normalized <- gsub("^\\d{2} \\d{2} \\d{2}\\s*", "", normalized)

    # Remove common lab/vendor identifiers
    normalized <- gsub("phillips\\s*", "", normalized, ignore.case = TRUE)
    normalized <- gsub("dow\\s*", "", normalized, ignore.case = TRUE)
    normalized <- gsub("agilent\\s*", "", normalized, ignore.case = TRUE)
    normalized <- gsub("waters\\s*", "", normalized, ignore.case = TRUE)

    # Remove alphanumeric batch identifiers at start (like ub08119)
    normalized <- gsub("^[a-z]{2}\\d+\\s*", "", normalized)

    # Clean up multiple spaces and trim
    normalized <- gsub("\\s+", " ", normalized)
    normalized <- trimws(normalized)

    normalized
  }, character(1), USE.NAMES = FALSE)
}

#' Match Sample Names to Recognized Standards
#'
#' Attempts to match sample names to recognized calibration standards using
#' pattern matching on normalized names.
#'
#' @param sample_names Character vector of sample names to match.
#' @param standards Optional named list of standards with regex patterns.
#'   If `NULL` (default), uses [get_recognized_standards()].
#' @param return_all Logical. If `TRUE`, returns all matches including
#'   non-matches as NA. If `FALSE` (default), returns only matched entries.
#'
#' @return A tibble with columns:
#'   - `sample_name`: Original sample name
#'   - `normalized_name`: Cleaned/normalized name
#'   - `matched_standard`: Canonical standard name (or NA if no match)
#'   - `standard_type`: Type of standard (PS, PMMA, etc.)
#'
#' @examples
#' samples <- c("20231215_PS_A", "Unknown Sample", "Std B", "PMMA-C")
#' match_standards(samples)
#'
#' # Return all samples including non-matches
#' match_standards(samples, return_all = TRUE)
#'
#' @family standards
#' @export
match_standards <- function(
    sample_names,
    standards = NULL,
    return_all = FALSE) {
  if (is.null(standards)) {
    standards <- get_recognized_standards()
  }

  # Handle empty input
  if (length(sample_names) == 0) {
    return(tibble::tibble(
      sample_name = character(0),
      normalized_name = character(0),
      matched_standard = character(0),
      standard_type = character(0)
    ))
  }

  # Normalize all sample names
  normalized <- normalize_standard_name(sample_names)

  # Try to match each sample to a standard
  matches <- purrr::map2(
    sample_names,
    normalized,
    function(orig, norm) {
      matched <- NA_character_

      for (std_name in names(standards)) {
        pattern <- standards[[std_name]]
        if (grepl(pattern, norm, perl = TRUE)) {
          matched <- std_name
          break
        }
      }

      tibble::tibble(
        sample_name = orig,
        normalized_name = norm,
        matched_standard = matched,
        standard_type = if (is.na(matched)) NA_character_ else get_standard_type(matched)
      )
    }
  )

  result <- dplyr::bind_rows(matches)

  if (!return_all) {
    result <- dplyr::filter(result, !is.na(.data$matched_standard))
  }

  result
}

#' Check if Sample Names Contain Standards
#'
#' Quick check to determine if any sample names match recognized standards.
#'
#' @param sample_names Character vector of sample names.
#'
#' @return Logical vector indicating whether each name matches a standard.
#'
#' @examples
#' is_standard(c("PS A", "Unknown", "PMMA-B"))
#'
#' @family standards
#' @export
is_standard <- function(sample_names) {
  matches <- match_standards(sample_names, return_all = TRUE)
  !is.na(matches$matched_standard)
}

#' Auto-Detect Standards in a Data Frame
#'
#' Scans a data frame for columns that might contain standard names and
#' attempts to match them to recognized standards.
#'
#' @param data A data frame containing sample information.
#' @param name_col Name of the column containing sample names. If `NULL`,
#'   attempts to auto-detect from common column names.
#' @param add_columns Logical. If `TRUE`, adds `matched_standard` and
#'   `standard_type` columns to the data frame. If `FALSE`, returns only
#'   the matching results.
#'
#' @return If `add_columns = TRUE`, the original data frame with added
#'   columns. If `FALSE`, a tibble with matching results.
#'
#' @examples
#' \dontrun{
#' # Auto-detect and add standard columns
#' data <- data.frame(
#'   sample_name = c("PS A", "Unknown", "PMMA B"),
#'   signal = c(100, 200, 150)
#' )
#' auto_detect_standards(data)
#' }
#'
#' @family standards
#' @export
auto_detect_standards <- function(
    data,
    name_col = NULL,
    add_columns = TRUE) {
  # Try to find the name column
  if (is.null(name_col)) {
    possible_cols <- c(
      "sample_name", "name", "sample", "id", "sample_id",
      "Name", "Sample", "Sample_Name", "SampleName"
    )
    found_cols <- intersect(possible_cols, names(data))

    if (length(found_cols) == 0) {
      cli::cli_abort(
        c(
          "Could not auto-detect sample name column.",
          "i" = "Specify {.arg name_col} explicitly."
        )
      )
    }

    name_col <- found_cols[1]
    cli::cli_inform("Using column {.val {name_col}} for standard detection.")
  }

  if (!name_col %in% names(data)) {
    cli::cli_abort("Column {.val {name_col}} not found in data.")
  }

  sample_names <- data[[name_col]]
  matches <- match_standards(sample_names, return_all = TRUE)

  if (add_columns) {
    data$matched_standard <- matches$matched_standard
    data$standard_type <- matches$standard_type
    return(tibble::as_tibble(data))
  }

  matches
}
