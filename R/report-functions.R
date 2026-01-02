# ==============================================================================
# SEC Report Generation Functions
#
# Quarto-based report generation for SEC analysis results
# ==============================================================================

#' Generate SEC Analysis Report
#'
#' Creates a comprehensive report from SEC analysis results using Quarto
#' templates. Reports can be generated in HTML, PDF, or Word format.
#'
#' @param data A data frame containing SEC results with measure columns,
#'   typically from a baked recipe.
#' @param template Report template to use: `"standard"` (default), `"detailed"`,
#'   or `"qc"`. See Details.
#' @param output_format Output format: `"html"` (default), `"pdf"`, or `"docx"`.
#' @param output_file Path for the output file. If `NULL` (default), creates
#'   a file in the current working directory with an auto-generated name.
#' @param title Report title. Default varies by template.
#' @param author Report author. Default is empty.
#' @param sample_id Column name containing sample identifiers.
#' @param include_plots Logical. Include visualizations? Default is `TRUE`.
#' @param include_slice_table Logical. Include slice-by-slice data table?
#'
#'   Default is `FALSE`. Only used with `"detailed"` template.
#' @param specs Named list of QC specifications for the `"qc"` template.
#'   See Details.
#' @param open Logical. Open the report after generation? Default is `TRUE`
#'   for interactive sessions.
#' @param quiet Logical. Suppress Quarto output messages? Default is `FALSE`.
#'
#' @return Invisibly returns the path to the generated report file.
#'
#' @details
#' ## Templates
#'
#' Three report templates are available:
#'
#' **`"standard"`** (default):
#' - Molecular weight summary table
#' - Chromatogram overlay
#' - Molecular weight distribution plot
#'
#' **`"detailed"`**:
#' - All content from standard template
#' - Multi-detector overlay
#' - Conformation plot (if MALS data available)
#' - Optional slice data table
#' - Analysis metadata
#'
#' **`"qc"`** (Quality Control):
#' - System suitability metrics with pass/fail status
#' - Plate count, asymmetry, and resolution tests
#' - Calibration verification
#' - Acceptance criteria summary
#'
#' ## QC Specifications
#'
#' For the `"qc"` template, you can provide custom specifications:
#' ```r
#' specs = list(
#'   plate_count_min = 10000,
#'   asymmetry_min = 0.8,
#'   asymmetry_max = 1.5,
#'   resolution_min = 1.5,
#'   recovery_min = 95,
#'   recovery_max = 105
#' )
#' ```
#'
#' ## Requirements
#'
#' - Quarto must be installed on your system
#' - The `quarto` R package must be installed
#' - For PDF output, a LaTeX distribution is required
#'
#' @family sec-export
#' @export
#'
#' @examples
#' \dontrun{
#' library(recipes)
#' library(measure)
#'
#' # Process SEC data
#' processed <- recipe(~., data = sec_triple_detect) |>
#'   step_measure_input_long(ri, location = vars(time), col_name = "ri") |>
#'   step_sec_baseline() |>
#'   step_sec_conventional_cal(standards = ps_standards) |>
#'   step_sec_mw_averages() |>
#'   prep() |>
#'   bake(new_data = NULL)
#'
#' # Generate standard HTML report
#' measure_sec_report(processed)
#'
#' # Generate detailed PDF report
#' measure_sec_report(
#'   processed,
#'   template = "detailed",
#'   output_format = "pdf",
#'   title = "Polymer Characterization Results",
#'   author = "Lab Analyst"
#' )
#'
#' # Generate QC report with custom specifications
#' measure_sec_report(
#'   sec_system_suitability,
#'   template = "qc",
#'   specs = list(plate_count_min = 15000)
#' )
#' }
measure_sec_report <- function(
  data,
  template = c("standard", "detailed", "qc"),
  output_format = c("html", "pdf", "docx"),
  output_file = NULL,
  title = NULL,
  author = "",
  sample_id = NULL,
  include_plots = TRUE,
  include_slice_table = FALSE,
  specs = NULL,
  open = interactive(),
  quiet = FALSE
) {
  # Check for quarto package
  if (!rlang::is_installed("quarto")) {
    cli::cli_abort(c(
      "Package {.pkg quarto} is required for report generation.",
      "i" = "Install it with {.code install.packages(\"quarto\")}"
    ))
  }

  # Check that Quarto CLI is available
  if (!nzchar(Sys.which("quarto"))) {
    cli::cli_abort(c(
      "Quarto CLI is not installed or not found in PATH.",
      "i" = "Install Quarto from {.url https://quarto.org/docs/get-started/}"
    ))
  }

  # Validate inputs
  if (!is.data.frame(data)) {
    cli::cli_abort("{.arg data} must be a data frame.")
  }

  template <- match.arg(template)
  output_format <- match.arg(output_format)

  # Set default title based on template

  if (is.null(title)) {
    title <- switch(
      template,
      standard = "SEC Analysis Report",
      detailed = "SEC Analysis Report - Detailed",
      qc = "SEC System Suitability Report"
    )
  }

  # Get template file path
  template_file <- get_sec_template(template)

  # Create output directory
  output_dir <- tempdir()

  # Copy template to temp directory for rendering
  temp_qmd <- file.path(output_dir, paste0("sec_report_", template, ".qmd"))
  if (!file.copy(template_file, temp_qmd, overwrite = TRUE)) {
    cli::cli_abort(c(
      "Failed to copy template to temporary directory.",
      "i" = "Template: {.file {template_file}}",
      "i" = "Destination: {.file {temp_qmd}}",
      "i" = "Check disk space and permissions."
    ))
  }

  # Determine output file name
  if (is.null(output_file)) {
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    output_file <- file.path(
      getwd(),
      paste0("sec_report_", template, "_", timestamp, ".", output_format)
    )
  }

  # Ensure output directory exists
  output_dir_final <- dirname(output_file)
  if (!dir.exists(output_dir_final)) {
    dir.create(output_dir_final, recursive = TRUE)
  }

  # Build execute_params
  execute_params <- list(
    title = title,
    author = author,
    data = data,
    sample_id = sample_id,
    include_plots = include_plots
  )

  # Add template-specific params
  if (template == "detailed") {
    execute_params$include_slice_table <- include_slice_table
    execute_params$include_summary <- TRUE
  } else if (template == "qc") {
    execute_params$specs <- specs
  } else {
    execute_params$include_summary <- TRUE
  }

  # Render the report
  cli::cli_progress_step("Generating {template} report...")

  tryCatch(
    {
      quarto::quarto_render(
        input = temp_qmd,
        output_format = output_format,
        output_file = basename(output_file),
        execute_params = execute_params,
        quiet = quiet
      )

      # Move output to final destination
      rendered_file <- file.path(
        output_dir,
        paste0("sec_report_", template, ".", output_format)
      )

      if (!file.exists(rendered_file)) {
        cli::cli_abort(c(
          "Report rendering succeeded but output file not found.",
          "i" = "Expected file: {.file {rendered_file}}",
          "i" = "This may indicate a Quarto version incompatibility."
        ))
      }

      if (!file.copy(rendered_file, output_file, overwrite = TRUE)) {
        cli::cli_abort(c(
          "Failed to copy rendered report to final destination.",
          "i" = "Source: {.file {rendered_file}}",
          "i" = "Destination: {.file {output_file}}",
          "i" = "Check disk space and permissions."
        ))
      }
      unlink(rendered_file)
    },
    error = function(e) {
      cli::cli_abort(c(
        "Report generation failed.",
        "i" = "Error: {e$message}",
        "i" = "Ensure Quarto is properly installed and the data format is correct."
      ))
    }
  )

  # Clean up temp files
  unlink(temp_qmd)

  cli::cli_alert_success("Report saved to {.file {output_file}}")

  # Open if requested
  if (open && file.exists(output_file)) {
    utils::browseURL(output_file)
  }

  invisible(output_file)
}


#' Get Path to SEC Report Template
#'
#' Returns the file path to a built-in SEC report template.
#'
#' @param template Template name: `"standard"`, `"detailed"`, or `"qc"`.
#'
#' @return Character string with the file path to the template.
#'
#' @keywords internal
#' @export
get_sec_template <- function(template = c("standard", "detailed", "qc")) {
  template <- match.arg(template)

  template_name <- paste0("sec_report_", template, ".qmd")

  template_path <- system.file(
    "templates",
    template_name,
    package = "measure.sec"
  )

  if (!nzchar(template_path) || !file.exists(template_path)) {
    cli::cli_abort(c(
      "Template {.val {template}} not found.",
      "i" = "Available templates: standard, detailed, qc"
    ))
  }

  template_path
}


#' List Available SEC Report Templates
#'
#' Lists all available report templates with descriptions.
#'
#' @return A tibble with template names and descriptions.
#'
#' @export
#'
#' @examples
#' list_sec_templates()
list_sec_templates <- function() {
  tibble::tibble(
    template = c("standard", "detailed", "qc"),
    description = c(
      "Summary table, chromatogram, and MWD plot",
      "All plots, multi-detector view, optional slice data",
      "System suitability with pass/fail metrics"
    ),
    formats = rep("html, pdf, docx", 3)
  )
}
