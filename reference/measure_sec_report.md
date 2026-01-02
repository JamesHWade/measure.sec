# Generate SEC Analysis Report

Creates a comprehensive report from SEC analysis results using Quarto
templates. Reports can be generated in HTML, PDF, or Word format.

## Usage

``` r
measure_sec_report(
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
)
```

## Arguments

- data:

  A data frame containing SEC results with measure columns, typically
  from a baked recipe.

- template:

  Report template to use: `"standard"` (default), `"detailed"`, or
  `"qc"`. See Details.

- output_format:

  Output format: `"html"` (default), `"pdf"`, or `"docx"`.

- output_file:

  Path for the output file. If `NULL` (default), creates a file in the
  current working directory with an auto-generated name.

- title:

  Report title. Default varies by template.

- author:

  Report author. Default is empty.

- sample_id:

  Column name containing sample identifiers.

- include_plots:

  Logical. Include visualizations? Default is `TRUE`.

- include_slice_table:

  Logical. Include slice-by-slice data table?

  Default is `FALSE`. Only used with `"detailed"` template.

- specs:

  Named list of QC specifications for the `"qc"` template. See Details.

- open:

  Logical. Open the report after generation? Default is `TRUE` for
  interactive sessions.

- quiet:

  Logical. Suppress Quarto output messages? Default is `FALSE`.

## Value

Invisibly returns the path to the generated report file.

## Details

### Templates

Three report templates are available:

**`"standard"`** (default):

- Molecular weight summary table

- Chromatogram overlay

- Molecular weight distribution plot

**`"detailed"`**:

- All content from standard template

- Multi-detector overlay

- Conformation plot (if MALS data available)

- Optional slice data table

- Analysis metadata

**`"qc"`** (Quality Control):

- System suitability metrics with pass/fail status

- Plate count, asymmetry, and resolution tests

- Calibration verification

- Acceptance criteria summary

### QC Specifications

For the `"qc"` template, you can provide custom specifications:

    specs = list(
      plate_count_min = 10000,
      asymmetry_min = 0.8,
      asymmetry_max = 1.5,
      resolution_min = 1.5,
      recovery_min = 95,
      recovery_max = 105
    )

### Requirements

- Quarto must be installed on your system

- The `quarto` R package must be installed

- For PDF output, a LaTeX distribution is required

## See also

Other sec-export:
[`measure_sec_compare()`](https://jameshwade.github.io/measure-sec/reference/measure_sec_compare.md),
[`measure_sec_slice_table()`](https://jameshwade.github.io/measure-sec/reference/measure_sec_slice_table.md),
[`measure_sec_summary_table()`](https://jameshwade.github.io/measure-sec/reference/measure_sec_summary_table.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(recipes)
library(measure)

# Process SEC data
processed <- recipe(~., data = sec_triple_detect) |>
  step_measure_input_long(ri, location = vars(time), col_name = "ri") |>
  step_sec_baseline() |>
  step_sec_conventional_cal(standards = ps_standards) |>
  step_sec_mw_averages() |>
  prep() |>
  bake(new_data = NULL)

# Generate standard HTML report
measure_sec_report(processed)

# Generate detailed PDF report
measure_sec_report(
  processed,
  template = "detailed",
  output_format = "pdf",
  title = "Polymer Characterization Results",
  author = "Lab Analyst"
)

# Generate QC report with custom specifications
measure_sec_report(
  sec_system_suitability,
  template = "qc",
  specs = list(plate_count_min = 15000)
)
} # }
```
