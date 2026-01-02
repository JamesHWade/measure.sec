# Extract Slice-by-Slice SEC Data

Extracts point-by-point (slice) data from SEC analysis results, creating
a long-format table suitable for export or further analysis.

## Usage

``` r
measure_sec_slice_table(
  data,
  measures = NULL,
  sample_id = NULL,
  include_location = TRUE,
  pivot = FALSE
)
```

## Arguments

- data:

  A data frame containing SEC results with measure columns.

- measures:

  Character vector of measure column names to include. If `NULL`,
  includes all measure columns found.

- sample_id:

  Column name containing sample identifiers. If `NULL`, uses row
  numbers.

- include_location:

  Logical. Include the location (time/volume) column? Default is `TRUE`.

- pivot:

  Logical. Pivot measures to wide format (one column per measure)?
  Default is `FALSE` (long format).

## Value

A tibble with slice-by-slice data:

- sample_id:

  Sample identifier

- slice:

  Slice index (1, 2, 3, ...)

- location:

  Elution time or volume

- measure:

  Measure column name (if pivot = FALSE)

- value:

  Signal value (if pivot = FALSE)

- \<measure_names\>:

  Individual measure columns (if pivot = TRUE)

## Details

This function extracts the raw slice data from processed SEC
chromatograms, making it easy to:

- Export to CSV/Excel for external analysis

- Create custom plots

- Perform slice-level calculations

- Compare samples point-by-point

**Typical Slice Data Columns:**

- Retention time/volume (location)

- Concentration (from RI or UV)

- Molecular weight (from calibration or MALS)

- Intrinsic viscosity (from viscometer)

- Radius of gyration (from MALS angles)

- Composition (from UV/RI ratio)

## See also

Other sec-export:
[`measure_sec_compare()`](https://jameshwade.github.io/measure-sec/reference/measure_sec_compare.md),
[`measure_sec_report()`](https://jameshwade.github.io/measure-sec/reference/measure_sec_report.md),
[`measure_sec_summary_table()`](https://jameshwade.github.io/measure-sec/reference/measure_sec_summary_table.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(recipes)
library(measure)

# Process SEC data
prepped <- recipe(~., data = sec_triple_detect) |>
  step_measure_input_long(ri, location = vars(time), col_name = "ri") |>
  step_sec_baseline() |>
  prep() |>
  bake(new_data = NULL)

# Extract slice table (long format)
slices <- measure_sec_slice_table(prepped, measures = "ri")

# Extract slice table (wide format)
slices_wide <- measure_sec_slice_table(
  prepped,
  measures = c("ri", "mw"),
  pivot = TRUE
)

# Export to CSV
write.csv(slices, "sec_slices.csv", row.names = FALSE)
} # }
```
