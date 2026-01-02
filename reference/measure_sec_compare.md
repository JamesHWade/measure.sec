# Compare Multiple SEC Samples

Compares SEC analysis results across multiple samples, providing a
summary table of key metrics and optional overlay plots.

## Usage

``` r
measure_sec_compare(
  ...,
  samples = NULL,
  metrics = c("mw_averages", "mwd", "branching"),
  plot = TRUE,
  reference = 1,
  digits = 2
)
```

## Arguments

- ...:

  Data frames containing SEC results. Each should be a baked recipe
  output with measure columns (MW, concentration, etc.).

- samples:

  Character vector of sample names. If `NULL`, uses names from `...` or
  generates sequential names ("Sample 1", "Sample 2", etc.).

- metrics:

  Character vector specifying which metrics to compare. Options:

  - `"mw_averages"`: Mn, Mw, Mz, dispersity

  - `"mwd"`: Molecular weight distribution overlap

  - `"branching"`: Branching metrics (if available)

  Default includes all available metrics.

- plot:

  Logical. Generate MWD comparison plot? Default is `TRUE`. Note: Plot
  is only generated when `"mwd"` is included in `metrics`.

- reference:

  Integer or character. Which sample to use as reference for percent
  differences. Default is `1` (first sample).

- digits:

  Integer. Number of decimal places for numeric values. Default is `2`.
  Must be non-negative.

## Value

A list of class `sec_comparison` containing:

- summary:

  Tibble with comparison metrics for all samples

- differences:

  Tibble with absolute and percent differences vs reference

- plot:

  ggplot2 object (if `plot = TRUE` and `"mwd"` in `metrics`), otherwise
  `NULL`

- samples:

  Character vector of sample names

- reference:

  Name of reference sample

## Details

This function is useful for:

- Batch-to-batch comparison

- Stability studies

- Process optimization

- Quality control

**Input Data Handling:**

When input data frames contain multiple rows (e.g., multiple
injections), numeric metrics are averaged across rows. For single-row
data frames (typical for processed SEC results), values are used
directly.

The function recognizes molecular weight columns with either naming
convention: prefixed (`mw_mn`, `mw_mw`, `mw_mz`, `mw_dispersity`) or
standard (`Mn`, `Mw`, `Mz`, `dispersity`).

**Comparison Metrics:**

*MW Averages:* Compares Mn, Mw, Mz, and dispersity with absolute and
percent differences from the reference sample.

*MWD Overlay:* Creates overlaid MWD plots showing distribution
differences. Useful for detecting bimodality, tailing, or distribution
shifts.

*Branching:* Compares branching index and frequency if available from
triple-detection data.

## See also

Other sec-export:
[`measure_sec_report()`](https://jameshwade.github.io/measure-sec/reference/measure_sec_report.md),
[`measure_sec_slice_table()`](https://jameshwade.github.io/measure-sec/reference/measure_sec_slice_table.md),
[`measure_sec_summary_table()`](https://jameshwade.github.io/measure-sec/reference/measure_sec_summary_table.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Compare three batches
comparison <- measure_sec_compare(
  batch1_data,
  batch2_data,
  batch3_data,
  samples = c("Batch 1", "Batch 2", "Batch 3")
)

# View summary table
comparison$summary

# View differences from reference
comparison$differences

# Display overlay plot
comparison$plot

# Compare only MW averages without plot
comparison <- measure_sec_compare(
  sample_a, sample_b,
  metrics = "mw_averages",
  plot = FALSE
)
} # }
```
