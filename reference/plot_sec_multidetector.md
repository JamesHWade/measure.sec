# Plot Multi-Detector SEC Overlay

Creates an overlay plot of multiple SEC detectors, optionally normalized
and aligned.

## Usage

``` r
plot_sec_multidetector(
  data,
  detectors,
  sample_id = NULL,
  samples = NULL,
  normalize = TRUE,
  x_label = "Elution Time (min)",
  facet = FALSE,
  ...
)
```

## Arguments

- data:

  A data frame containing SEC results with multiple detector measure
  columns.

- detectors:

  Character vector of detector column names to include. Common values:
  `c("ri", "uv", "mals", "visc")`.

- sample_id:

  Column name containing sample identifiers. If `NULL`, plots all
  samples or auto-detects.

- samples:

  Character vector of specific sample IDs to plot. If `NULL`, plots all
  samples.

- normalize:

  Logical. Normalize each detector to 0-1 range for comparison? Default
  is `TRUE`.

- x_label:

  Label for x-axis. Default is "Elution Time (min)".

- facet:

  Logical. Create separate panel for each sample? Default is `FALSE`
  (overlay).

- ...:

  Additional arguments passed to
  [`ggplot2::geom_line()`](https://ggplot2.tidyverse.org/reference/geom_path.html).

## Value

A ggplot2 object.

## Details

Multi-detector overlay plots are essential for:

- Verifying detector alignment after delay correction

- Identifying composition drift in copolymers (UV/RI differences)

- Detecting aggregates (MALS response higher than expected from RI)

- Checking for baseline issues across detectors

When `normalize = TRUE` (default), each detector signal is scaled to 0-1
range, making it easy to compare peak shapes and positions across
detectors with very different response magnitudes.

## See also

Other sec-visualization:
[`plot_sec()`](https://jameshwade.github.io/measure-sec/reference/plot_sec.md),
[`plot_sec_calibration()`](https://jameshwade.github.io/measure-sec/reference/plot_sec_calibration.md),
[`plot_sec_chromatogram()`](https://jameshwade.github.io/measure-sec/reference/plot_sec_chromatogram.md),
[`plot_sec_conformation()`](https://jameshwade.github.io/measure-sec/reference/plot_sec_conformation.md),
[`plot_sec_mwd()`](https://jameshwade.github.io/measure-sec/reference/plot_sec_mwd.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Multi-detector overlay for single sample
plot_sec_multidetector(
  processed_data,
  detectors = c("ri", "uv", "mals"),
  samples = "PMMA-1"
)

# Faceted by sample
plot_sec_multidetector(
  processed_data,
  detectors = c("ri", "uv"),
  facet = TRUE
)
} # }
```
