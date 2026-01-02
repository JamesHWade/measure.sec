# Plot SEC Chromatogram

Creates a chromatogram plot from SEC data showing detector signal vs
elution time/volume.

## Usage

``` r
plot_sec_chromatogram(
  data,
  measures = NULL,
  sample_id = NULL,
  x_label = "Elution Time (min)",
  y_label = "Signal",
  normalize = FALSE,
  facet_by = c("none", "measure", "sample"),
  color_by = c("sample", "measure"),
  ...
)
```

## Arguments

- data:

  A data frame containing SEC results with measure columns, or a tibble
  from
  [`measure_sec_slice_table()`](https://jameshwade.github.io/measure-sec/reference/measure_sec_slice_table.md).

- measures:

  Character vector of measure column names to plot. If `NULL`, plots all
  measure columns found.

- sample_id:

  Column name containing sample identifiers. If `NULL`, attempts to
  auto-detect or uses row numbers.

- x_label:

  Label for x-axis. Default is "Elution Time (min)".

- y_label:

  Label for y-axis. Default is "Signal".

- normalize:

  Logical. Normalize signals to 0-1 range for comparison? Default is
  `FALSE`.

- facet_by:

  How to facet the plot. One of:

  - `"none"`: All on single plot (default)

  - `"measure"`: Separate panel per detector/measure

  - `"sample"`: Separate panel per sample

- color_by:

  What to map to color aesthetic. One of `"sample"` (default) or
  `"measure"`.

- ...:

  Additional arguments passed to
  [`ggplot2::geom_line()`](https://ggplot2.tidyverse.org/reference/geom_path.html).

## Value

A ggplot2 object.

## Details

This is the fundamental SEC visualization showing raw or processed
chromatographic data. Works with both:

- Processed recipe output (data frames with measure_list columns)

- Slice tables from
  [`measure_sec_slice_table()`](https://jameshwade.github.io/measure-sec/reference/measure_sec_slice_table.md)

## See also

Other sec-visualization:
[`autoplot.sec_results()`](https://jameshwade.github.io/measure-sec/reference/autoplot.sec_results.md),
[`plot_sec()`](https://jameshwade.github.io/measure-sec/reference/plot_sec.md),
[`plot_sec_calibration()`](https://jameshwade.github.io/measure-sec/reference/plot_sec_calibration.md),
[`plot_sec_composition()`](https://jameshwade.github.io/measure-sec/reference/plot_sec_composition.md),
[`plot_sec_conformation()`](https://jameshwade.github.io/measure-sec/reference/plot_sec_conformation.md),
[`plot_sec_multidetector()`](https://jameshwade.github.io/measure-sec/reference/plot_sec_multidetector.md),
[`plot_sec_mwd()`](https://jameshwade.github.io/measure-sec/reference/plot_sec_mwd.md),
[`sec_results()`](https://jameshwade.github.io/measure-sec/reference/sec_results.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(recipes)
library(measure)

# Process SEC data
processed <- recipe(~., data = sec_triple_detect) |>
  step_measure_input_long(
    ri_signal,
    location = vars(elution_time),
    col_name = "ri"
  ) |>
  step_sec_baseline(measures = "ri") |>
  prep() |>
  bake(new_data = NULL)

# Basic chromatogram
plot_sec_chromatogram(processed, measures = "ri")

# Normalized overlay of multiple samples
plot_sec_chromatogram(processed, measures = "ri", normalize = TRUE)

# Faceted by sample
plot_sec_chromatogram(processed, facet_by = "sample")
} # }
```
