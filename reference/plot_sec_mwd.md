# Plot Molecular Weight Distribution

Creates a molecular weight distribution (MWD) plot showing differential
or cumulative distribution.

## Usage

``` r
plot_sec_mwd(
  data,
  mw_col = "mw",
  concentration_col = NULL,
  sample_id = NULL,
  type = c("differential", "cumulative", "both"),
  show_averages = TRUE,
  log_mw = TRUE,
  x_label = NULL,
  y_label = NULL,
  ...
)
```

## Arguments

- data:

  A data frame containing SEC results with calibrated MW data.

- mw_col:

  Name of the molecular weight measure column. Default is `"mw"`.

- concentration_col:

  Name of the concentration/signal column used for weighting. If `NULL`,
  uses the first available measure column.

- sample_id:

  Column name containing sample identifiers.

- type:

  Type of distribution to plot:

  - `"differential"`: dW/d(log M) vs log M (default)

  - `"cumulative"`: Cumulative weight fraction vs log M

  - `"both"`: Both on faceted plot

- show_averages:

  Logical. Show vertical lines for Mn, Mw, Mz? Default is `TRUE`.
  Requires these columns to be present in data.

- log_mw:

  Logical. Use log10(MW) on x-axis? Default is `TRUE`.

- x_label:

  Label for x-axis. Default is auto-generated based on `log_mw`.

- y_label:

  Label for y-axis. Default is auto-generated based on `type`.

- ...:

  Additional arguments passed to
  [`ggplot2::geom_line()`](https://ggplot2.tidyverse.org/reference/geom_path.html).

## Value

A ggplot2 object.

## Details

The MWD plot is the standard way to visualize polymer molecular weight
distributions. The differential distribution (dW/d log M) shows the
weight fraction of polymer at each molecular weight.

When `show_averages = TRUE`, vertical dashed lines are added for:

- Mn (number-average): leftmost

- Mw (weight-average): middle

- Mz (z-average): rightmost

## See also

Other sec-visualization:
[`autoplot.sec_results()`](https://jameshwade.github.io/measure-sec/reference/autoplot.sec_results.md),
[`plot_sec()`](https://jameshwade.github.io/measure-sec/reference/plot_sec.md),
[`plot_sec_calibration()`](https://jameshwade.github.io/measure-sec/reference/plot_sec_calibration.md),
[`plot_sec_chromatogram()`](https://jameshwade.github.io/measure-sec/reference/plot_sec_chromatogram.md),
[`plot_sec_composition()`](https://jameshwade.github.io/measure-sec/reference/plot_sec_composition.md),
[`plot_sec_conformation()`](https://jameshwade.github.io/measure-sec/reference/plot_sec_conformation.md),
[`plot_sec_multidetector()`](https://jameshwade.github.io/measure-sec/reference/plot_sec_multidetector.md),
[`sec_results()`](https://jameshwade.github.io/measure-sec/reference/sec_results.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# After SEC processing with calibration
plot_sec_mwd(processed_data)

# Cumulative distribution
plot_sec_mwd(processed_data, type = "cumulative")

# Without MW average lines
plot_sec_mwd(processed_data, show_averages = FALSE)
} # }
```
