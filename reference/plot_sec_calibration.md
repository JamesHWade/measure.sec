# Plot SEC Calibration Curve

Creates a calibration curve plot from SEC calibration data or a prepped
recipe with conventional calibration.

## Usage

``` r
plot_sec_calibration(
  data,
  retention_col = "retention_time",
  mw_col = "log_mp",
  show_residuals = FALSE,
  show_r_squared = TRUE,
  ...
)
```

## Arguments

- data:

  Either:

  - A data frame of calibration standards with retention and log_mw
    columns

  - A prepped recipe containing a conventional calibration step

- retention_col:

  Name of retention time/volume column. Default is `"retention_time"`.

- mw_col:

  Name of log MW column. Default is `"log_mp"`.

- show_residuals:

  Logical. Show residual plot below? Default is `FALSE`.

- show_r_squared:

  Logical. Show R-squared value? Default is `TRUE`.

- ...:

  Additional arguments passed to
  [`ggplot2::geom_point()`](https://ggplot2.tidyverse.org/reference/geom_point.html).

## Value

A ggplot2 object (or patchwork of plots if show_residuals = TRUE).

## See also

Other sec-visualization:
[`plot_sec()`](https://jameshwade.github.io/measure-sec/reference/plot_sec.md),
[`plot_sec_chromatogram()`](https://jameshwade.github.io/measure-sec/reference/plot_sec_chromatogram.md),
[`plot_sec_conformation()`](https://jameshwade.github.io/measure-sec/reference/plot_sec_conformation.md),
[`plot_sec_multidetector()`](https://jameshwade.github.io/measure-sec/reference/plot_sec_multidetector.md),
[`plot_sec_mwd()`](https://jameshwade.github.io/measure-sec/reference/plot_sec_mwd.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Plot calibration from standards data
plot_sec_calibration(sec_ps_standards)

# With residuals panel
plot_sec_calibration(sec_ps_standards, show_residuals = TRUE)
} # }
```
