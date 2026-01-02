# Quick SEC Plot

Creates an automatic plot appropriate for SEC analysis results.
Dispatches to the most appropriate plot type based on available data.

## Usage

``` r
plot_sec(
  data,
  type = c("auto", "chromatogram", "mwd", "multidetector", "conformation"),
  ...
)
```

## Arguments

- data:

  A data frame containing SEC results with measure columns.

- type:

  Type of plot to create. One of:

  - `"auto"`: Automatically detect (default)

  - `"chromatogram"`: Basic chromatogram

  - `"mwd"`: Molecular weight distribution

  - `"multidetector"`: Multi-detector overlay

  - `"conformation"`: Rg-MW or eta-MW plot

- ...:

  Additional arguments passed to the underlying plot function.

## Value

A ggplot2 object.

## Details

When `type = "auto"` (default), the function chooses the plot type based
on available data:

- If MW column present: MWD plot

- If multiple detectors: Multi-detector overlay

- Otherwise: Basic chromatogram

This is a convenience wrapper that dispatches to the specific plot
functions:
[`plot_sec_chromatogram()`](https://jameshwade.github.io/measure-sec/reference/plot_sec_chromatogram.md),
[`plot_sec_mwd()`](https://jameshwade.github.io/measure-sec/reference/plot_sec_mwd.md),
[`plot_sec_multidetector()`](https://jameshwade.github.io/measure-sec/reference/plot_sec_multidetector.md),
or
[`plot_sec_conformation()`](https://jameshwade.github.io/measure-sec/reference/plot_sec_conformation.md).

## See also

Other sec-visualization:
[`autoplot.sec_results()`](https://jameshwade.github.io/measure-sec/reference/autoplot.sec_results.md),
[`plot_sec_calibration()`](https://jameshwade.github.io/measure-sec/reference/plot_sec_calibration.md),
[`plot_sec_chromatogram()`](https://jameshwade.github.io/measure-sec/reference/plot_sec_chromatogram.md),
[`plot_sec_composition()`](https://jameshwade.github.io/measure-sec/reference/plot_sec_composition.md),
[`plot_sec_conformation()`](https://jameshwade.github.io/measure-sec/reference/plot_sec_conformation.md),
[`plot_sec_multidetector()`](https://jameshwade.github.io/measure-sec/reference/plot_sec_multidetector.md),
[`plot_sec_mwd()`](https://jameshwade.github.io/measure-sec/reference/plot_sec_mwd.md),
[`sec_results()`](https://jameshwade.github.io/measure-sec/reference/sec_results.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(ggplot2)

# Auto-detect plot type
plot_sec(processed_sec_data)

# Specific plot type
plot_sec(processed_sec_data, type = "mwd")
} # }
```
