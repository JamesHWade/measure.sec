# Automatic Plot for SEC Results

Creates a ggplot2 visualization appropriate for SEC/GPC analysis
results. Automatically selects the best plot type based on available
data, or allows explicit selection via the `type` argument.

## Usage

``` r
# S3 method for class 'sec_results'
autoplot(
  object,
  type = c("auto", "chromatogram", "mwd", "conformation", "composition"),
  overlay_mw = TRUE,
  detectors = c("ri", "uv", "mals"),
  log_scale = c("x", "none", "y", "both"),
  ...
)
```

## Arguments

- object:

  An `sec_results` object created by
  [`sec_results()`](https://jameshwade.github.io/measure-sec/reference/sec_results.md).

- type:

  Type of plot to create. One of:

  - `"auto"`: Automatically detect best plot type (default)

  - `"chromatogram"`: Basic chromatogram (signal vs time)

  - `"mwd"`: Molecular weight distribution

  - `"conformation"`: Rg-MW or eta-MW scaling plot

  - `"composition"`: UV/RI ratio or composition plot

- overlay_mw:

  Logical. For chromatogram plots, overlay molecular weight on secondary
  y-axis? Default is `TRUE` when MW data is available.

- detectors:

  Character vector of detector columns to plot for multi-detector
  overlays. Default is `c("ri", "uv", "mals")`.

- log_scale:

  Character. Apply log scale to axes. Options:

  - `"x"`: Log scale on x-axis (default for MWD plots)

  - `"y"`: Log scale on y-axis

  - `"both"`: Log scale on both axes

  - `"none"`: No log scaling

- ...:

  Additional arguments passed to the underlying plot function.

## Value

A ggplot2 object.

## Details

When `type = "auto"` (default), the plot type is selected based on
available data:

1.  If `mw` column present: MWD plot

2.  If multiple detectors: Multi-detector overlay

3.  Otherwise: Basic chromatogram

The resulting ggplot2 object can be further customized with standard
ggplot2 functions like `+ theme_bw()` or `+ labs()`.

## See also

Other sec-visualization:
[`plot_sec()`](https://jameshwade.github.io/measure-sec/reference/plot_sec.md),
[`plot_sec_calibration()`](https://jameshwade.github.io/measure-sec/reference/plot_sec_calibration.md),
[`plot_sec_chromatogram()`](https://jameshwade.github.io/measure-sec/reference/plot_sec_chromatogram.md),
[`plot_sec_conformation()`](https://jameshwade.github.io/measure-sec/reference/plot_sec_conformation.md),
[`plot_sec_multidetector()`](https://jameshwade.github.io/measure-sec/reference/plot_sec_multidetector.md),
[`plot_sec_mwd()`](https://jameshwade.github.io/measure-sec/reference/plot_sec_mwd.md),
[`sec_results()`](https://jameshwade.github.io/measure-sec/reference/sec_results.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(ggplot2)

# Create sec_results object
results <- sec_results(processed_sec_data)

# Auto-detect best plot type
autoplot(results)

# Specific plot types
autoplot(results, type = "chromatogram")
autoplot(results, type = "mwd", show_averages = TRUE)
autoplot(results, type = "conformation")

# Customize with ggplot2
autoplot(results, type = "mwd") +
  theme_classic() +
  labs(title = "Molecular Weight Distribution")
} # }
```
