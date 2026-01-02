# Create SEC Results Object

Constructor for the `sec_results` class, which wraps processed SEC/GPC
data and enables ggplot2's
[`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
functionality for automatic visualization.

## Usage

``` r
sec_results(data, sample_id = NULL)
```

## Arguments

- data:

  A data frame containing processed SEC results with measure columns.
  Typically the output from
  [`bake()`](https://recipes.tidymodels.org/reference/bake.html) on a
  prepped SEC recipe.

- sample_id:

  Optional. Column name containing sample identifiers. If `NULL`,
  auto-detection is attempted.

## Value

An object of class `sec_results` (inherits from `tbl_df`).

## Details

The `sec_results` class provides a unified interface for SEC/GPC data
that enables:

- Automatic plot selection via
  [`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)

- Integration with ggplot2 theming

- Summary statistics access

**Expected Data Structure:** The input data should contain measure
columns (list columns with `location` and `value` components). Common
measure columns include:

- `ri`, `uv`, `mals` - Detector signals

- `mw` - Molecular weight from calibration

- `concentration` - Concentration profile

- `intrinsic_visc` - Intrinsic viscosity

- `rg` - Radius of gyration

## See also

Other sec-visualization:
[`autoplot.sec_results()`](https://jameshwade.github.io/measure-sec/reference/autoplot.sec_results.md),
[`plot_sec()`](https://jameshwade.github.io/measure-sec/reference/plot_sec.md),
[`plot_sec_calibration()`](https://jameshwade.github.io/measure-sec/reference/plot_sec_calibration.md),
[`plot_sec_chromatogram()`](https://jameshwade.github.io/measure-sec/reference/plot_sec_chromatogram.md),
[`plot_sec_conformation()`](https://jameshwade.github.io/measure-sec/reference/plot_sec_conformation.md),
[`plot_sec_multidetector()`](https://jameshwade.github.io/measure-sec/reference/plot_sec_multidetector.md),
[`plot_sec_mwd()`](https://jameshwade.github.io/measure-sec/reference/plot_sec_mwd.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(recipes)
library(measure)
library(ggplot2)

# Process SEC data
processed <- recipe(~., data = sec_triple_detect) |>
  step_measure_input_long(ri_signal, location = vars(elution_time), col_name = "ri") |>
  step_sec_baseline(measures = "ri") |>
  step_sec_conventional_cal(standards = ps_standards) |>
  prep() |>
  bake(new_data = NULL)

# Wrap as sec_results
results <- sec_results(processed, sample_id = "sample_id")

# Use autoplot for automatic visualization
autoplot(results)
autoplot(results, type = "mwd")
autoplot(results, type = "chromatogram", normalize = TRUE)
} # }
```
