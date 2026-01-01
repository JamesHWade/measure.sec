# Calculate Molecular Weight Averages for SEC/GPC

`step_sec_mw_averages()` creates a *specification* of a recipe step that
calculates molecular weight averages from size exclusion chromatography
data.

## Usage

``` r
step_sec_mw_averages(
  recipe,
  measures = NULL,
  calibration = NULL,
  integration_range = NULL,
  output_cols = c("mn", "mw", "mz", "mp", "dispersity"),
  prefix = "mw_",
  role = "predictor",
  trained = FALSE,
  skip = FALSE,
  id = recipes::rand_id("sec_mw_averages")
)
```

## Arguments

- recipe:

  A recipe object.

- measures:

  An optional character vector of measure column names.

- calibration:

  Calibration method for converting x-axis to log(MW). Can be:

  - `NULL` (default): Assumes x-axis is already log10(MW)

  - A numeric vector of length 2: Linear calibration
    `c(slope, intercept)` where `log10(MW) = slope * x + intercept`

  - `"auto"`: Estimate from data range (assumes typical polymer range)

- integration_range:

  Optional numeric vector `c(min, max)` specifying the x-axis range for
  integration. If `NULL`, uses full range.

- output_cols:

  Character vector of metrics to calculate. Default includes all:
  `c("mn", "mw", "mz", "mp", "dispersity")`.

- prefix:

  Prefix for output column names. Default is `"mw_"`.

- role:

  Role for generated columns. Default is `"predictor"`.

- trained:

  Logical indicating if the step has been trained.

- skip:

  Logical. Should the step be skipped when baking?

- id:

  Unique step identifier.

## Value

An updated recipe with the new step added.

## Details

This step calculates standard molecular weight averages from SEC/GPC
data:

|        |                     |                                   |
|--------|---------------------|-----------------------------------|
| Metric | Formula             | Description                       |
| Mn     | sum(w) / sum(w/M)   | Number-average molecular weight   |
| Mw     | sum(wM) / sum(w)    | Weight-average molecular weight   |
| Mz     | sum(wM^2) / sum(wM) | Z-average molecular weight        |
| Mp     | M at peak maximum   | Peak molecular weight             |
| D      | Mw/Mn               | Dispersity (polydispersity index) |

The detector signal is assumed to be proportional to weight
concentration. For RI detection, this is typically valid. For UV
detection, response factors may need to be applied first.

**Prerequisites:**

- Data should be baseline corrected

- X-axis should represent retention time/volume or log(MW)

- Integration limits should exclude solvent peaks

## See also

Other sec-chromatography: [`step_sec_baseline()`](step_sec_baseline.md),
[`step_sec_detector_delay()`](step_sec_detector_delay.md),
[`step_sec_mw_distribution()`](step_sec_mw_distribution.md),
[`step_sec_mw_fractions()`](step_sec_mw_fractions.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(recipes)
library(measure)

# Assuming x-axis is already calibrated to log10(MW)
rec <- recipe(~., data = sec_triple_detect) |>
  step_measure_input_wide(starts_with("signal_")) |>
  step_sec_baseline() |>
  step_sec_mw_averages() |>
  prep()
} # }
```
