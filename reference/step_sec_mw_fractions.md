# Calculate Molecular Weight Fractions for SEC/GPC

`step_sec_mw_fractions()` creates a *specification* of a recipe step
that calculates weight fractions above and below specified molecular
weight cutoffs.

## Usage

``` r
step_sec_mw_fractions(
  recipe,
  measures = NULL,
  cutoffs = c(1000, 10000, 1e+05),
  calibration = NULL,
  integration_range = NULL,
  prefix = "frac_",
  role = "predictor",
  trained = FALSE,
  skip = FALSE,
  id = recipes::rand_id("sec_mw_fractions")
)
```

## Arguments

- recipe:

  A recipe object.

- measures:

  An optional character vector of measure column names.

- cutoffs:

  Numeric vector of MW cutoff values. For each cutoff, the step
  calculates the weight fraction below and above that value.

- calibration:

  Calibration method for converting x-axis to log(MW). See
  [`step_sec_mw_averages()`](https://jameshwade.github.io/measure-sec/reference/step_sec_mw_averages.md)
  for details.

- integration_range:

  Optional numeric vector `c(min, max)` specifying the x-axis range for
  integration. If `NULL`, uses full range.

- prefix:

  Prefix for output column names. Default is `"frac_"`.

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

For each cutoff value `C`, this step calculates:

- `frac_below_C`: Weight fraction with MW \< C

- `frac_above_C`: Weight fraction with MW \>= C

These fractions sum to 1.0 and are useful for characterizing polymer
distributions. Common cutoffs include:

- 1000 Da for oligomer content

- 10000 Da for low MW fraction

- 100000 Da for high MW fraction

## See also

Other sec-chromatography:
[`step_sec_band_broadening()`](https://jameshwade.github.io/measure-sec/reference/step_sec_band_broadening.md),
[`step_sec_baseline()`](https://jameshwade.github.io/measure-sec/reference/step_sec_baseline.md),
[`step_sec_detector_delay()`](https://jameshwade.github.io/measure-sec/reference/step_sec_detector_delay.md),
[`step_sec_mw_averages()`](https://jameshwade.github.io/measure-sec/reference/step_sec_mw_averages.md),
[`step_sec_mw_distribution()`](https://jameshwade.github.io/measure-sec/reference/step_sec_mw_distribution.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(recipes)
library(measure)

# Calculate fractions at multiple cutoffs
rec <- recipe(~., data = sec_triple_detect) |>
  step_measure_input_wide(starts_with("signal_")) |>
  step_sec_baseline() |>
  step_sec_mw_fractions(cutoffs = c(1000, 10000, 100000)) |>
  prep()
} # }
```
