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
  id = recipes::rand_id("sec_mw_fractions"),
  mw_col = NULL,
  mw_scale = NULL
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
  [`step_sec_mw_averages()`](https://jameshwade.github.io/measure.sec/reference/step_sec_mw_averages.md)
  for details.

- integration_range:

  Optional numeric vector `c(min, max)` specifying the x-axis range for
  integration. If `NULL`, uses a row-specific `.integration_window` when
  present, otherwise the full measured range.

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

- mw_col:

  Optional per-slice molecular-weight measure column. A matching curve
  produced by
  [`step_sec_conventional_cal()`](https://jameshwade.github.io/measure.sec/reference/step_sec_conventional_cal.md)
  is detected automatically.

- mw_scale:

  Scale of `mw_col`: `"log10"` or `"absolute"`. Required for untagged
  curves; see
  [`step_sec_mw_averages()`](https://jameshwade.github.io/measure.sec/reference/step_sec_mw_averages.md).

## Value

An updated recipe with the new step added.

## Details

For each cutoff value `C`, this step calculates:

- `frac_below_C`: Weight fraction with MW \< C

- `frac_above_C`: Weight fraction with MW \>= C

Slice weights use trapezoidal integration on the elution axis. Regions
in `.excluded_regions` marked `integration` or `both` are omitted and
the retained intervals are integrated separately; baseline-only
exclusions are ignored. Explicit `integration_range` values override
row-specific `.integration_window` metadata. The cutoff is applied to
the sampled molecular-weight values, so accuracy near a cutoff depends
on the sampling resolution. These fractions sum to 1.0 and are useful
for characterizing polymer distributions. Common cutoffs include:

- 1000 Da for oligomer content

- 10000 Da for low MW fraction

- 100000 Da for high MW fraction

## See also

Other sec-chromatography:
[`step_sec_band_broadening()`](https://jameshwade.github.io/measure.sec/reference/step_sec_band_broadening.md),
[`step_sec_baseline()`](https://jameshwade.github.io/measure.sec/reference/step_sec_baseline.md),
[`step_sec_detector_delay()`](https://jameshwade.github.io/measure.sec/reference/step_sec_detector_delay.md),
[`step_sec_exclude_regions()`](https://jameshwade.github.io/measure.sec/reference/step_sec_exclude_regions.md),
[`step_sec_integration_window()`](https://jameshwade.github.io/measure.sec/reference/step_sec_integration_window.md),
[`step_sec_mw_averages()`](https://jameshwade.github.io/measure.sec/reference/step_sec_mw_averages.md),
[`step_sec_mw_distribution()`](https://jameshwade.github.io/measure.sec/reference/step_sec_mw_distribution.md),
[`step_sec_peaks_deconvolve()`](https://jameshwade.github.io/measure.sec/reference/step_sec_peaks_deconvolve.md),
[`step_sec_peaks_detect()`](https://jameshwade.github.io/measure.sec/reference/step_sec_peaks_detect.md),
[`step_sec_peaks_refine()`](https://jameshwade.github.io/measure.sec/reference/step_sec_peaks_refine.md)

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
