# SEC/GPC Baseline Correction

`step_sec_baseline()` creates a *specification* of a recipe step that
applies baseline correction optimized for Gel Permeation Chromatography
(GPC) or Size Exclusion Chromatography (SEC) data. This method estimates
the baseline by interpolating between baseline regions at the start and
end of the chromatogram.

## Usage

``` r
step_sec_baseline(
  recipe,
  measures = NULL,
  left_frac = 0.05,
  right_frac = 0.05,
  method = "linear",
  role = NA,
  trained = FALSE,
  skip = FALSE,
  id = recipes::rand_id("sec_baseline")
)
```

## Arguments

- recipe:

  A recipe object. The step will be added to the sequence of operations
  for this recipe.

- measures:

  An optional character vector of measure column names to process. If
  `NULL` (the default), all measure columns (columns with class
  `measure_list`) will be processed.

- left_frac:

  Fraction of points from the beginning to use as the left baseline
  region. Default is `0.05` (first 5% of data points).

- right_frac:

  Fraction of points from the end to use as the right baseline region.
  Default is `0.05` (last 5% of data points).

- method:

  Method for baseline estimation. One of:

  - `"linear"` (default): Linear interpolation between left and right
    means

  - `"median"`: Uses median of baseline regions (more robust to
    outliers)

  - `"spline"`: Smooth spline through baseline regions

- role:

  Not used by this step since no new variables are created.

- trained:

  A logical to indicate if the quantities for preprocessing have been
  estimated.

- skip:

  A logical. Should the step be skipped when the recipe is baked?

- id:

  A character string that is unique to this step to identify it.

## Value

An updated version of `recipe` with the new step added to the sequence
of any existing operations.

## Details

GPC/SEC chromatograms typically have distinct baseline regions at the
beginning and end where no polymer elutes. This step leverages this
characteristic by:

1.  Identifying baseline regions at the start and end of the
    chromatogram

2.  Computing a representative baseline value for each region (mean or
    median)

3.  Interpolating between these values to estimate the full baseline

4.  Subtracting the estimated baseline from the signal

The `left_frac` and `right_frac` parameters control how much of the
chromatogram is considered "baseline". Choose values that:

- Include only the flat, signal-free regions

- Exclude any polymer peaks or system peaks

- Are large enough to average out noise

Unlike general-purpose baseline methods like ALS or polynomial fitting,
this approach is specifically designed for the characteristic shape of
GPC/SEC chromatograms and is computationally very fast.

**No selectors should be supplied to this step function**. The data
should be in the internal format produced by
[`measure::step_measure_input_wide()`](https://jameshwade.github.io/measure/reference/step_measure_input_wide.html)
or
[`measure::step_measure_input_long()`](https://jameshwade.github.io/measure/reference/step_measure_input_long.html).

## Tidying

When you
[`tidy()`](https://recipes.tidymodels.org/reference/tidy.recipe.html)
this step, a tibble with columns `terms`, `left_frac`, `right_frac`,
`method`, and `id` is returned.

## See also

[`measure::step_measure_baseline_als()`](https://jameshwade.github.io/measure/reference/step_measure_baseline_als.html)
for general-purpose baseline correction,
[`measure::step_measure_detrend()`](https://jameshwade.github.io/measure/reference/step_measure_detrend.html)
for simple trend removal.

Other sec-chromatography:
[`step_sec_band_broadening()`](step_sec_band_broadening.md),
[`step_sec_detector_delay()`](step_sec_detector_delay.md),
[`step_sec_mw_averages()`](step_sec_mw_averages.md),
[`step_sec_mw_distribution()`](step_sec_mw_distribution.md),
[`step_sec_mw_fractions()`](step_sec_mw_fractions.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(recipes)
library(measure)

# SEC baseline correction with default settings
rec <- recipe(~., data = sec_triple_detect) |>
  step_measure_input_wide(starts_with("signal_")) |>
  step_sec_baseline() |>
  prep()

# Using median method for robustness to outliers
rec <- recipe(~., data = sec_triple_detect) |>
  step_measure_input_wide(starts_with("signal_")) |>
  step_sec_baseline(left_frac = 0.1, right_frac = 0.1, method = "median") |>
  prep()
} # }
```
