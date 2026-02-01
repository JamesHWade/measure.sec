# Refine SEC/GPC Peak Boundaries

`step_sec_peaks_refine()` creates a *specification* of a recipe step
that tightens peak boundaries set by
[`step_sec_peaks_detect()`](https://jameshwade.github.io/measure-sec/reference/step_sec_peaks_detect.md).
The default `"height_fraction"` method sets boundaries where the
corrected signal drops below a fraction of each peak's apex height,
matching the approach used by industry GPC software (e.g., Waters
Empower's 0.5% default).

## Usage

``` r
step_sec_peaks_refine(
  recipe,
  method = "height_fraction",
  cutoff = 0.005,
  peaks_col = ".peaks",
  measures_col = ".measures",
  role = NA,
  trained = FALSE,
  skip = FALSE,
  id = recipes::rand_id("sec_peaks_refine")
)
```

## Arguments

- recipe:

  A recipe object. The step will be added to the sequence of operations
  for this recipe.

- method:

  Boundary refinement method. Currently only `"height_fraction"` is
  supported:

  - `"height_fraction"` (default): Set boundaries where the corrected
    signal drops below `cutoff * apex_height`.

- cutoff:

  Numeric threshold for the chosen method. For `"height_fraction"`, this
  is the fraction of apex height below which the signal is considered
  outside the peak. Default is `0.005` (0.5%). Common values:

  - `0.005`: 0.5% of peak height (Empower default, most common)

  - `0.01`: 1% (more aggressive trimming)

  - `0.001`: 0.1% (looser, closer to baseline-return)

- peaks_col:

  Name of the peaks column to refine. Default is `".peaks"`.

- measures_col:

  Name of the measure column containing the baseline-corrected signal.
  Default is `".measures"`. If not found, the first available measure
  column is used automatically.

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
of any existing operations. The `.peaks` column will be updated with
refined `left_base` and `right_base` values. Two new columns are added
to each `peaks_tbl`: `original_left_base` and `original_right_base`
preserving the boundaries from
[`step_sec_peaks_detect()`](https://jameshwade.github.io/measure-sec/reference/step_sec_peaks_detect.md).

## Details

Wide peak boundaries from baseline-return detection cause systematic
errors in molecular weight calculations:

- **Mn bias**: Low-signal noise at high elution volumes contributes
  small M_i values, pulling Mn downward.

- **PDI inflation**: With Mw stable but Mn depressed, dispersity (Mw/Mn)
  increases artificially.

- **Noise integration**: Baseline noise fluctuating around zero gets
  included when `corrected > 0`.

The `"height_fraction"` method addresses this by walking inward from
each boundary until the signal exceeds a fraction of the peak's apex
height.

This step should run after
[`step_sec_baseline()`](https://jameshwade.github.io/measure-sec/reference/step_sec_baseline.md)
and
[`step_sec_peaks_detect()`](https://jameshwade.github.io/measure-sec/reference/step_sec_peaks_detect.md).
The signal in `measures_col` should be baseline-corrected for best
results.

**Edge cases handled:**

- Peak with apex height \<= 0: boundaries unchanged

- Fewer than 5 data points within boundaries: boundaries unchanged

- No points above threshold: boundaries unchanged

- Refined boundaries are guaranteed to never be wider than originals

## Tidying

When you
[`tidy()`](https://recipes.tidymodels.org/reference/tidy.recipe.html)
this step, a tibble with columns `method`, `cutoff`, `peaks_col`,
`measures_col`, and `id` is returned.

## See also

[`step_sec_peaks_detect()`](https://jameshwade.github.io/measure-sec/reference/step_sec_peaks_detect.md)
for peak detection,
[`step_sec_baseline()`](https://jameshwade.github.io/measure-sec/reference/step_sec_baseline.md)
for baseline correction.

Other sec-chromatography:
[`step_sec_band_broadening()`](https://jameshwade.github.io/measure-sec/reference/step_sec_band_broadening.md),
[`step_sec_baseline()`](https://jameshwade.github.io/measure-sec/reference/step_sec_baseline.md),
[`step_sec_detector_delay()`](https://jameshwade.github.io/measure-sec/reference/step_sec_detector_delay.md),
[`step_sec_exclude_regions()`](https://jameshwade.github.io/measure-sec/reference/step_sec_exclude_regions.md),
[`step_sec_integration_window()`](https://jameshwade.github.io/measure-sec/reference/step_sec_integration_window.md),
[`step_sec_mw_averages()`](https://jameshwade.github.io/measure-sec/reference/step_sec_mw_averages.md),
[`step_sec_mw_distribution()`](https://jameshwade.github.io/measure-sec/reference/step_sec_mw_distribution.md),
[`step_sec_mw_fractions()`](https://jameshwade.github.io/measure-sec/reference/step_sec_mw_fractions.md),
[`step_sec_peaks_deconvolve()`](https://jameshwade.github.io/measure-sec/reference/step_sec_peaks_deconvolve.md),
[`step_sec_peaks_detect()`](https://jameshwade.github.io/measure-sec/reference/step_sec_peaks_detect.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(recipes)
library(measure)

# Refine peaks with default 0.5% height fraction
rec <- recipe(~., data = sec_data) |>
  step_measure_input_long(ri_signal, location = vars(elution_volume)) |>
  step_sec_baseline() |>
  step_sec_peaks_detect() |>
  step_sec_peaks_refine() |>
  prep()

# More aggressive trimming at 1%
rec <- recipe(~., data = sec_data) |>
  step_measure_input_long(ri_signal, location = vars(elution_volume)) |>
  step_sec_baseline() |>
  step_sec_peaks_detect() |>
  step_sec_peaks_refine(cutoff = 0.01) |>
  prep()

# Looser boundaries at 0.1%
rec <- recipe(~., data = sec_data) |>
  step_measure_input_long(ri_signal, location = vars(elution_volume)) |>
  step_sec_baseline() |>
  step_sec_peaks_detect() |>
  step_sec_peaks_refine(cutoff = 0.001) |>
  prep()
} # }
```
