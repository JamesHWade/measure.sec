# SEC/GPC Peak Detection

`step_sec_peaks_detect()` creates a *specification* of a recipe step
that detects peaks in SEC/GPC chromatography data. The default algorithm
is `finderskeepers`, which uses LOESS smoothing with Iterative Soft
Thresholding (IST) and changepoint analysis for automatic threshold
detection.

## Usage

``` r
step_sec_peaks_detect(
  recipe,
  measures = NULL,
  algorithm = "finderskeepers",
  min_height = 1,
  min_distance = 0,
  loess_span = 0.01,
  ist_points = 50L,
  ist_nonlinearity = 5,
  snr_threshold = FALSE,
  role = NA,
  trained = FALSE,
  skip = FALSE,
  id = recipes::rand_id("sec_peaks_detect")
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

- algorithm:

  Peak detection algorithm. Currently only `"finderskeepers"` (the
  default) is supported. This SEC-optimized algorithm uses LOESS
  smoothing with Iterative Soft Thresholding and changepoint detection.
  For other algorithms, use
  [`measure::step_measure_peaks_detect()`](https://jameshwade.github.io/measure/reference/step_measure_peaks_detect.html)
  directly.

- min_height:

  Minimum peak height. For `"finderskeepers"`, this is the minimum
  height above baseline. For other algorithms with
  `snr_threshold = TRUE`, this is interpreted as a signal-to-noise
  ratio. Default is `1`.

- min_distance:

  Minimum distance between peaks in location units (e.g., mL). Only used
  for non-finderskeepers algorithms. Default is `0`.

- loess_span:

  LOESS smoothing span for `finderskeepers`. A value between 0 and 1
  controlling the smoothness of the fit. Default is `0.01` (minimal
  smoothing to preserve peak shapes).

- ist_points:

  Number of threshold levels for Iterative Soft Thresholding. Default is
  `50`. Higher values give finer threshold resolution.

- ist_nonlinearity:

  Nonlinearity parameter for IST threshold spacing. Default is `5`.
  Higher values concentrate thresholds near the baseline.

- snr_threshold:

  Logical. For non-finderskeepers algorithms, if `TRUE`, `min_height` is
  interpreted as a signal-to-noise ratio. Default is `FALSE`.

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

An updated version of `recipe` with the new step added to the

sequence of any existing operations. A new `.peaks` column will be added
containing detected peaks for each sample.

## Details

The `finderskeepers` algorithm is specifically designed for SEC/GPC
data:

1 . **LOESS smoothing**: Applies local polynomial regression to reduce
noise while preserving peak shapes (controlled by `loess_span`) 2.
**Iterative Soft Thresholding (IST)**: Creates a series of thresholds
with nonlinear spacing to detect changes in peak structure 3.
**Changepoint detection**: Uses
[`changepoint::cpt.mean()`](https://rdrr.io/pkg/changepoint/man/cpt.mean.html)
to automatically determine the optimal threshold for peak/baseline
separation 4. **Peak boundary detection**: Identifies peak start, apex,
and end points

This approach is robust to baseline drift and varying peak heights,
making it well-suited for polymer SEC chromatograms.

**Peak properties stored:**

- `peak_id`: Integer identifier

- `location`: X-axis position of peak apex (elution volume)

- `height`: Peak height above baseline

- `left_base`, `right_base`: X-axis positions of peak boundaries

- `area`: Initially NA; use
  [`measure::step_measure_peaks_integrate()`](https://jameshwade.github.io/measure/reference/step_measure_peaks_integrate.html)
  to calculate

## Tidying

When you
[`tidy()`](https://recipes.tidymodels.org/reference/tidy.recipe.html)
this step, a tibble with columns `terms`, `algorithm`, `min_height`, and
`id` is returned.

## See also

[`measure::step_measure_peaks_detect()`](https://jameshwade.github.io/measure/reference/step_measure_peaks_detect.html)
for general peak detection,
[`measure::step_measure_peaks_integrate()`](https://jameshwade.github.io/measure/reference/step_measure_peaks_integrate.html)
to calculate peak areas.

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
[`step_sec_peaks_refine()`](https://jameshwade.github.io/measure-sec/reference/step_sec_peaks_refine.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(recipes)
library(measure)

# SEC peak detection with finderskeepers (default)
rec <- recipe(~., data = sec_triple_detect) |>
  step_measure_input_long(ri_signal, location = vars(elution_time)) |>
  step_sec_peaks_detect() |>
  prep()

# Adjust sensitivity with min_height
rec <- recipe(~., data = sec_triple_detect) |>
  step_measure_input_long(ri_signal, location = vars(elution_time)) |>
  step_sec_peaks_detect(min_height = 5) |>
  prep()

# Adjust LOESS smoothing for noisy data
rec <- recipe(~., data = sec_triple_detect) |>
  step_measure_input_long(ri_signal, location = vars(elution_time)) |>
  step_sec_peaks_detect(loess_span = 0.05) |>
  prep()
} # }
```
