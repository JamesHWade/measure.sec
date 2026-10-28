# Generate Molecular Weight Distribution Curve

`step_sec_mw_distribution()` creates a *specification* of a recipe step
that generates molecular weight distribution curves from SEC/GPC data.

## Usage

``` r
step_sec_mw_distribution(
  recipe,
  measures = NULL,
  type = c("differential", "cumulative", "both"),
  calibration = NULL,
  n_points = 100L,
  mw_range = NULL,
  normalize = TRUE,
  role = NA,
  trained = FALSE,
  skip = FALSE,
  id = recipes::rand_id("sec_mw_distribution"),
  mw_col = NULL,
  mw_scale = NULL,
  integration_range = NULL
)
```

## Arguments

- recipe:

  A recipe object.

- measures:

  An optional character vector of measure column names.

- type:

  Type of distribution to generate:

  - `"differential"` (default): dW/d(log M) differential distribution

  - `"cumulative"`: Cumulative weight fraction distribution

  - `"both"`: Generate both distributions

- calibration:

  Calibration method for converting x-axis to log(MW). See
  [`step_sec_mw_averages()`](https://jameshwade.github.io/measure.sec/reference/step_sec_mw_averages.md)
  for details.

- n_points:

  Number of points in the output distribution. Default is 100. If
  `NULL`, uses the original data resolution.

- mw_range:

  Optional numeric vector `c(min, max)` specifying the MW range for the
  output distribution. If `NULL`, uses the range from data.

- normalize:

  Logical. Should the differential distribution be normalized to
  integrate to 1? Default is `TRUE`.

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

- integration_range:

  Optional limits on the original elution axis.

## Value

An updated recipe with the new step added.

## Details

This step transforms the raw chromatogram into standard MW distribution
representations:

**Differential Distribution (dW/d(log M)):** The weight fraction per
unit log(MW). This representation is preferred because the area under
the curve represents the weight fraction in that MW range.

**Cumulative Distribution:** The cumulative weight fraction from low to
high MW. Values range from 0 to 1. Weight is integrated on the original
elution axis before transforming to log(MW), so a nonlinear calibration
is handled without treating the signal as a density on the
molecular-weight axis. The selected MW curve must be strictly monotonic
and valid over the integration range. When `mw_range` restricts the
plotted range, it displays a portion of the full distribution; it does
not renormalize that portion to 1.

The selected measure column is replaced by the requested distribution.
With `type = "both"`, two additional columns, `<measure>_differential`
and `<measure>_cumulative`, are added and the selected column is
replaced by the differential curve. Their `location` values are
log10(MW), and their `value` entries are the corresponding distribution
values.

## See also

Other sec-chromatography:
[`step_sec_band_broadening()`](https://jameshwade.github.io/measure.sec/reference/step_sec_band_broadening.md),
[`step_sec_baseline()`](https://jameshwade.github.io/measure.sec/reference/step_sec_baseline.md),
[`step_sec_detector_delay()`](https://jameshwade.github.io/measure.sec/reference/step_sec_detector_delay.md),
[`step_sec_exclude_regions()`](https://jameshwade.github.io/measure.sec/reference/step_sec_exclude_regions.md),
[`step_sec_integration_window()`](https://jameshwade.github.io/measure.sec/reference/step_sec_integration_window.md),
[`step_sec_mw_averages()`](https://jameshwade.github.io/measure.sec/reference/step_sec_mw_averages.md),
[`step_sec_mw_fractions()`](https://jameshwade.github.io/measure.sec/reference/step_sec_mw_fractions.md),
[`step_sec_peaks_deconvolve()`](https://jameshwade.github.io/measure.sec/reference/step_sec_peaks_deconvolve.md),
[`step_sec_peaks_detect()`](https://jameshwade.github.io/measure.sec/reference/step_sec_peaks_detect.md),
[`step_sec_peaks_refine()`](https://jameshwade.github.io/measure.sec/reference/step_sec_peaks_refine.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(recipes)
library(measure)

# Generate differential MW distribution
rec <- recipe(~., data = sec_triple_detect) |>
  step_measure_input_wide(starts_with("signal_")) |>
  step_sec_baseline() |>
  step_sec_mw_distribution(type = "differential") |>
  prep()
} # }
```
