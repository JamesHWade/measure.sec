# Correct Inter-Detector Volume Delays

`step_sec_detector_delay()` creates a *specification* of a recipe step
that corrects for volume delays between detectors in multi-detector SEC
systems.

## Usage

``` r
step_sec_detector_delay(
  recipe,
  reference = NULL,
  targets = NULL,
  delay_volumes = NULL,
  delay_times = NULL,
  flow_rate = 1,
  method = c("shift", "interpolate"),
  role = NA,
  trained = FALSE,
  skip = FALSE,
  id = recipes::rand_id("sec_detector_delay")
)
```

## Arguments

- recipe:

  A recipe object.

- reference:

  Character name of the reference detector column (typically RI). All
  other detectors will be aligned to this reference.

- targets:

  Character vector of detector column names to shift. If `NULL`, all
  measure columns except the reference will be shifted.

- delay_volumes:

  Named numeric vector of delay volumes in mL. Names should match the
  target column names. Positive values indicate the detector sees the
  sample after the reference.

- delay_times:

  Named numeric vector of delay times in minutes. Alternative to
  `delay_volumes`. Requires `flow_rate` to be specified.

- flow_rate:

  Flow rate in mL/min. Required if using `delay_times`.

- method:

  Method for shifting signals:

  - `"shift"` (default): Simple index shift (fastest, slight edge
    effects)

  - `"interpolate"`: Linear interpolation (smoother, preserves signal
    shape)

- role:

  Not used by this step.

- trained:

  Logical indicating if the step has been trained.

- skip:

  Logical. Should the step be skipped when baking?

- id:

  Unique step identifier.

## Value

An updated recipe with the new step added.

## Details

In multi-detector SEC systems, detectors are connected in series and
separated by tubing. This causes each detector to see the same analyte
at different times. For accurate molecular weight calculations that
combine signals from multiple detectors (e.g., RI + MALS for absolute
MW), these delays must be corrected.

**Typical detector order and delays:**

- UV detector: Often first, minimal delay

- RI detector: Common reference detector

- MALS detector: Often has 0.1-0.3 mL delay from RI

- Viscometer: May have 0.2-0.5 mL delay

**Determining delay volumes:**

1.  Inject a narrow standard and record all detector signals

2.  Measure the time offset between peak maxima

3.  Convert to volume: delay_volume = time_offset × flow_rate

## See also

Other sec-chromatography:
[`step_sec_band_broadening()`](https://jameshwade.github.io/measure-sec/reference/step_sec_band_broadening.md),
[`step_sec_baseline()`](https://jameshwade.github.io/measure-sec/reference/step_sec_baseline.md),
[`step_sec_exclude_regions()`](https://jameshwade.github.io/measure-sec/reference/step_sec_exclude_regions.md),
[`step_sec_integration_window()`](https://jameshwade.github.io/measure-sec/reference/step_sec_integration_window.md),
[`step_sec_mw_averages()`](https://jameshwade.github.io/measure-sec/reference/step_sec_mw_averages.md),
[`step_sec_mw_distribution()`](https://jameshwade.github.io/measure-sec/reference/step_sec_mw_distribution.md),
[`step_sec_mw_fractions()`](https://jameshwade.github.io/measure-sec/reference/step_sec_mw_fractions.md),
[`step_sec_peaks_deconvolve()`](https://jameshwade.github.io/measure-sec/reference/step_sec_peaks_deconvolve.md),
[`step_sec_peaks_detect()`](https://jameshwade.github.io/measure-sec/reference/step_sec_peaks_detect.md),
[`step_sec_peaks_refine()`](https://jameshwade.github.io/measure-sec/reference/step_sec_peaks_refine.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(recipes)
library(measure)

# Correct UV and MALS signals relative to RI
rec <- recipe(~., data = sec_triple_detect) |>
  step_measure_input_long(ri_signal, location = vars(elution_time), col_name = "ri") |>
  step_measure_input_long(uv_signal, location = vars(elution_time), col_name = "uv") |>
  step_measure_input_long(mals_signal, location = vars(elution_time), col_name = "mals") |>
  step_sec_detector_delay(
    reference = "ri",
    delay_volumes = c(uv = -0.05, mals = 0.15)
  ) |>
  prep()
} # }
```
