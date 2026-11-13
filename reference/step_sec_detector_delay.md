# Correct Inter-Detector Elution Delays

`step_sec_detector_delay()` creates a specification of a recipe step
that aligns detector signals to a reference detector in multi-detector
SEC systems.

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
  id = recipes::rand_id("sec_detector_delay"),
  elution_units = "min"
)
```

## Arguments

- recipe:

  A recipe object.

- reference:

  Character name of the reference detector column. The reference is
  unchanged and cannot also be a target.

- targets:

  Character vector of detector columns to shift. If `NULL`, all measure
  columns except the reference are shifted.
  [`character()`](https://rdrr.io/r/base/character.html) selects no
  target columns.

- delay_volumes:

  Named numeric vector of delay volumes in mL. Names must identify
  target columns. Positive values indicate that the target detector sees
  the sample after the reference.

- delay_times:

  Named numeric vector of delay times in minutes. This is an alternative
  to `delay_volumes`; the values remain in minutes in the step object.

- flow_rate:

  Flow rate in mL/min. Required only when converting between time and
  volume axes. Any supplied value must be positive and finite.

- method:

  Method for shifting signals: `"shift"` (index shift with round-to-even
  quantization) or `"interpolate"` (piecewise-linear interpolation).

- role:

  Not used by this step.

- trained:

  Logical indicating if the step has been trained.

- skip:

  Logical. Should the step be skipped when baking?

- id:

  Unique step identifier.

- elution_units:

  Units of the measure `location` coordinates: `"min"` (minutes), `"s"`
  (seconds), `"mL"` (millilitres), or `"uL"` (microlitres). The default
  is `"min"`.

## Value

An updated recipe with the new step added.

## Details

Detector delay declarations retain their own dimensions: `delay_times`
are always minutes and `delay_volumes` are always millilitres,
regardless of `elution_units`. The step converts the declared delay to
the location axis units during baking. Time-to-time and volume-to-volume
conversions do not need a flow rate. Converting between time and volume
needs a positive finite `flow_rate` in mL/min.

A corrected target value at location `x` is evaluated at `x + delay`.
Thus, a positive delay moves target values toward lower locations. The
reference and every target curve in a current bake row must have exactly
the same numeric location grid. This step does not resample mismatched
grids.

The `"shift"` method requires a nonzero-delay grid to be uniformly
spaced; adjacent widths may differ by at most 1e-8 relative to the first
width. Delays are rounded to the nearest whole index using R's
round-to-even rule. A zero delay is an exact pass-through and does not
require uniform spacing. The `"interpolate"` method accepts irregular
monotonic grids, returns `NA` outside the sampled range, and propagates
missing values without bridging an internal `NA` gap. Neither method
changes the location coordinates.

`elution_units` is a caller declaration. This step does not infer units
from record metadata or add an axis-unit attribute to measure columns.

## See also

Other sec-chromatography:
[`step_sec_band_broadening()`](https://jameshwade.github.io/measure.sec/reference/step_sec_band_broadening.md),
[`step_sec_baseline()`](https://jameshwade.github.io/measure.sec/reference/step_sec_baseline.md),
[`step_sec_exclude_regions()`](https://jameshwade.github.io/measure.sec/reference/step_sec_exclude_regions.md),
[`step_sec_integration_window()`](https://jameshwade.github.io/measure.sec/reference/step_sec_integration_window.md),
[`step_sec_mw_averages()`](https://jameshwade.github.io/measure.sec/reference/step_sec_mw_averages.md),
[`step_sec_mw_distribution()`](https://jameshwade.github.io/measure.sec/reference/step_sec_mw_distribution.md),
[`step_sec_mw_fractions()`](https://jameshwade.github.io/measure.sec/reference/step_sec_mw_fractions.md),
[`step_sec_peaks_deconvolve()`](https://jameshwade.github.io/measure.sec/reference/step_sec_peaks_deconvolve.md),
[`step_sec_peaks_detect()`](https://jameshwade.github.io/measure.sec/reference/step_sec_peaks_detect.md),
[`step_sec_peaks_refine()`](https://jameshwade.github.io/measure.sec/reference/step_sec_peaks_refine.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(recipes)
library(measure)

# Correct UV and MALS signals relative to RI on a minute axis
rec <- recipe(~., data = sec_triple_detect) |>
  step_measure_input_long(ri_signal, location = vars(elution_time), col_name = "ri") |>
  step_measure_input_long(uv_signal, location = vars(elution_time), col_name = "uv") |>
  step_measure_input_long(mals_signal, location = vars(elution_time), col_name = "mals") |>
  step_sec_detector_delay(
    reference = "ri",
    delay_volumes = c(uv = -0.05, mals = 0.15),
    elution_units = "min"
  ) |>
  prep()
} # }
```
