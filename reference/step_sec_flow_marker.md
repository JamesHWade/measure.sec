# Flow Marker Correction for SEC/GPC

`step_sec_flow_marker()` creates a *specification* of a recipe step that
detects a flow marker peak (typically toluene or other small molecule)
and applies a linear correction to align retention times/volumes across
runs.

## Usage

``` r
step_sec_flow_marker(
  recipe,
  measures = NULL,
  marker_range = NULL,
  target_volume = NULL,
  auto_detect = TRUE,
  min_peak_height = NULL,
  store_correction = TRUE,
  role = NA,
  trained = FALSE,
  skip = FALSE,
  id = recipes::rand_id("sec_flow_marker")
)
```

## Arguments

- recipe:

  A recipe object.

- measures:

  Character vector of measure column names to correct. If `NULL`, uses
  all measure columns.

- marker_range:

  Numeric vector of length 2 specifying the expected elution range
  `c(min, max)` where the flow marker peak is expected. Required unless
  `auto_detect = TRUE` with `expected_volume` specified.

- target_volume:

  The target elution volume to align the flow marker to. If `NULL`
  (default), uses the first value of `marker_range`.

- auto_detect:

  Logical. If `TRUE`, automatically detect the flow marker peak within
  `marker_range`. If `FALSE`, uses the peak maximum within range.
  Default is `TRUE`.

- min_peak_height:

  Minimum peak height (in signal units) for a valid flow marker
  detection. Peaks below this threshold are ignored. Default is `NULL`
  (no minimum).

- store_correction:

  Logical. If `TRUE`, stores the correction factor as a column named
  `flow_marker_correction`. Default is `TRUE`.

- role:

  Role for generated columns.

- trained:

  Logical indicating if the step has been trained.

- skip:

  Logical. Should the step be skipped when baking?

- id:

  Unique step identifier.

## Value

An updated recipe with the new step added.

## Details

Flow marker correction compensates for small variations in flow rate
between runs. A flow marker is a small molecule (often toluene) that
elutes near the total permeation volume and provides a reference point
for alignment.

**Correction Algorithm:**

1.  Find the flow marker peak maximum within `marker_range`

2.  Calculate the shift: `correction = observed_volume - target_volume`

3.  Apply linear correction:
    `corrected_volume = original_volume - correction`

**Auto-Detection:**

When `auto_detect = TRUE`, the step uses second-derivative analysis to
find the sharpest peak in the specified range, which is typically the
flow marker. This is more robust than simply finding the maximum signal.

**Prerequisites:**

- Should be applied before calibration and MW calculations

- Best applied after baseline correction

## Examples

``` r
if (FALSE) { # \dontrun{
library(recipes)
library(measure)

# Apply flow marker correction with known range
rec <- recipe(~., data = sec_data) |>
  step_measure_input_long(ri, location = vars(elution_volume)) |>
  step_sec_baseline() |>
  step_sec_flow_marker(
    marker_range = c(18, 20),
    target_volume = 18.5
  ) |>
  step_sec_conventional_cal(standards = ps_standards) |>
  prep()

# Auto-detect flow marker with default target (first value of range)
rec <- recipe(~., data = sec_data) |>
  step_measure_input_long(ri, location = vars(elution_volume)) |>
  step_sec_flow_marker(marker_range = c(18, 20)) |>
  prep()
} # }
```
