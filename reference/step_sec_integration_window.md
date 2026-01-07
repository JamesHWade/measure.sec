# SEC/GPC Integration Window

`step_sec_integration_window()` creates a *specification* of a recipe
step that defines the integration window (start and end x-axis bounds)
for molecular weight calculations. This step adds an
`.integration_window` column that downstream steps like
[`step_sec_mw_averages()`](https://jameshwade.github.io/measure-sec/reference/step_sec_mw_averages.md)
can use.

## Usage

``` r
step_sec_integration_window(
  recipe,
  measures = NULL,
  start = NULL,
  end = NULL,
  auto_detect = TRUE,
  signal_threshold = 0.01,
  extend_beyond_cal = 0.5,
  calibration_range = NULL,
  min_window_width = 1,
  role = NA,
  trained = FALSE,
  skip = FALSE,
  id = recipes::rand_id("sec_integration_window")
)
```

## Arguments

- recipe:

  A recipe object. The step will be added to the sequence of operations
  for this recipe.

- measures:

  An optional character vector of measure column names to process. If
  `NULL` (the default), all measure columns will be processed.

- start:

  Numeric. Start of integration window (x-axis value, e.g., mL). If
  `NULL` and `auto_detect = TRUE`, determined automatically from data.

- end:

  Numeric. End of integration window (x-axis value, e.g., mL). If `NULL`
  and `auto_detect = TRUE`, determined automatically from data.

- auto_detect:

  Logical. If `TRUE` (default), automatically determine window bounds
  from the data when `start` or `end` is `NULL`.

- signal_threshold:

  Numeric between 0 and 1. When auto-detecting, the fraction of maximum
  signal to use as threshold for defining significant signal region.
  Default is `0.01` (1% of max).

- extend_beyond_cal:

  Numeric. Fraction of calibration range to extend beyond the maximum
  calibration point for capturing low MW species. Default is `0.5` (50%
  extension). Only used when calibration info is available via
  `calibration_range`.

- calibration_range:

  Optional numeric vector `c(min, max)` specifying the calibration range
  in x-axis units. When provided, auto-detection respects these bounds
  and applies `extend_beyond_cal`.

- min_window_width:

  Numeric. Minimum window width to ensure valid integration. Default is
  `1.0` (1 mL for SEC).

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
of any existing operations. An `.integration_window` column will be
added containing a tibble with `start` and `end` for each sample.

## Details

The integration window defines the x-axis region used for molecular
weight calculations. For SEC/GPC, this typically corresponds to elution
volume.

**Auto-detection algorithm:**

1.  Find significant signal region (above `signal_threshold` of max)

2.  Extend slightly beyond signal boundaries

3.  If `calibration_range` provided, constrain start to calibration
    minimum

4.  Allow extension up to `extend_beyond_cal` beyond calibration maximum
    to capture low MW species

5.  Ensure minimum window width

**Output format:** The `.integration_window` column contains a list of
tibbles, one per row, each with columns:

- `start`: Start of integration window

- `end`: End of integration window

## Tidying

When you
[`tidy()`](https://recipes.tidymodels.org/reference/tidy.recipe.html)
this step, a tibble with columns `terms`, `start`, `end`, `auto_detect`,
and `id` is returned.

## See also

[`step_sec_mw_averages()`](https://jameshwade.github.io/measure-sec/reference/step_sec_mw_averages.md)
which uses the integration window,
[`step_sec_conventional_cal()`](https://jameshwade.github.io/measure-sec/reference/step_sec_conventional_cal.md)
for calibration.

Other sec-chromatography:
[`step_sec_band_broadening()`](https://jameshwade.github.io/measure-sec/reference/step_sec_band_broadening.md),
[`step_sec_baseline()`](https://jameshwade.github.io/measure-sec/reference/step_sec_baseline.md),
[`step_sec_detector_delay()`](https://jameshwade.github.io/measure-sec/reference/step_sec_detector_delay.md),
[`step_sec_exclude_regions()`](https://jameshwade.github.io/measure-sec/reference/step_sec_exclude_regions.md),
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

# Auto-detect integration window
rec <- recipe(~., data = sec_data) |>
  step_measure_input_long(ri_signal, location = vars(elution_volume)) |>
  step_sec_integration_window() |>
  prep()

# Explicit window bounds
rec <- recipe(~., data = sec_data) |>
  step_measure_input_long(ri_signal, location = vars(elution_volume)) |>
  step_sec_integration_window(start = 8.0, end = 18.0) |>
  prep()

# With calibration constraints
rec <- recipe(~., data = sec_data) |>
  step_measure_input_long(ri_signal, location = vars(elution_volume)) |>
  step_sec_integration_window(
    calibration_range = c(9.0, 16.0),
    extend_beyond_cal = 0.5
  ) |>
  prep()
} # }
```
