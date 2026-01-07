# SEC/GPC Peak Deconvolution

`step_sec_peaks_deconvolve()` creates a *specification* of a recipe step
that resolves overlapping peaks in SEC/GPC chromatograms using curve
fitting. This is a thin wrapper around
[`measure::step_measure_peaks_deconvolve()`](https://jameshwade.github.io/measure/reference/step_measure_peaks_deconvolve.html)
with SEC-optimized defaults.

## Usage

``` r
step_sec_peaks_deconvolve(
  recipe,
  model = "emg",
  optimizer = "auto",
  max_iter = 500L,
  quality_threshold = 0.9,
  smart_init = TRUE,
  constrain_positions = TRUE,
  peaks_col = ".peaks",
  measures_col = ".measures",
  role = NA,
  trained = FALSE,
  skip = FALSE,
  id = recipes::rand_id("sec_peaks_deconvolve")
)
```

## Arguments

- recipe:

  A recipe object. The step will be added to the sequence of operations
  for this recipe.

- model:

  Peak model to use:

  - `"emg"` (default): Exponentially Modified Gaussian, recommended for
    SEC chromatograms with peak tailing

  - `"gaussian"`: Symmetric Gaussian, for well-behaved symmetric peaks

  - `"bigaussian"`: Bi-Gaussian for flexible asymmetry

  EMG is the default because chromatographic peaks typically exhibit
  tailing due to mass transfer kinetics and column effects.

- optimizer:

  Optimization method:

  - `"auto"` (default): Selects based on problem complexity

  - `"lbfgsb"`: L-BFGS-B (fast, local optimization)

  - `"multistart"`: Multiple starts for robustness

  - `"nelder_mead"`: Derivative-free simplex method

- max_iter:

  Maximum iterations for optimization. Default is `500`.

- quality_threshold:

  Minimum R-squared to accept fit. Default is `0.9` (stricter than base
  measure default of 0.8 for analytical quality).

- smart_init:

  Logical. Use smart initialization based on peak properties. Default is
  `TRUE`. Highly recommended for SEC data.

- constrain_positions:

  Logical. Enforce that peak centers maintain their relative ordering
  during optimization. Default is `TRUE`.

- peaks_col:

  Name of the peaks column. Default is `".peaks"`.

- measures_col:

  Name of the measures column containing the chromatogram. Default is
  `".measures"` but often `"ri"`, `"uv"`, etc. in SEC workflows.

- role:

  Not used by this step.

- trained:

  A logical to indicate if the step has been trained.

- skip:

  A logical. Should the step be skipped when baking?

- id:

  A character string that is unique to this step.

## Value

An updated version of `recipe` with the new step added. The `.peaks`
column will be updated with:

- Refined peak parameters from curve fitting

- `fit_r_squared`: R-squared of the overall fit

- `area`: Integrated area under the fitted curve (analytical
  integration)

## Details

Peak deconvolution is essential for SEC/GPC analysis when peaks overlap,
which is common for:

- Bimodal molecular weight distributions

- Aggregate/monomer/fragment species in protein SEC

- Copolymer component separation

**EMG Model (Default):**

The Exponentially Modified Gaussian (EMG) is ideal for SEC because it
accounts for the characteristic tailing observed in chromatographic
peaks:

\$\$EMG(x) = h \cdot \exp(\sigma^2/(2\tau^2) + (c-x)/\tau) \cdot
\Phi((x-c)/\sigma - \sigma/\tau)\$\$

where \\h\\ is height, \\c\\ is center, \\\sigma\\ is Gaussian width,
and \\\tau\\ is the exponential decay parameter (tailing).

**Quality Assessment:**

The `quality_threshold` parameter sets the minimum acceptable R-squared
for fits. The default of 0.9 is stricter than the base measure default,
appropriate for quantitative analytical work.

## Tidying

When you
[`tidy()`](https://recipes.tidymodels.org/reference/tidy.recipe.html)
this step, a tibble with columns `model`, `optimizer`,
`quality_threshold`, and `id` is returned.

## See also

[`step_sec_peaks_detect()`](https://jameshwade.github.io/measure-sec/reference/step_sec_peaks_detect.md)
for SEC-optimized peak detection,
[`measure::step_measure_peaks_deconvolve()`](https://jameshwade.github.io/measure/reference/step_measure_peaks_deconvolve.html)
for the underlying implementation.

Other sec-chromatography:
[`step_sec_band_broadening()`](https://jameshwade.github.io/measure-sec/reference/step_sec_band_broadening.md),
[`step_sec_baseline()`](https://jameshwade.github.io/measure-sec/reference/step_sec_baseline.md),
[`step_sec_detector_delay()`](https://jameshwade.github.io/measure-sec/reference/step_sec_detector_delay.md),
[`step_sec_exclude_regions()`](https://jameshwade.github.io/measure-sec/reference/step_sec_exclude_regions.md),
[`step_sec_integration_window()`](https://jameshwade.github.io/measure-sec/reference/step_sec_integration_window.md),
[`step_sec_mw_averages()`](https://jameshwade.github.io/measure-sec/reference/step_sec_mw_averages.md),
[`step_sec_mw_distribution()`](https://jameshwade.github.io/measure-sec/reference/step_sec_mw_distribution.md),
[`step_sec_mw_fractions()`](https://jameshwade.github.io/measure-sec/reference/step_sec_mw_fractions.md),
[`step_sec_peaks_detect()`](https://jameshwade.github.io/measure-sec/reference/step_sec_peaks_detect.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(recipes)
library(measure)

# SEC peak deconvolution with EMG model (default)
rec <- recipe(~., data = sec_data) |>
  step_measure_input_long(ri_signal, location = vars(elution_volume)) |>
  step_sec_peaks_detect() |>
  step_sec_peaks_deconvolve() |>
  prep()

# Use Gaussian model for symmetric peaks
rec <- recipe(~., data = sec_data) |>
  step_measure_input_long(ri_signal, location = vars(elution_volume)) |>
  step_sec_peaks_detect() |>
  step_sec_peaks_deconvolve(model = "gaussian") |>
  prep()

# Use multistart optimizer for complex multi-peak scenarios
rec <- recipe(~., data = sec_data) |>
  step_measure_input_long(ri_signal, location = vars(elution_volume)) |>
  step_sec_peaks_detect() |>
  step_sec_peaks_deconvolve(optimizer = "multistart") |>
  prep()
} # }
```
