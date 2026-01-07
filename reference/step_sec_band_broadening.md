# Band Broadening Correction for SEC

`step_sec_band_broadening()` creates a *specification* of a recipe step
that corrects for axial dispersion (band broadening) in SEC
chromatograms. This improves the accuracy of molecular weight
distribution measurements.

## Usage

``` r
step_sec_band_broadening(
  recipe,
  measures = NULL,
  method = c("tung", "emg"),
  sigma = NULL,
  calibration_peak = NULL,
  tau = NULL,
  iterations = 1,
  damping = 0.5,
  role = NA,
  trained = FALSE,
  skip = FALSE,
  id = recipes::rand_id("sec_band_broadening")
)
```

## Arguments

- recipe:

  A recipe object.

- measures:

  Character vector of measure column names to process. If `NULL`, all
  measure columns will be processed.

- method:

  Correction method. One of:

  - `"tung"` (default): Tung's linear correction

  - `"emg"`: Exponentially Modified Gaussian deconvolution

- sigma:

  Spreading parameter (standard deviation of the instrumental broadening
  function) in the same units as the location axis (typically minutes or
  mL). If `NULL`, must provide `calibration_peak`.

- calibration_peak:

  A `measure_tbl` or data frame with `location` and `value` columns
  representing a narrow standard peak used to estimate sigma.

- tau:

  Exponential time constant for EMG method. If `NULL` with EMG method,
  estimated from `calibration_peak`.

- iterations:

  Number of iterations for iterative correction. Default is 1 for Tung's
  method (single pass). Higher values may improve correction but can
  introduce instability.

- damping:

  Damping factor (0-1) to prevent over-correction and instability.
  Default is 0.5. Lower values are more conservative.

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

Band broadening in SEC occurs due to:

- Axial diffusion during elution

- Non-ideal column packing

- Extra-column volume (tubing, connections, detector cell)

This causes:

- Artificially broadened peaks

- Underestimated Mn (number-average MW)

- Overestimated dispersity (Mw/Mn)

**Tung's Method** (default):

The observed chromatogram F(V) is related to the true distribution W(V)
by: \$\$F(V) = \int W(V') G(V - V') dV'\$\$

where G is a Gaussian spreading function with standard deviation sigma.
Tung's linear correction approximates: \$\$W(V) \approx F(V) - \sigma^2
\frac{d^2 F(V)}{dV^2}\$\$

**EMG Method**:

Models band broadening as convolution with an Exponentially Modified
Gaussian, which better handles asymmetric peak shapes caused by tailing.

**Sigma Determination**:

The spreading parameter sigma should be determined from a narrow
molecular weight standard (e.g., polystyrene with PDI \< 1.05). Use
[`estimate_sigma()`](https://jameshwade.github.io/measure-sec/reference/estimate_sigma.md)
to calculate sigma from such a standard.

## Note

- Correction is applied to the signal, not to molecular weight values

- Large corrections (\> 50% change in peak width) may indicate
  unreliable sigma or poor chromatographic conditions

- This step preserves the area under the curve (mass conservation)

## References

Tung, L. H. (1966). Method of calculating molecular weight distribution
function from gel permeation chromatograms. Journal of Applied Polymer
Science, 10(3), 375-385.

## See also

[`estimate_sigma()`](https://jameshwade.github.io/measure-sec/reference/estimate_sigma.md)
for determining the spreading parameter

Other sec-chromatography:
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

# Using a known sigma value
rec <- recipe(~., data = sec_data) |>
  step_measure_input_long(signal, location = vars(time), col_name = "ri") |>
  step_sec_baseline() |>
  step_sec_band_broadening(sigma = 0.05) |>
  prep()

# Estimating sigma from a narrow standard
narrow_std <- estimate_sigma(narrow_standard_peak)
rec <- recipe(~., data = sec_data) |>
  step_measure_input_long(signal, location = vars(time), col_name = "ri") |>
  step_sec_baseline() |>
  step_sec_band_broadening(sigma = narrow_std$sigma) |>
  prep()
} # }
```
