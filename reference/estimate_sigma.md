# Estimate Spreading Parameter from Narrow Standard

Estimates the instrumental spreading parameter (sigma) from a narrow
molecular weight standard peak. This sigma value can then be used in
[`step_sec_band_broadening()`](https://jameshwade.github.io/measure-sec/reference/step_sec_band_broadening.md)
to correct for band broadening.

## Usage

``` r
estimate_sigma(peak, method = c("gaussian", "fwhm", "moments"))
```

## Arguments

- peak:

  A `measure_tbl` or data frame with `location` and `value` columns
  representing a chromatographic peak from a narrow MW standard.

- method:

  Method for sigma estimation:

  - `"gaussian"` (default): Fit a Gaussian and extract sigma

  - `"fwhm"`: Calculate from full width at half maximum (sigma = FWHM /
    2.355)

  - `"moments"`: Calculate from second moment of the peak

## Value

A list with components:

- sigma:

  The estimated spreading parameter (Gaussian std dev)

- tau:

  Exponential tail parameter (for EMG model), NA if not applicable

- fwhm:

  Full width at half maximum

- asymmetry:

  Peak asymmetry factor (\> 1 indicates tailing)

- method:

  The method used for estimation

## Details

The spreading parameter sigma represents the standard deviation of the
instrumental broadening function, assumed to be approximately Gaussian.

For best results, use a narrow polydispersity standard (PDI \< 1.05) run
under the same conditions as your samples.

**Method Details:**

- **gaussian**: Fits a Gaussian function to the peak using nonlinear
  least squares. Most accurate for symmetric peaks.

- **fwhm**: Uses the relationship sigma = FWHM / 2.355 for a Gaussian.
  Fast and robust but assumes symmetric peak.

- **moments**: Calculates sigma from the second central moment. Accounts
  for asymmetry but sensitive to baseline.

## Examples

``` r
if (FALSE) { # \dontrun{
# From a measure_tbl
sigma_result <- estimate_sigma(narrow_standard_peak)
sigma_result$sigma

# Use in band broadening correction
rec <- recipe(~., data = sec_data) |>
  step_sec_band_broadening(sigma = sigma_result$sigma)
} # }
```
