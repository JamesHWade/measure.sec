# Diode Array Detector Processing for SEC

`step_sec_dad()` creates a *specification* of a recipe step that
processes diode array detector (DAD) signals across multiple
wavelengths. It can apply wavelength-specific extinction coefficients
and optionally compute ratios to a reference wavelength.

## Usage

``` r
step_sec_dad(
  recipe,
  measures = NULL,
  wavelengths = c(254, 280, 220),
  extinction_coefs = NULL,
  reference_wavelength = NULL,
  output_prefix = "uv",
  path_length = 1,
  role = NA,
  trained = FALSE,
  skip = FALSE,
  id = recipes::rand_id("sec_dad")
)
```

## Arguments

- recipe:

  A recipe object.

- measures:

  Character vector of DAD/UV measure columns. If `NULL`, the step
  searches for measure columns containing "dad" or "uv".

- wavelengths:

  Numeric vector of wavelengths (nm) aligned with `measures`. If `NULL`,
  attempts to parse wavelengths from `measures` names.

- extinction_coefs:

  Extinction coefficients by wavelength. Accepts a named numeric vector
  (names are wavelengths), an unnamed vector aligned with `wavelengths`,
  or a data frame with columns `wavelength` and `extinction_coef`.

- reference_wavelength:

  Optional wavelength used as denominator for ratio calculations.

- output_prefix:

  Prefix used to name output columns (e.g., `uv_254`).

- path_length:

  Path length of the flow cell in cm. Default is 1.0.

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

For each wavelength, the signal can be normalized using the Beer-Lambert
law:

\$\$A = \varepsilon \times c \times l\$\$

where A is absorbance (AU), epsilon is the extinction coefficient, and l
is the path length. When `reference_wavelength` is provided, the step
additionally creates ratio columns for each wavelength vs the reference.

## See also

Other sec-detectors:
[`step_sec_concentration()`](https://jameshwade.github.io/measure-sec/reference/step_sec_concentration.md),
[`step_sec_dls()`](https://jameshwade.github.io/measure-sec/reference/step_sec_dls.md),
[`step_sec_intrinsic_visc()`](https://jameshwade.github.io/measure-sec/reference/step_sec_intrinsic_visc.md),
[`step_sec_lals()`](https://jameshwade.github.io/measure-sec/reference/step_sec_lals.md),
[`step_sec_mals()`](https://jameshwade.github.io/measure-sec/reference/step_sec_mals.md),
[`step_sec_rals()`](https://jameshwade.github.io/measure-sec/reference/step_sec_rals.md),
[`step_sec_ri()`](https://jameshwade.github.io/measure-sec/reference/step_sec_ri.md),
[`step_sec_uv()`](https://jameshwade.github.io/measure-sec/reference/step_sec_uv.md),
[`step_sec_viscometer()`](https://jameshwade.github.io/measure-sec/reference/step_sec_viscometer.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(recipes)
library(measure)

rec <- recipe(~., data = sec_triple_detect) |>
  step_measure_input_long(uv_254, location = vars(elution_time), col_name = "uv_254") |>
  step_measure_input_long(uv_280, location = vars(elution_time), col_name = "uv_280") |>
  step_sec_dad(
    measures = c("uv_254", "uv_280"),
    wavelengths = c(254, 280),
    extinction_coefs = c(`254` = 1.2, `280` = 1.0),
    reference_wavelength = 280
  ) |>
  prep()
} # }
```
