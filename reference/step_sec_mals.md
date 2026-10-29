# Multi-Angle Light Scattering Processing for SEC

`step_sec_mals()` creates a *specification* of a recipe step that
processes one light-scattering detector signal using a single-angle
approximation. Multi-angle extrapolation and radius of gyration are not
yet implemented.

## Usage

``` r
step_sec_mals(
  recipe,
  mals_col = NULL,
  concentration_col = NULL,
  dn_dc = NULL,
  dn_dc_column = NULL,
  wavelength = 658,
  solvent_ri = 1.333,
  angles = 90,
  formalism = c("zimm", "debye", "berry"),
  calibration_constant = NULL,
  output_mw = "mw_mals",
  output_rg = NULL,
  min_signal = 0.01,
  role = NA,
  trained = FALSE,
  skip = FALSE,
  id = recipes::rand_id("sec_mals"),
  concentration_unit = NULL
)
```

## Arguments

- recipe:

  A recipe object.

- mals_col:

  Name of the MALS detector measure column. For single-angle detectors
  (RALS/LALS), this is the only signal needed.

- concentration_col:

  Name of the concentration measure column (from
  [`step_sec_concentration()`](https://jameshwade.github.io/measure.sec/reference/step_sec_concentration.md)
  or similar).

- dn_dc:

  Refractive index increment (mL/g). May be signed but must be finite
  and nonzero because the optical constant uses its square. Required for
  absolute MW unless `dn_dc_column` is provided.

- dn_dc_column:

  Column containing sample-specific finite, nonzero dn/dc values; either
  sign is valid.

- wavelength:

  Laser wavelength in nm. Default is 658 (common for MALS).

- solvent_ri:

  Solvent refractive index. Default is 1.333 (water). Common values:
  water = 1.333, THF = 1.407, toluene = 1.497.

- angles:

  One detection angle in degrees, default 90. The calculation assumes
  the particle form factor is approximately 1 at this angle.

- formalism:

  Reserved for future angular extrapolation. Only the default `"zimm"`
  is accepted; no angular fit is performed by this approximation.

- calibration_constant:

  Factor converting the detector signal to excess Rayleigh ratio
  (cm^-1). If `NULL`, results are relative and uncalibrated. The
  appropriate factor depends on the instrument's signal convention.

- output_mw:

  Name for the molecular weight output column. Default is `"mw_mals"`.

- output_rg:

  Must be `NULL`; radius of gyration is not yet implemented.

- min_signal:

  Minimum signal threshold (as fraction of max) below which MW is set to
  NA. Default is 0.01.

- role:

  Role for generated columns.

- trained:

  Logical indicating if the step has been trained.

- skip:

  Logical. Should the step be skipped when baking?

- id:

  Unique step identifier.

- concentration_unit:

  Unit of the concentration curve: `"mg/mL"` or `"g/mL"`. When omitted,
  the training curve must record its unit. If prep infers the unit from
  that tag, baked curves must retain it. An explicit unit permits
  untagged baked curves, but any recorded tag must agree.

## Value

An updated recipe with a single-angle MW estimate at each elution point.
With no calibration constant the output is explicitly uncalibrated.

## Details

Calibrated light scattering can estimate molecular weight without
polymer standards. This approximation also requires an instrument
calibration constant and a concentration signal with known units. The
relationship is:

\$\$\frac{K \cdot c}{R(\theta)} = \frac{1}{M_w} \cdot P(\theta) + 2A_2
c\$\$

where:

- K is the optical constant

- c is the concentration

- R(theta) is the excess Rayleigh ratio

- Mw is the weight-average molecular weight

- P(theta) is the particle scattering function

- A2 is the second virial coefficient

The optical constant K is calculated as: \$\$K = \frac{4\pi^2 n_0^2
(dn/dc)^2}{N_A \lambda^4}\$\$

The implementation uses a single-angle approximation with P(theta) ~ 1;
it does not calculate an angular extrapolation or Rg.

## Note

Accurate results require:

- Known and accurate dn/dc value

- Calibrated instrument (calibration_constant from toluene)

- Accurate concentration from RI detector

- Clean baseline and aligned detectors

## See also

Other sec-detectors:
[`step_sec_concentration()`](https://jameshwade.github.io/measure.sec/reference/step_sec_concentration.md),
[`step_sec_dad()`](https://jameshwade.github.io/measure.sec/reference/step_sec_dad.md),
[`step_sec_dls()`](https://jameshwade.github.io/measure.sec/reference/step_sec_dls.md),
[`step_sec_intrinsic_visc()`](https://jameshwade.github.io/measure.sec/reference/step_sec_intrinsic_visc.md),
[`step_sec_lals()`](https://jameshwade.github.io/measure.sec/reference/step_sec_lals.md),
[`step_sec_rals()`](https://jameshwade.github.io/measure.sec/reference/step_sec_rals.md),
[`step_sec_ri()`](https://jameshwade.github.io/measure.sec/reference/step_sec_ri.md),
[`step_sec_uv()`](https://jameshwade.github.io/measure.sec/reference/step_sec_uv.md),
[`step_sec_viscometer()`](https://jameshwade.github.io/measure.sec/reference/step_sec_viscometer.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(recipes)
library(measure)

# Single-angle (90 degree) light scattering
rec <- recipe(~., data = sec_triple_detect) |>
  step_measure_input_long(ri_signal, location = vars(elution_time), col_name = "ri") |>
  step_measure_input_long(mals_signal, location = vars(elution_time), col_name = "mals") |>
  step_sec_baseline() |>
  step_sec_ri(dn_dc = 0.185) |>
  step_sec_concentration(detector = "ri", injection_mass = 0.2, flow_rate = 1.0) |>
  step_sec_mals(
    mals_col = "mals",
    concentration_col = "ri",
    dn_dc = 0.185,
    wavelength = 658,
    angles = 90
  ) |>
  prep()
} # }
```
