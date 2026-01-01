# Multi-Angle Light Scattering Processing for SEC

`step_sec_mals()` creates a *specification* of a recipe step that
processes multi-angle light scattering (MALS) detector signals to
determine absolute molecular weight and radius of gyration at each
elution slice.

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
  output_rg = "rg",
  min_signal = 0.01,
  role = NA,
  trained = FALSE,
  skip = FALSE,
  id = recipes::rand_id("sec_mals")
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
  [`step_sec_concentration()`](step_sec_concentration.md) or similar).

- dn_dc:

  Refractive index increment (mL/g). Required for absolute MW.

- dn_dc_column:

  Column containing sample-specific dn/dc values.

- wavelength:

  Laser wavelength in nm. Default is 658 (common for MALS).

- solvent_ri:

  Solvent refractive index. Default is 1.333 (water). Common values:
  water = 1.333, THF = 1.407, toluene = 1.497.

- angles:

  Numeric vector of detection angles in degrees. For single-angle
  detectors, provide just one value (e.g., `90` for RALS). Default
  assumes a 90-degree detector.

- formalism:

  Angular extrapolation method for multi-angle data: `"zimm"` (default),
  `"debye"`, or `"berry"`.

- calibration_constant:

  MALS instrument calibration constant. If `NULL`, results are in
  relative units. Obtain from toluene standard calibration.

- output_mw:

  Name for the molecular weight output column. Default is `"mw_mals"`.

- output_rg:

  Name for the radius of gyration output column. Default is `"rg"`. Only
  calculated if multiple angles are provided.

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

## Value

An updated recipe with the new step added, containing absolute MW (and
optionally Rg) at each elution point.

## Details

Light scattering provides absolute molecular weight without calibration
standards. The fundamental relationship is:

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

**Formalisms for Angular Extrapolation:**

- Zimm: Kc/R vs sin^2(theta/2) - best for random coils

- Debye: Kc/R vs sin^2(theta/2) - similar to Zimm

- Berry: sqrt(Kc/R) vs sin^2(theta/2) - better for large particles

**Single-Angle vs Multi-Angle:**

- Single angle (RALS/LALS): Provides Mw only, assumes P(theta) ~ 1

- Multi-angle (MALS): Provides both Mw and Rg from angular dependence

## Note

Accurate results require:

- Known and accurate dn/dc value

- Calibrated instrument (calibration_constant from toluene)

- Accurate concentration from RI detector

- Clean baseline and aligned detectors

## See also

Other sec-detectors:
[`step_sec_concentration()`](step_sec_concentration.md),
[`step_sec_dad()`](step_sec_dad.md),
[`step_sec_dls()`](step_sec_dls.md),
[`step_sec_intrinsic_visc()`](step_sec_intrinsic_visc.md),
[`step_sec_lals()`](step_sec_lals.md),
[`step_sec_rals()`](step_sec_rals.md),
[`step_sec_ri()`](step_sec_ri.md), [`step_sec_uv()`](step_sec_uv.md),
[`step_sec_viscometer()`](step_sec_viscometer.md)

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
