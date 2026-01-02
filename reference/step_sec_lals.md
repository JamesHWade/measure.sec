# Low-Angle Light Scattering Processing for SEC

`step_sec_lals()` creates a *specification* of a recipe step that
processes low-angle light scattering (LALS) signals for absolute
molecular weight.

## Usage

``` r
step_sec_lals(
  recipe,
  measures = NULL,
  concentration_col = NULL,
  angle = 7,
  laser_wavelength = 670,
  dn_dc = NULL,
  solvent_ri = 1.333,
  optical_constant = NULL,
  calibration_constant = NULL,
  output_mw = "mw_lals",
  min_signal = 0.01,
  role = NA,
  trained = FALSE,
  skip = FALSE,
  id = recipes::rand_id("sec_lals")
)
```

## Arguments

- recipe:

  A recipe object.

- measures:

  Character vector of LALS measure columns. If `NULL`, the step searches
  for measure columns containing "lals".

- concentration_col:

  Name of the concentration measure column (from
  [`step_sec_concentration()`](https://jameshwade.github.io/measure-sec/reference/step_sec_concentration.md)
  or similar).

- angle:

  Detection angle in degrees (must be \< 20 for LALS).

- laser_wavelength:

  Laser wavelength in nm.

- dn_dc:

  Refractive index increment (mL/g). Required unless `optical_constant`
  is provided.

- solvent_ri:

  Solvent refractive index. Default is 1.333 (water).

- optical_constant:

  Optional optical constant K; overrides dn/dc.

- calibration_constant:

  LALS instrument calibration constant. If `NULL`, results are in
  relative units.

- output_mw:

  Name for the molecular weight output column.

- min_signal:

  Minimum signal threshold (as fraction of max) below which MW is set to
  NA.

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

LALS provides absolute MW using a low-angle detector (e.g., ~7 degrees).
It assumes P(theta) ~ 1 for small particles and is most accurate when
angular dependence is minimal.

**When to use LALS vs MALS:**

- **LALS**: Preferred for smaller molecules (Rg \< ~10 nm) or when
  multi-angle data is not available. Single-angle measurement is faster
  and simpler but cannot determine Rg.

- **MALS**: Required for large molecules where angular dependence is
  significant. Provides both Mw and Rg from extrapolation to zero angle.

## See also

Other sec-detectors:
[`step_sec_concentration()`](https://jameshwade.github.io/measure-sec/reference/step_sec_concentration.md),
[`step_sec_dad()`](https://jameshwade.github.io/measure-sec/reference/step_sec_dad.md),
[`step_sec_dls()`](https://jameshwade.github.io/measure-sec/reference/step_sec_dls.md),
[`step_sec_intrinsic_visc()`](https://jameshwade.github.io/measure-sec/reference/step_sec_intrinsic_visc.md),
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
  step_measure_input_long(ri_signal, location = vars(elution_time), col_name = "ri") |>
  step_measure_input_long(lals_signal, location = vars(elution_time), col_name = "lals") |>
  step_sec_concentration(detector = "ri", injection_mass = 0.2, flow_rate = 1.0) |>
  step_sec_lals(measures = "lals", concentration_col = "ri", dn_dc = 0.185) |>
  prep()
} # }
```
