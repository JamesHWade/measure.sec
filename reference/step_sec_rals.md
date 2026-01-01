# Right-Angle Light Scattering Processing for SEC

`step_sec_rals()` creates a *specification* of a recipe step that
processes right-angle light scattering (RALS) signals for absolute
molecular weight.

## Usage

``` r
step_sec_rals(
  recipe,
  measures = NULL,
  concentration_col = NULL,
  angle = 90,
  laser_wavelength = 670,
  dn_dc = NULL,
  solvent_ri = 1.333,
  optical_constant = NULL,
  calibration_constant = NULL,
  output_mw = "mw_rals",
  min_signal = 0.01,
  role = NA,
  trained = FALSE,
  skip = FALSE,
  id = recipes::rand_id("sec_rals")
)
```

## Arguments

- recipe:

  A recipe object.

- measures:

  Character vector of RALS measure columns. If `NULL`, the step searches
  for measure columns containing "rals".

- concentration_col:

  Name of the concentration measure column (from
  [`step_sec_concentration()`](step_sec_concentration.md) or similar).

- angle:

  Detection angle in degrees. Default is 90.

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

  RALS instrument calibration constant. If `NULL`, results are in
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

RALS provides absolute MW using a 90-degree detector. It is most
accurate for smaller molecules where angular dependence is minimal (Rg
\<\< lambda/20).

## See also

Other sec-detectors:
[`step_sec_concentration()`](step_sec_concentration.md),
[`step_sec_dad()`](step_sec_dad.md),
[`step_sec_dls()`](step_sec_dls.md),
[`step_sec_intrinsic_visc()`](step_sec_intrinsic_visc.md),
[`step_sec_lals()`](step_sec_lals.md),
[`step_sec_mals()`](step_sec_mals.md),
[`step_sec_ri()`](step_sec_ri.md), [`step_sec_uv()`](step_sec_uv.md),
[`step_sec_viscometer()`](step_sec_viscometer.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(recipes)
library(measure)

rec <- recipe(~., data = sec_triple_detect) |>
  step_measure_input_long(ri_signal, location = vars(elution_time), col_name = "ri") |>
  step_measure_input_long(rals_signal, location = vars(elution_time), col_name = "rals") |>
  step_sec_concentration(detector = "ri", injection_mass = 0.2, flow_rate = 1.0) |>
  step_sec_rals(measures = "rals", concentration_col = "ri", dn_dc = 0.185) |>
  prep()
} # }
```
