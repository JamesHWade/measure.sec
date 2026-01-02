# Dynamic Light Scattering Processing for SEC

`step_sec_dls()` creates a *specification* of a recipe step that
processes dynamic light scattering (DLS) correlation data to estimate
diffusion coefficient and hydrodynamic radius (Rh).

## Usage

``` r
step_sec_dls(
  recipe,
  measures = NULL,
  temperature = 25,
  viscosity = NULL,
  laser_wavelength = 633,
  angle = 90,
  solvent_ri = 1.333,
  rg_col = "rg",
  output_cols = c("rh", "diffusion_coef"),
  role = NA,
  trained = FALSE,
  skip = FALSE,
  id = recipes::rand_id("sec_dls")
)
```

## Arguments

- recipe:

  A recipe object.

- measures:

  Character vector of DLS measure columns. If `NULL`, the step searches
  for measure columns containing "dls" or "corr".

- temperature:

  Temperature in degrees Celsius.

- viscosity:

  Solvent viscosity in mPa\*s. If `NULL`, uses a water approximation
  based on temperature.

- laser_wavelength:

  Laser wavelength in nm. Default is 633.

- angle:

  Scattering angle in degrees. Default is 90.

- solvent_ri:

  Solvent refractive index. Default is 1.333 (water).

- rg_col:

  Optional Rg measure column from MALS for Rg/Rh ratio.

- output_cols:

  Output column names for Rh and diffusion coefficient.

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

The DLS correlation function g2(tau) is approximated as: \$\$g2(\tau) =
1 + \beta e^{-2 \Gamma \tau}\$\$ where Gamma is the decay rate. The
diffusion coefficient is calculated as D = Gamma / q^2 with q determined
from the scattering angle and wavelength. Rh is then obtained from the
Stokes-Einstein equation.

## See also

Other sec-detectors:
[`step_sec_concentration()`](https://jameshwade.github.io/measure-sec/reference/step_sec_concentration.md),
[`step_sec_dad()`](https://jameshwade.github.io/measure-sec/reference/step_sec_dad.md),
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
  step_measure_input_long(dls_corr, location = vars(lag_time), col_name = "dls") |>
  step_sec_dls(measures = "dls", temperature = 25, viscosity = 0.89) |>
  prep()
} # }
```
