# Differential Viscometer Processing for SEC

`step_sec_viscometer()` creates a *specification* of a recipe step that
processes differential viscometer signals to calculate specific
viscosity at each elution point.

## Usage

``` r
step_sec_viscometer(
  recipe,
  dp_col = NULL,
  ip_col = NULL,
  output_col = "specific_visc",
  viscometer_constant = 1,
  min_signal = 0.01,
  role = NA,
  trained = FALSE,
  skip = FALSE,
  id = recipes::rand_id("sec_viscometer")
)
```

## Arguments

- recipe:

  A recipe object.

- dp_col:

  Name of the differential pressure (DP) measure column.

- ip_col:

  Name of the inlet pressure (IP) measure column. If `NULL`, assumes a
  single-capillary viscometer where DP is directly proportional to
  specific viscosity.

- output_col:

  Name for the specific viscosity output column. Default is
  `"specific_visc"`.

- viscometer_constant:

  Instrument calibration constant. Default is 1.0. Obtain from viscosity
  standard calibration.

- min_signal:

  Minimum signal threshold (as fraction of max) below which viscosity is
  set to NA. Default is 0.01.

- role:

  Role for generated columns.

- trained:

  Logical indicating if the step has been trained.

- skip:

  Logical. Should the step be skipped when baking?

- id:

  Unique step identifier.

## Value

An updated recipe with the new step added, containing specific viscosity
at each elution point.

## Details

Differential viscometers measure the pressure difference between a
sample and reference capillary to determine solution viscosity. The
specific viscosity is calculated from the differential pressure (DP) and
inlet pressure (IP):

\$\$\eta\_{sp} = \frac{4 \cdot DP}{IP - 2 \cdot DP}\$\$

For single-capillary viscometers or when IP is not available:
\$\$\eta\_{sp} = K \cdot DP\$\$

where K is a calibration constant.

**Viscometry in SEC:**

- Provides specific viscosity at each MW slice

- Combined with concentration gives intrinsic viscosity \[eta\]

- Used for universal calibration: log(\[eta\] \* M) vs retention

- Essential for branching analysis (g' = \[eta\]\_branched /
  \[eta\]\_linear)

## Note

For intrinsic viscosity calculation, use
[`step_sec_intrinsic_visc()`](step_sec_intrinsic_visc.md) after this
step.

## See also

[`step_sec_intrinsic_visc()`](step_sec_intrinsic_visc.md) for intrinsic
viscosity calculation

Other sec-detectors:
[`step_sec_concentration()`](step_sec_concentration.md),
[`step_sec_dad()`](step_sec_dad.md),
[`step_sec_dls()`](step_sec_dls.md),
[`step_sec_intrinsic_visc()`](step_sec_intrinsic_visc.md),
[`step_sec_lals()`](step_sec_lals.md),
[`step_sec_mals()`](step_sec_mals.md),
[`step_sec_rals()`](step_sec_rals.md),
[`step_sec_ri()`](step_sec_ri.md), [`step_sec_uv()`](step_sec_uv.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(recipes)
library(measure)

# Process differential viscometer data
rec <- recipe(~., data = sec_visc_data) |>
  step_measure_input_long(dp_signal, location = vars(elution_time), col_name = "dp") |>
  step_measure_input_long(ip_signal, location = vars(elution_time), col_name = "ip") |>
  step_sec_baseline() |>
  step_sec_viscometer(dp_col = "dp", ip_col = "ip") |>
  prep()
} # }
```
