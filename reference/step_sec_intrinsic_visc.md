# Calculate Intrinsic Viscosity for SEC

`step_sec_intrinsic_visc()` creates a *specification* of a recipe step
that calculates intrinsic viscosity (\[eta\]) from specific viscosity
and concentration at each elution point.

## Usage

``` r
step_sec_intrinsic_visc(
  recipe,
  specific_visc_col = NULL,
  concentration_col = NULL,
  output_col = "intrinsic_visc",
  min_concentration = 1e-06,
  units = c("dL/g", "mL/g"),
  role = NA,
  trained = FALSE,
  skip = FALSE,
  id = recipes::rand_id("sec_intrinsic_visc")
)
```

## Arguments

- recipe:

  A recipe object.

- specific_visc_col:

  Name of the specific viscosity measure column (from
  [`step_sec_viscometer()`](step_sec_viscometer.md)).

- concentration_col:

  Name of the concentration measure column.

- output_col:

  Name for the intrinsic viscosity output column. Default is
  `"intrinsic_visc"`.

- min_concentration:

  Minimum concentration threshold below which intrinsic viscosity is set
  to NA. Default is 1e-6 mg/mL.

- units:

  Output units for intrinsic viscosity. Default is `"dL/g"`. Common
  alternatives are `"mL/g"` (multiply by 100).

- role:

  Role for generated columns.

- trained:

  Logical indicating if the step has been trained.

- skip:

  Logical. Should the step be skipped when baking?

- id:

  Unique step identifier.

## Value

An updated recipe with the new step added, containing intrinsic
viscosity at each elution point.

## Details

Intrinsic viscosity is defined as the limit of reduced viscosity at
infinite dilution:

\$\$\[\eta\] = \lim\_{c \to 0} \frac{\eta\_{sp}}{c}\$\$

In SEC, each elution slice has very low concentration, so the
approximation \[eta\] = eta_sp / c is valid.

**Applications of Intrinsic Viscosity:**

- Universal Calibration: log(\[eta\] \* M) is linear with retention
  volume, allowing calibration transfer between polymer types

- Mark-Houwink Equation: \[eta\] = K \* M^a, where K and a are polymer-
  and solvent-specific constants

- Branching Analysis: g' = \[eta\]\_branched / \[eta\]\_linear provides
  information about long-chain branching

- Polymer Conformation: The scaling exponent in \[eta\] vs M reveals
  chain conformation (coil, rod, sphere)

**Typical Intrinsic Viscosity Values (dL/g):**

- Polystyrene in THF (MW 100k): ~0.5

- PEG in water (MW 10k): ~0.2

- Proteins: 0.03-0.05 (globular), 0.2-1.0 (denatured)

## See also

[`step_sec_viscometer()`](step_sec_viscometer.md) for specific viscosity
calculation

Other sec-detectors:
[`step_sec_concentration()`](step_sec_concentration.md),
[`step_sec_dad()`](step_sec_dad.md),
[`step_sec_dls()`](step_sec_dls.md),
[`step_sec_lals()`](step_sec_lals.md),
[`step_sec_mals()`](step_sec_mals.md),
[`step_sec_rals()`](step_sec_rals.md),
[`step_sec_ri()`](step_sec_ri.md), [`step_sec_uv()`](step_sec_uv.md),
[`step_sec_viscometer()`](step_sec_viscometer.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(recipes)
library(measure)

# Calculate intrinsic viscosity
rec <- recipe(~., data = sec_visc_data) |>
  step_measure_input_long(dp_signal, location = vars(elution_time), col_name = "dp") |>
  step_measure_input_long(ri_signal, location = vars(elution_time), col_name = "ri") |>
  step_sec_baseline() |>
  step_sec_viscometer(dp_col = "dp") |>
  step_sec_ri(dn_dc = 0.185) |>
  step_sec_concentration(detector = "ri", injection_mass = 0.2, flow_rate = 1.0) |>
  step_sec_intrinsic_visc(
    specific_visc_col = "specific_visc",
    concentration_col = "ri"
  ) |>
  prep()
} # }
```
