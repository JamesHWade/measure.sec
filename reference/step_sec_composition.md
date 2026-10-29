# Calculate Copolymer Composition from Detector Signals

`step_sec_composition()` creates a *specification* of a recipe step that
calculates the weight fraction of components in a copolymer or blend
using UV and RI detector signals with known response factors.

## Usage

``` r
step_sec_composition(
  recipe,
  uv_col = NULL,
  ri_col = NULL,
  component_a_uv,
  component_a_ri,
  component_b_uv,
  component_b_ri,
  output_col = "composition_a",
  min_signal = 0.01,
  clip = TRUE,
  role = NA,
  trained = FALSE,
  skip = FALSE,
  id = recipes::rand_id("sec_composition"),
  input_basis = c("detector_response", "normalized_concentration"),
  uv_concentration_unit = NULL,
  ri_concentration_unit = NULL
)
```

## Arguments

- recipe:

  A recipe object.

- uv_col:

  Name of the UV detector measure column.

- ri_col:

  Name of the RI detector measure column.

- component_a_uv, component_b_uv:

  UV response coefficients for components A and B, calibrated in the
  selected input-channel basis per g/mL.

- component_a_ri, component_b_ri:

  RI response coefficients for components A and B, calibrated in the
  selected input-channel basis per g/mL.

- output_col:

  Name for the output composition column. Default is `"composition_a"`.
  Contains weight fraction of component A (0-1).

- min_signal:

  Minimum signal threshold (as fraction of max) below which composition
  is set to `NA`. Default is 0.01. The threshold is a fraction of the
  maximum absolute RI value; its effective units are native RI-response
  units for raw mode and g/mL after conversion for normalized mode.

- clip:

  Logical. Clip composition values to \[0, 1\] range? Default is `TRUE`.

- role:

  Role for generated columns.

- trained:

  Logical indicating if the step has been trained.

- skip:

  Logical. Should the step be skipped when baking?

- id:

  Unique step identifier.

- input_basis:

  Input interpretation: `"detector_response"` (default) uses untagged
  raw UV and RI detector responses; `"normalized_concentration"` uses
  apparent concentration channels and converts them to g/mL before
  calculating composition.

- uv_concentration_unit, ri_concentration_unit:

  Units of the respective curves when
  `input_basis = "normalized_concentration"`. Each can be omitted when
  the training curve records `sec_concentration_unit`; an explicit unit
  must agree with any recorded unit. If prep infers a unit from its
  training tag, the corresponding baked curve must retain that tag. An
  explicit unit permits untagged baked curves, but any recorded tag must
  agree. Both units are required for untagged training curves. These
  arguments must be `NULL` for raw detector responses.

## Value

An updated recipe with the new step added.

## Details

For a two-component system, the detector signals are:

\$\$UV = u_A c_A + u_B c_B, \qquad RI = r_A c_A + r_B c_B\$\$

where `c_A` and `c_B` are the component mass concentrations in g/mL. The
coefficients `u_A`, `u_B`, `r_A`, and `r_B` are the response
coefficients for the selected `input_basis`, measured in the same
input-channel basis per g/mL. In raw detector-response mode these are
calibrated detector responses per g/mL; in normalized-concentration mode
they are calibrated apparent channel responses per g/mL. They are not
inferred from physical extinction coefficients or `dn/dc` values.

Composition is calculated by solving the two-response linear system and
converting the recovered component concentrations to a mass fraction:

\$\$c_A = \frac{UV r_B - u_B RI}{u_A r_B - u_B r_A}, \qquad c_B =
\frac{u_A RI - UV r_A}{u_A r_B - u_B r_A}, \qquad w_A = \frac{c_A}{c_A +
c_B}\$\$

The response matrix must be identifiable: the component response vectors
cannot be proportional. A zero coefficient is allowed when the full
matrix remains identifiable.

`min_signal` is a fraction of the maximum absolute RI value. The
effective cutoff is therefore in raw RI detector units for
`"detector_response"`, or in g/mL after unit conversion for
`"normalized_concentration"`.

**Common Applications:**

- Styrene-acrylate copolymers (styrene is UV-active)

- Block copolymers with different chromophore content

- PEGylated proteins (protein at 280nm, PEG is UV-transparent)

- Polymer blends with known compositions

Calibrate each response coefficient from the measured channel output of
a pure-component standard at known concentrations in g/mL. A literature
extinction coefficient or `dn/dc` value alone is not a substitute;
account for path length, instrument calibration, and signal scaling when
relating it to the actual input-channel values.

## Note

The response coefficients must be calibrated for the selected input
basis and the same g/mL concentration denominator. Raw detector
responses and normalized concentration channels are separate models;
switching between them requires recalibrating all four coefficients. If
both channels have been normalized to the same mass concentration
response for both components, the response matrix is singular and
composition is not identifiable from those channels.

## See also

Other sec-composition:
[`step_sec_uv_ri_ratio()`](https://jameshwade.github.io/measure.sec/reference/step_sec_uv_ri_ratio.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(recipes)
library(measure)

# Calibrate the response coefficients with pure-component standards.
rec <- recipe(~., data = copolymer_data) |>
  step_measure_input_long(ri_signal, location = vars(elution_time), col_name = "ri") |>
  step_measure_input_long(uv_signal, location = vars(elution_time), col_name = "uv") |>
  step_sec_baseline() |>
  step_sec_composition(
    uv_col = "uv",
    ri_col = "ri",
    component_a_uv = 2.0,    # measured UV response per g/mL
    component_a_ri = 1.0,    # measured RI response per g/mL
    component_b_uv = 0.5,
    component_b_ri = 3.0,
    output_col = "styrene_fraction"
  ) |>
  prep()
} # }
```
