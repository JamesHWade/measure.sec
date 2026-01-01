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
  id = recipes::rand_id("sec_composition")
)
```

## Arguments

- recipe:

  A recipe object.

- uv_col:

  Name of the UV detector measure column.

- ri_col:

  Name of the RI detector measure column.

- component_a_uv:

  UV response factor for component A (extinction coefficient in
  mL/(mg\*cm) or relative units).

- component_a_ri:

  RI response factor for component A (dn/dc in mL/g or relative units).

- component_b_uv:

  UV response factor for component B.

- component_b_ri:

  RI response factor for component B.

- output_col:

  Name for the output composition column. Default is `"composition_a"`.
  Contains weight fraction of component A (0-1).

- min_signal:

  Minimum signal threshold (as fraction of max) below which composition
  is set to NA. Default is 0.01.

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

## Value

An updated recipe with the new step added.

## Details

For a two-component system, the detector signals are:

\$\$UV = \varepsilon_A \cdot c_A + \varepsilon_B \cdot c_B\$\$ \$\$RI =
(dn/dc)\_A \cdot c_A + (dn/dc)\_B \cdot c_B\$\$

where c_A and c_B are the concentrations of components A and B.

The weight fraction of component A is calculated by solving this system:

\$\$w_A = \frac{R\_{obs} - R_B}{R_A - R_B}\$\$

where R_obs is the observed UV/RI ratio, and R_A, R_B are the ratios for
pure components (e/dn/dc).

**Common Applications:**

- Styrene-acrylate copolymers (styrene is UV-active)

- Block copolymers with different chromophore content

- PEGylated proteins (protein at 280nm, PEG is UV-transparent)

- Polymer blends with known compositions

**Example Response Factors:**

- Polystyrene: UV (254nm) ~ 1.0, dn/dc ~ 0.185

- PMMA: UV (254nm) ~ 0.01, dn/dc ~ 0.084

- PEG: UV (280nm) ~ 0, dn/dc ~ 0.135

- Proteins: UV (280nm) ~ 1.0, dn/dc ~ 0.185

## Note

Response factors must be in consistent units. The absolute values don't
matter as long as the ratios are correct for pure components.

## See also

Other sec-composition:
[`step_sec_uv_ri_ratio()`](step_sec_uv_ri_ratio.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(recipes)
library(measure)

# Styrene-MMA copolymer composition
rec <- recipe(~., data = copolymer_data) |>
  step_measure_input_long(ri_signal, location = vars(elution_time), col_name = "ri") |>
  step_measure_input_long(uv_signal, location = vars(elution_time), col_name = "uv") |>
  step_sec_baseline() |>
  step_sec_composition(
    uv_col = "uv",
    ri_col = "ri",
    component_a_uv = 1.0,    # Styrene (UV-active)
    component_a_ri = 0.185,  # Styrene dn/dc
    component_b_uv = 0.01,   # MMA (weak UV)
    component_b_ri = 0.084,  # MMA dn/dc
    output_col = "styrene_fraction"
  ) |>
  prep()
} # }
```
