# RI Detector Processing for SEC

`step_sec_ri()` creates a *specification* of a recipe step that
processes refractive index (RI) detector signals for SEC analysis,
including application of dn/dc (refractive index increment) values.

## Usage

``` r
step_sec_ri(
  recipe,
  measures = NULL,
  dn_dc = NULL,
  dn_dc_column = NULL,
  instrument_constant = 1,
  output_col = NULL,
  role = NA,
  trained = FALSE,
  skip = FALSE,
  id = recipes::rand_id("sec_ri")
)
```

## Arguments

- recipe:

  A recipe object.

- measures:

  Character vector of RI detector column names to process. If `NULL`,
  will look for columns containing "ri" in the name.

- dn_dc:

  Refractive index increment (mL/g). Can be:

  - A single numeric value applied to all samples

  - `NULL` to skip dn/dc normalization (signal remains in detector
    units)

- dn_dc_column:

  Character name of a column containing sample-specific dn/dc values.
  Overrides `dn_dc` if provided.

- instrument_constant:

  RI detector instrument constant. Default is 1.0 (no adjustment). This
  converts raw detector response to refractive index units.

- output_col:

  Name for the output column. Default is to modify in place.

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

The refractive index (RI) detector is the most common concentration
detector in SEC/GPC. The detector response is proportional to the
product of concentration and dn/dc:

\$\$RI\_{signal} = K \times c \times (dn/dc)\$\$

where:

- K is the instrument constant

- c is the concentration (mg/mL)

- dn/dc is the refractive index increment (mL/g)

**Common dn/dc values (in water at 25°C, 633 nm):**

- Polystyrene in THF: 0.185 mL/g

- PMMA in THF: 0.084 mL/g

- PEG in water: 0.135 mL/g

- Proteins: ~0.185 mL/g

- DNA: ~0.170 mL/g

## See also

Other sec-detectors:
[`step_sec_concentration()`](https://jameshwade.github.io/measure-sec/reference/step_sec_concentration.md),
[`step_sec_dad()`](https://jameshwade.github.io/measure-sec/reference/step_sec_dad.md),
[`step_sec_dls()`](https://jameshwade.github.io/measure-sec/reference/step_sec_dls.md),
[`step_sec_intrinsic_visc()`](https://jameshwade.github.io/measure-sec/reference/step_sec_intrinsic_visc.md),
[`step_sec_lals()`](https://jameshwade.github.io/measure-sec/reference/step_sec_lals.md),
[`step_sec_mals()`](https://jameshwade.github.io/measure-sec/reference/step_sec_mals.md),
[`step_sec_rals()`](https://jameshwade.github.io/measure-sec/reference/step_sec_rals.md),
[`step_sec_uv()`](https://jameshwade.github.io/measure-sec/reference/step_sec_uv.md),
[`step_sec_viscometer()`](https://jameshwade.github.io/measure-sec/reference/step_sec_viscometer.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(recipes)
library(measure)

# Apply fixed dn/dc value
rec <- recipe(~., data = sec_triple_detect) |>
  step_measure_input_long(ri_signal, location = vars(elution_time), col_name = "ri") |>
  step_sec_ri(dn_dc = 0.185) |>
  prep()

# Use sample-specific dn/dc from a column
rec <- recipe(~., data = sec_triple_detect) |>
  step_measure_input_long(ri_signal, location = vars(elution_time), col_name = "ri") |>
  step_sec_ri(dn_dc_column = "dn_dc") |>
  prep()
} # }
```
