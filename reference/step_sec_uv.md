# UV Detector Processing for SEC

`step_sec_uv()` creates a *specification* of a recipe step that
processes UV detector signals for SEC analysis, including application of
extinction coefficients for concentration determination.

## Usage

``` r
step_sec_uv(
  recipe,
  measures = NULL,
  extinction_coef = NULL,
  extinction_column = NULL,
  wavelength = 280,
  path_length = 1,
  output_col = NULL,
  role = NA,
  trained = FALSE,
  skip = FALSE,
  id = recipes::rand_id("sec_uv")
)
```

## Arguments

- recipe:

  A recipe object.

- measures:

  Character vector of UV detector column names to process. If `NULL`,
  will look for columns containing "uv" in the name.

- extinction_coef:

  Extinction coefficient in mL/(mg*cm) or L/(g*cm). Can be:

  - A single numeric value applied to all samples

  - `NULL` to skip normalization (signal remains in AU)

- extinction_column:

  Character name of a column containing sample-specific extinction
  coefficients. Overrides `extinction_coef` if provided.

- wavelength:

  UV detection wavelength in nm. For documentation only.

- path_length:

  Path length of the flow cell in cm. Default is 1.0.

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

The UV detector measures absorbance according to the Beer-Lambert law:

\$\$A = \varepsilon \times c \times l\$\$

where:

- A is absorbance (AU)

- epsilon is the molar extinction coefficient in mL/(mg\*cm)

- c is the concentration (mg/mL)

- l is the path length (cm)

This step can convert UV absorbance to concentration-proportional
signals by dividing by the extinction coefficient and path length.

**Common UV applications in SEC:**

- Proteins at 280 nm (aromatic amino acids)

- Nucleic acids at 260 nm

- Conjugated polymers

- UV-active end groups or labels

**UV vs RI for concentration:**

- UV is more sensitive for chromophore-containing analytes

- UV response depends on chemical composition (may vary with MW)

- RI is more universal but less sensitive

- For accurate MW, combine both detectors

## See also

Other sec-detectors:
[`step_sec_concentration()`](step_sec_concentration.md),
[`step_sec_dad()`](step_sec_dad.md),
[`step_sec_dls()`](step_sec_dls.md),
[`step_sec_intrinsic_visc()`](step_sec_intrinsic_visc.md),
[`step_sec_lals()`](step_sec_lals.md),
[`step_sec_mals()`](step_sec_mals.md),
[`step_sec_rals()`](step_sec_rals.md),
[`step_sec_ri()`](step_sec_ri.md),
[`step_sec_viscometer()`](step_sec_viscometer.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(recipes)
library(measure)

# Apply fixed extinction coefficient
rec <- recipe(~., data = sec_triple_detect) |>
  step_measure_input_long(uv_signal, location = vars(elution_time), col_name = "uv") |>
  step_sec_uv(extinction_coef = 1.0, wavelength = 280) |>
  prep()

# Use sample-specific extinction coefficients
rec <- recipe(~., data = sec_triple_detect) |>
  step_measure_input_long(uv_signal, location = vars(elution_time), col_name = "uv") |>
  step_sec_uv(extinction_column = "ext_coef") |>
  prep()
} # }
```
