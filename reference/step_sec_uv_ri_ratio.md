# Calculate UV/RI Ratio for Composition Analysis

`step_sec_uv_ri_ratio()` creates a *specification* of a recipe step that
calculates the ratio of UV to RI detector signals at each elution point.
This ratio is useful for detecting compositional heterogeneity in
copolymers and conjugates.

## Usage

``` r
step_sec_uv_ri_ratio(
  recipe,
  uv_col = NULL,
  ri_col = NULL,
  output_col = "uv_ri_ratio",
  min_signal = 0.01,
  smooth = TRUE,
  smooth_window = 5,
  role = NA,
  trained = FALSE,
  skip = FALSE,
  id = recipes::rand_id("sec_uv_ri_ratio")
)
```

## Arguments

- recipe:

  A recipe object.

- uv_col:

  Name of the UV detector measure column.

- ri_col:

  Name of the RI detector measure column.

- output_col:

  Name for the output ratio column. Default is `"uv_ri_ratio"`.

- min_signal:

  Minimum signal threshold (as fraction of max) below which ratio is set
  to NA. Default is 0.01 (1%). Prevents noisy ratios in baseline.

- smooth:

  Logical. Apply smoothing to the ratio? Default is `TRUE`.

- smooth_window:

  Window size for smoothing (number of points). Default is 5.

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

The UV/RI ratio provides information about chemical composition across
the molecular weight distribution:

\$\$Ratio(t) = \frac{UV(t)}{RI(t)} = \frac{\varepsilon \cdot
c(t)}{(dn/dc) \cdot c(t) \cdot K}\$\$

Since concentration cancels out, the ratio reflects the relative
detector response factors, which depend on chemical composition.

**Applications:**

- Copolymer composition drift with molecular weight

- Block copolymer characterization

- PEGylation analysis (protein-PEG conjugates)

- Detection of chemical heterogeneity

- End-group analysis with UV labels

**Interpretation:**

- Constant ratio: Uniform composition across MW

- Increasing ratio with MW: More chromophore in higher MW species

- Decreasing ratio with MW: Less chromophore in higher MW species

## Note

Both UV and RI signals should be baseline-corrected and properly aligned
(using
[`step_sec_detector_delay()`](https://jameshwade.github.io/measure-sec/reference/step_sec_detector_delay.md))
before calculating the ratio.

## See also

Other sec-composition:
[`step_sec_composition()`](https://jameshwade.github.io/measure-sec/reference/step_sec_composition.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(recipes)
library(measure)

# Calculate UV/RI ratio for copolymer
rec <- recipe(~., data = sec_copolymer) |>
  step_measure_input_long(ri_signal, location = vars(elution_time), col_name = "ri") |>
  step_measure_input_long(uv_signal, location = vars(elution_time), col_name = "uv") |>
  step_sec_detector_delay(reference = "ri", delay_volumes = c(uv = 0.05)) |>
  step_sec_baseline() |>
  step_sec_uv_ri_ratio(uv_col = "uv", ri_col = "ri") |>
  prep()
} # }
```
