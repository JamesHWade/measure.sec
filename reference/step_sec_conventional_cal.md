# Conventional Calibration for SEC Using Narrow Standards

`step_sec_conventional_cal()` creates a *specification* of a recipe step
that fits a calibration curve from narrow molecular weight standards and
applies it to convert elution time/volume to molecular weight.

## Usage

``` r
step_sec_conventional_cal(
  recipe,
  measures = NULL,
  standards = NULL,
  calibration = NULL,
  fit_type = c("cubic", "quadratic", "linear", "fifth", "gam"),
  extrapolation = c("warn", "none", "linear"),
  output_col = "mw",
  log_output = TRUE,
  role = NA,
  trained = FALSE,
  skip = FALSE,
  id = recipes::rand_id("sec_conventional_cal")
)
```

## Arguments

- recipe:

  A recipe object.

- measures:

  Character vector of measure columns to apply calibration to. If
  `NULL`, uses all measure columns.

- standards:

  A data frame containing calibration standards with columns:

  - `location` (or `time`, `volume`, `retention`): Elution position

  - `log_mw` (or `mw`): Molecular weight (will be log-transformed if
    `mw`) Required unless `calibration` is provided.

- calibration:

  A pre-loaded calibration object from
  [`load_sec_calibration()`](https://jameshwade.github.io/measure-sec/reference/load_sec_calibration.md).
  When provided, skips fitting and uses the saved calibration directly.
  Takes precedence over `standards`.

- fit_type:

  Type of fit for the calibration curve:

  - `"cubic"` (default): Third-order polynomial

  - `"quadratic"`: Second-order polynomial

  - `"linear"`: First-order (linear) fit

  - `"fifth"`: Fifth-order polynomial

  - `"gam"`: Generalized Additive Model with cubic splines (requires
    mgcv)

- extrapolation:

  How to handle data outside the calibration range:

  - `"warn"` (default): Extrapolate but warn

  - `"none"`: Return NA for out-of-range values

  - `"linear"`: Use linear extrapolation at boundaries

- output_col:

  Name for the output molecular weight column. Default is `"mw"`.

- log_output:

  Logical. If `TRUE` (default), output column contains log10(MW). If
  `FALSE`, output contains MW in Daltons.

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

This step performs conventional (also called relative) SEC calibration
using narrow dispersity standards of known molecular weight. The
calibration curve relates elution position to log(MW):

\$\$\log\_{10}(M) = f(V_e)\$\$

where f is a polynomial function and V_e is the elution volume or time.

**Calibration Curve Fitting:**

The calibration is fit using orthogonal polynomials for numerical
stability. At least 3 standards are required for cubic fits, 4 for
quadratic, etc.

**Important Considerations:**

- Standards should bracket the MW range of interest

- Calibration is polymer-specific (different polymers have different
  hydrodynamic volumes at the same MW)

- For cross-polymer comparisons, use universal calibration instead
  ([`step_sec_universal_cal`](https://jameshwade.github.io/measure-sec/reference/step_sec_universal_cal.md))

**Fit Quality Metrics:** The
[`tidy()`](https://generics.r-lib.org/reference/tidy.html) method
returns calibration coefficients and R-squared values for assessing fit
quality. R² \> 0.999 is typical for good calibrations.

## See also

Other sec-calibration:
[`load_sec_calibration()`](https://jameshwade.github.io/measure-sec/reference/load_sec_calibration.md),
[`save_sec_calibration()`](https://jameshwade.github.io/measure-sec/reference/save_sec_calibration.md),
[`sec_calibration_standards`](https://jameshwade.github.io/measure-sec/reference/sec_calibration_standards.md),
[`sec_pmma_standards`](https://jameshwade.github.io/measure-sec/reference/sec_pmma_standards.md),
[`sec_ps_standards`](https://jameshwade.github.io/measure-sec/reference/sec_ps_standards.md),
[`step_sec_broad_standard()`](https://jameshwade.github.io/measure-sec/reference/step_sec_broad_standard.md),
[`step_sec_universal_cal()`](https://jameshwade.github.io/measure-sec/reference/step_sec_universal_cal.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(recipes)
library(measure)

# Create calibration standards data
ps_standards <- data.frame(
  retention = c(12.5, 13.2, 14.1, 15.0, 16.2, 17.5),
  log_mw = c(6.0, 5.5, 5.0, 4.5, 4.0, 3.5)
)

# Apply conventional calibration
rec <- recipe(~., data = polymer_data) |>
  step_measure_input_long(ri, location = vars(time), col_name = "ri") |>
  step_sec_baseline() |>
  step_sec_conventional_cal(
    standards = ps_standards,
    fit_type = "cubic"
  ) |>
  prep()

# Check calibration quality
tidy(rec, number = 3)
} # }
```
