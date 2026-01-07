# Broad Standard Calibration for SEC/GPC

`step_sec_broad_standard()` creates a *specification* of a recipe step
that fits a calibration curve using a polydisperse (broad) molecular
weight standard with known Mn and Mw values.

## Usage

``` r
step_sec_broad_standard(
  recipe,
  measures = NULL,
  broad_standard = NULL,
  known_mn = NULL,
  known_mw = NULL,
  fit_type = c("linear", "quadratic"),
  method = c("hamielec", "integral"),
  reference_mwd = NULL,
  integration_range = NULL,
  extrapolation = c("warn", "none"),
  output_col = "mw",
  log_output = TRUE,
  role = NA,
  trained = FALSE,
  skip = FALSE,
  id = recipes::rand_id("sec_broad_standard")
)
```

## Arguments

- recipe:

  A recipe object.

- measures:

  Character vector of measure columns to apply calibration to. If
  `NULL`, uses all measure columns.

- broad_standard:

  A data frame containing the broad standard chromatogram:

  - `location` (or `time`, `volume`, `retention`, `elution_time`,
    `elution_volume`): Elution position

  - `value` (or `signal`, `response`, `intensity`, `ri`, `uv`): Detector
    response

- known_mn:

  Known number-average molecular weight (Mn) in Daltons.

- known_mw:

  Known weight-average molecular weight (Mw) in Daltons.

- fit_type:

  Type of calibration curve:

  - `"linear"` (default): log10(M) = C1 + C2\*V (classic Hamielec)

  - `"quadratic"`: log10(M) = C1 + C2*V + C3*V^2

- method:

  Calibration method:

  - `"hamielec"` (default): Optimize to match Mn and Mw

  - `"integral"`: Use cumulative MWD matching (requires `reference_mwd`)

- reference_mwd:

  Optional data frame for integral method containing the known
  cumulative molecular weight distribution of the broad standard:

  - `mw`: Molecular weight values (Daltons)

  - `cumulative`: Cumulative weight fraction (0 to 1)

  Required when `method = "integral"`. Can be obtained from the
  standard's certificate of analysis or determined by light
  scattering/viscometry.

- integration_range:

  Optional numeric vector `c(min, max)` specifying the elution range to
  use for the broad standard. If `NULL`, auto-detects peak region.

- extrapolation:

  How to handle data outside the calibration range:

  - `"warn"` (default): Extrapolate but warn

  - `"none"`: Return NA for out-of-range values

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

This step implements broad standard calibration methods for SEC/GPC,
which use a single polydisperse standard with known molecular weight
averages instead of multiple narrow standards.

**Hamielec Method:**

Based on Balke, Hamielec, LeClair, and Pearce (1969). The original
method assumes a linear calibration relationship, though this
implementation also supports a quadratic extension for better fit over
wide MW ranges:

\$\$\log\_{10}(M) = C_1 + C_2 \cdot V\$\$

The algorithm simultaneously optimizes coefficients C1 and C2 (and C3
for quadratic fits) using Nelder-Mead optimization to minimize the
squared relative errors between calculated and known Mn and Mw values.

**Integral/Cumulative Match Method:**

The integral method matches the entire cumulative MWD shape rather than
just Mn and Mw. This requires a reference cumulative distribution
(`reference_mwd`) typically obtained from the standard's certificate or
measured by light scattering/viscometry. The algorithm optimizes
calibration coefficients to minimize the sum of squared differences
between the calculated and reference cumulative distributions.

This method is more robust than Hamielec because it uses the full
distribution shape, not just two moments. It's particularly useful when
the standard's MWD shape is well-characterized.

**When to Use Broad Standard Calibration:**

- QC labs running the same polymer type repeatedly

- When narrow standards aren't available for your polymer

- When you have a well-characterized in-house reference material

- Provides "absolute" molecular weights for the same polymer type

**Limitations:**

- Results are only valid for polymers with similar hydrodynamic behavior

- Linear calibration may not fit well over very wide MW ranges

- Requires well-characterized broad standard (accurate Mn and Mw)

- Integral method requires full cumulative MWD data

## References

Balke, S.T., Hamielec, A.E., LeClair, B.P., and Pearce, S.L. (1969). Gel
permeation chromatography. *Industrial & Engineering Chemistry Product
Research and Development*, 8(1), 54-57.

## See also

Other sec-calibration:
[`load_sec_calibration()`](https://jameshwade.github.io/measure-sec/reference/load_sec_calibration.md),
[`save_sec_calibration()`](https://jameshwade.github.io/measure-sec/reference/save_sec_calibration.md),
[`sec_calibration_standards`](https://jameshwade.github.io/measure-sec/reference/sec_calibration_standards.md),
[`sec_pmma_standards`](https://jameshwade.github.io/measure-sec/reference/sec_pmma_standards.md),
[`sec_ps_standards`](https://jameshwade.github.io/measure-sec/reference/sec_ps_standards.md),
[`step_sec_conventional_cal()`](https://jameshwade.github.io/measure-sec/reference/step_sec_conventional_cal.md),
[`step_sec_universal_cal()`](https://jameshwade.github.io/measure-sec/reference/step_sec_universal_cal.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(recipes)
library(measure)

data(sec_triple_detect)

# Create broad standard chromatogram data
broad_std <- data.frame(
  time = seq(10, 20, by = 0.1),
  signal = dnorm(seq(10, 20, by = 0.1), mean = 15, sd = 1.5)
)

# Apply broad standard calibration
rec <- recipe(~., data = sec_triple_detect) |>
  step_measure_input_long(
    ri_signal,
    location = vars(elution_time),
    col_name = "ri"
  ) |>
  step_sec_baseline() |>
  step_sec_broad_standard(
    broad_standard = broad_std,
    known_mn = 50000,
    known_mw = 150000
  ) |>
  prep()

# Check calibration results
tidy(rec, number = 3)
} # }
```
