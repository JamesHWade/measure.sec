# Universal Calibration for SEC

`step_sec_universal_cal()` creates a *specification* of a recipe step
that applies universal calibration to determine molecular weight from
intrinsic viscosity and retention data.

## Usage

``` r
step_sec_universal_cal(
  recipe,
  measures = NULL,
  calibration = NULL,
  calibration_col = NULL,
  intrinsic_visc_col = NULL,
  K_sample,
  a_sample,
  K_standard = 0.000114,
  a_standard = 0.716,
  output_col = "mw_universal",
  role = NA,
  trained = FALSE,
  skip = FALSE,
  id = recipes::rand_id("sec_universal_cal")
)
```

## Arguments

- recipe:

  A recipe object.

- measures:

  Character vector of measure columns to apply calibration to. If
  `NULL`, uses all measure columns.

- calibration:

  A calibration object or data frame containing the universal
  calibration curve (log(\[eta\]M) vs retention).

- calibration_col:

  If calibration is a column name in the data, specify it here.

- intrinsic_visc_col:

  Column containing intrinsic viscosity values. Required for converting
  between polymer types.

- K_sample:

  Mark-Houwink K parameter for the sample polymer.

- a_sample:

  Mark-Houwink a (alpha) exponent for the sample polymer.

- K_standard:

  Mark-Houwink K for the calibration standard polymer. Default is
  0.000114 (polystyrene in THF).

- a_standard:

  Mark-Houwink a for the calibration standard. Default is 0.716
  (polystyrene in THF).

- output_col:

  Name for the output molecular weight column. Default is
  `"mw_universal"`.

- role:

  Role for generated columns.

- trained:

  Logical indicating if the step has been trained.

- skip:

  Logical. Should the step be skipped when baking?

- id:

  Unique step identifier.

## Value

An updated recipe with molecular weight calculated via universal
calibration.

## Details

Universal calibration is based on the principle that polymers with the
same hydrodynamic volume elute at the same retention time, regardless of
chemical structure. The hydrodynamic volume is proportional to \[eta\]M:

\$\$V_h \propto \[\eta\] \cdot M\$\$

**The Universal Calibration Curve:** \$\$\log(\[\eta\] \cdot
M)\_{sample} = \log(\[\eta\] \cdot M)\_{standard}\$\$

At the same retention volume, using Mark-Houwink equations: \$\$\[\eta\]
= K \cdot M^a\$\$

We can solve for sample MW: \$\$M\_{sample} = \left(\frac{K\_{std} \cdot
M\_{std}^{1+a\_{std}}}{K\_{sample}}\right)^{\frac{1}{1+a\_{sample}}}\$\$

**Mark-Houwink Parameters (THF, 25C):**

- Polystyrene: K = 0.000114, a = 0.716

- PMMA: K = 0.000128, a = 0.690

- Polyisoprene: K = 0.000251, a = 0.728

- Polybutadiene: K = 0.000457, a = 0.693

## Note

Universal calibration requires:

- Known Mark-Houwink parameters for both standard and sample

- Calibration with narrow standards (typically polystyrene)

- Same solvent and temperature for all measurements

For absolute MW determination, consider using MALS detection instead.

## See also

Other sec-calibration:
[`sec_calibration_standards`](sec_calibration_standards.md),
[`sec_pmma_standards`](sec_pmma_standards.md),
[`sec_ps_standards`](sec_ps_standards.md),
[`step_sec_conventional_cal()`](step_sec_conventional_cal.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(recipes)
library(measure)

# Apply universal calibration to convert PS calibration to PMMA
rec <- recipe(~., data = pmma_data) |>
  step_measure_input_long(ri, location = vars(time), col_name = "ri") |>
  step_sec_baseline() |>
  step_sec_universal_cal(
    calibration = ps_calibration,
    K_sample = 0.000128,      # PMMA
    a_sample = 0.690,
    K_standard = 0.000114,    # PS (default)
    a_standard = 0.716
  ) |>
  prep()
} # }
```
