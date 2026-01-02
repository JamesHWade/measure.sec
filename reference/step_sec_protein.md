# Protein SEC Analysis Workflow

`step_sec_protein()` creates a *specification* of a recipe step that
provides a streamlined workflow for protein SEC analysis, combining
baseline correction, aggregate quantitation, and optionally oligomer
analysis in a single step.

## Usage

``` r
step_sec_protein(
  recipe,
  measures = NULL,
  type = c("native", "denaturing"),
  monomer_mw = NULL,
  monomer_start = NULL,
  monomer_end = NULL,
  extinction_coef = NULL,
  aggregate_threshold = 0.001,
  baseline_method = c("linear", "median", "spline"),
  baseline_left_frac = 0.05,
  baseline_right_frac = 0.05,
  include_oligomer = NULL,
  output_prefix = "protein_",
  role = NA,
  trained = FALSE,
  skip = FALSE,
  id = recipes::rand_id("sec_protein")
)
```

## Arguments

- recipe:

  A recipe object.

- measures:

  Character vector of measure columns to analyze. If `NULL`, analyzes
  all measure columns.

- type:

  Analysis type:

  - `"native"` (default): Native conditions (aqueous buffer,
    non-denaturing)

  - `"denaturing"`: Denaturing conditions (e.g., with SDS or guanidine)

- monomer_mw:

  Expected monomer molecular weight in Da. Required for oligomer
  analysis if `include_oligomer = TRUE`.

- monomer_start:

  Start of the monomer peak region (in location units). If `NULL`,
  automatically determined.

- monomer_end:

  End of the monomer peak region. If `NULL`, automatically determined.

- extinction_coef:

  Extinction coefficient for UV-based concentration. If `NULL`, signal
  remains in raw units.

- aggregate_threshold:

  Minimum fraction of signal to report as aggregate/fragment. Default is
  0.001 (0.1%).

- baseline_method:

  Method for baseline correction. One of `"linear"` (default),
  `"median"`, or `"spline"`.

- baseline_left_frac:

  Fraction of chromatogram start for baseline. Default is 0.05.

- baseline_right_frac:

  Fraction of chromatogram end for baseline. Default is 0.05.

- include_oligomer:

  Logical. Include detailed oligomer analysis? Default is `TRUE` if
  `monomer_mw` is provided.

- output_prefix:

  Prefix for output columns. Default is `"protein_"`.

- role:

  Role for generated columns.

- trained:

  Logical indicating if the step has been trained.

- skip:

  Logical. Should the step be skipped when baking?

- id:

  Unique step identifier.

## Value

An updated recipe with new columns:

- protein_hmws_pct:

  Percent high molecular weight species

- protein_monomer_pct:

  Percent monomer

- protein_lmws_pct:

  Percent low molecular weight species

- protein_main_start:

  Start of main peak region

- protein_main_end:

  End of main peak region

If `include_oligomer = TRUE` and `monomer_mw` is provided:

- protein_monomer_oligo_pct:

  Percent monomer (from oligomer analysis)

- protein_dimer_pct:

  Percent dimer

- protein_trimer_pct:

  Percent trimer

- protein_hmw_oligo_pct:

  Percent HMW oligomers

- protein_lmw_oligo_pct:

  Percent fragments

- protein_species_count:

  Number of detected species

## Details

This step provides a convenient "one-stop" workflow for protein SEC
analysis, suitable for biopharmaceutical characterization. It combines:

1.  **Baseline correction**: SEC-optimized linear or median baseline

2.  **Aggregate quantitation**: HMWS/monomer/LMWS percentages

3.  **Oligomer analysis** (optional): Detailed species identification

**Native vs Denaturing:**

- **Native**: Preserves quaternary structure; use for oligomer analysis

- **Denaturing**: Disrupts non-covalent interactions; use for covalent
  aggregate detection

**Regulatory Context:** Aggregate analysis is critical for
biopharmaceutical characterization:

- ICH Q6B requires aggregate content specification

- USP \<129\> provides guidance on aggregate testing

- Typical specifications: HMWS \< 5%, Monomer \> 95%

**For More Control:** For advanced analysis or custom workflows, use the
individual steps:

- [`step_sec_baseline()`](https://jameshwade.github.io/measure-sec/reference/step_sec_baseline.md)
  for baseline correction

- [`step_sec_uv()`](https://jameshwade.github.io/measure-sec/reference/step_sec_uv.md)
  for UV signal processing

- [`step_sec_aggregates()`](https://jameshwade.github.io/measure-sec/reference/step_sec_aggregates.md)
  for HMWS/monomer/LMWS

- [`step_sec_oligomer()`](https://jameshwade.github.io/measure-sec/reference/step_sec_oligomer.md)
  for detailed species analysis

## See also

[`step_sec_aggregates()`](https://jameshwade.github.io/measure-sec/reference/step_sec_aggregates.md),
[`step_sec_oligomer()`](https://jameshwade.github.io/measure-sec/reference/step_sec_oligomer.md)

Other sec-protein:
[`step_sec_aggregates()`](https://jameshwade.github.io/measure-sec/reference/step_sec_aggregates.md),
[`step_sec_oligomer()`](https://jameshwade.github.io/measure-sec/reference/step_sec_oligomer.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(recipes)
library(measure)

# Basic protein SEC workflow
rec <- recipe(~., data = mab_data) |>
  step_measure_input_long(uv280, location = vars(time), col_name = "uv") |>
  step_sec_protein(monomer_mw = 150000) |>
  prep()

# Native mAb analysis with oligomer detection
rec <- recipe(~., data = mab_data) |>
  step_measure_input_long(uv280, location = vars(time), col_name = "uv") |>
  step_sec_protein(
    type = "native",
    monomer_mw = 150000,
    extinction_coef = 1.4,
    include_oligomer = TRUE
  ) |>
  prep()

# Denaturing conditions (SDS-SEC)
rec <- recipe(~., data = sds_sec_data) |>
  step_measure_input_long(uv280, location = vars(time), col_name = "uv") |>
  step_sec_protein(type = "denaturing", monomer_mw = 150000) |>
  prep()
} # }
```
