# Oligomeric Species Analysis for Protein SEC

`step_sec_oligomer()` creates a *specification* of a recipe step that
identifies and quantifies individual oligomeric species (monomer, dimer,
trimer, etc.) in protein SEC chromatograms.

## Usage

``` r
step_sec_oligomer(
  recipe,
  measures = NULL,
  monomer_mw = NULL,
  peak_detection = c("auto", "manual"),
  peaks = NULL,
  species = c("monomer", "dimer", "trimer", "hmw", "lmw"),
  mw_tolerance = 0.15,
  mw_column = NULL,
  min_area_pct = 0.1,
  output_prefix = "oligo_",
  role = NA,
  trained = FALSE,
  skip = FALSE,
  id = recipes::rand_id("sec_oligomer")
)
```

## Arguments

- recipe:

  A recipe object.

- measures:

  Character vector of measure columns to analyze. If `NULL`, analyzes
  all measure columns.

- monomer_mw:

  Expected monomer molecular weight in Da. Required for species
  assignment. Typical range: 10,000 - 1,000,000 Da for proteins.

- peak_detection:

  Method for peak detection:

  - `"auto"` (default): Automatic peak detection using derivative
    analysis

  - `"manual"`: Use peaks specified in `peaks` argument

- peaks:

  For manual mode, a list or data frame with peak definitions. Each peak
  should have `start` and `end` retention times.

- species:

  Character vector of species to identify. Default includes `"monomer"`,
  `"dimer"`, `"trimer"`, `"hmw"`, `"lmw"`.

- mw_tolerance:

  Tolerance for MW-based species assignment as a fraction. Default is
  0.15 (15%). A peak is assigned to "dimer" if its MW is within 15% of
  2x the monomer MW.

- mw_column:

  Optional name of a molecular weight measure column (from MALS or
  LALS). If provided, uses MW for species assignment; otherwise uses
  retention time patterns.

- min_area_pct:

  Minimum peak area percentage to report. Peaks below this threshold are
  grouped into "other". Default is 0.1 (0.1%).

- output_prefix:

  Prefix for output columns. Default is `"oligo_"`.

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

- oligo_monomer_pct:

  Percent monomer

- oligo_dimer_pct:

  Percent dimer

- oligo_trimer_pct:

  Percent trimer

- oligo_hmw_pct:

  Percent high molecular weight (\> trimer)

- oligo_lmw_pct:

  Percent low molecular weight (fragments)

- oligo_species_count:

  Number of detected species

- oligo_monomer_mw:

  Observed monomer MW (if mw_column provided)

- oligo_dimer_mw:

  Observed dimer MW (if mw_column provided)

## Details

This step extends [`step_sec_aggregates()`](step_sec_aggregates.md) by
providing detailed species identification rather than just
HMWS/monomer/LMWS classification.

**Species Assignment:**

- Monomer: MW within `mw_tolerance` of `monomer_mw`

- Dimer: MW within `mw_tolerance` of 2x `monomer_mw`

- Trimer: MW within `mw_tolerance` of 3x `monomer_mw`

- HMW: MW \> 3x `monomer_mw`

- LMW: MW \< `monomer_mw` (fragments, clips)

**Peak Detection:**

When using `"auto"` detection, peaks are identified using:

1.  Signal derivative analysis

2.  Local maxima identification

3.  Valley detection for peak boundaries

**Without MW Data:**

If no MW column is available, species are assigned based on retention
time:

- Largest peak is assumed to be monomer

- Earlier-eluting peaks are HMW/oligomers

- Later-eluting peaks are LMW/fragments

## Note

For best results:

- Provide `monomer_mw` for accurate species assignment

- Use with MALS/LALS data for MW-based assignment

- Ensure good chromatographic resolution between species

## See also

[`step_sec_aggregates()`](step_sec_aggregates.md) for simpler
HMWS/monomer/LMWS analysis

Other sec-protein: [`step_sec_aggregates()`](step_sec_aggregates.md),
[`step_sec_protein()`](step_sec_protein.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(recipes)
library(measure)

# Analyze IgG oligomers (monomer MW ~150 kDa)
rec <- recipe(~., data = igg_data) |>
  step_measure_input_long(uv280, location = vars(time), col_name = "uv") |>
  step_sec_baseline() |>
  step_sec_oligomer(monomer_mw = 150000) |>
  prep()

# With MALS-derived MW for accurate assignment
rec <- recipe(~., data = igg_mals_data) |>
  step_measure_input_long(uv280, location = vars(time), col_name = "uv") |>
  step_sec_mals() |>
  step_sec_oligomer(monomer_mw = 150000, mw_column = "mw_mals") |>
  prep()
} # }
```
