# Quantify Protein Aggregates and Fragments in SEC

`step_sec_aggregates()` creates a *specification* of a recipe step that
quantifies high molecular weight species (HMWS/aggregates), monomers,
and low molecular weight species (LMWS/fragments) from protein SEC
chromatograms.

## Usage

``` r
step_sec_aggregates(
  recipe,
  measures = NULL,
  monomer_start = NULL,
  monomer_end = NULL,
  method = c("tallest", "manual"),
  hmws_threshold = 0.001,
  include_main_peak = TRUE,
  output_prefix = "purity_",
  role = NA,
  trained = FALSE,
  skip = FALSE,
  id = recipes::rand_id("sec_aggregates")
)
```

## Arguments

- recipe:

  A recipe object.

- measures:

  Character vector of measure columns to analyze. If `NULL`, analyzes
  all measure columns.

- monomer_start:

  Start of the monomer peak region (in location units, typically
  minutes). If `NULL`, automatically determined.

- monomer_end:

  End of the monomer peak region. If `NULL`, automatically determined.

- method:

  Method for peak boundary detection when `monomer_start` or
  `monomer_end` is `NULL`:

  - `"tallest"` (default): Uses the tallest peak as monomer

  - `"manual"`: Requires explicit boundaries

- hmws_threshold:

  Minimum fraction of monomer height to consider as HMWS signal. Default
  is 0.001 (0.1%). Below this, signal is considered baseline.

- include_main_peak:

  Logical. Include the main peak boundaries in output? Default is TRUE.

- output_prefix:

  Prefix for output columns. Default is `"purity_"`. Creates columns:
  `{prefix}hmws`, `{prefix}monomer`, `{prefix}lmws`.

- role:

  Role for generated columns.

- trained:

  Logical indicating if the step has been trained.

- skip:

  Logical. Should the step be skipped when baking?

- id:

  Unique step identifier.

## Value

An updated recipe with new columns containing aggregate percentages:

- purity_hmws:

  Percent high molecular weight species (aggregates)

- purity_monomer:

  Percent monomer (main peak)

- purity_lmws:

  Percent low molecular weight species (fragments)

- purity_main_start:

  Start of main peak region (if include_main_peak)

- purity_main_end:

  End of main peak region (if include_main_peak)

## Details

Aggregate analysis is critical for biopharmaceutical characterization:

- HMWS (High Molecular Weight Species): Elute before the monomer peak.
  Includes dimers, trimers, and higher-order aggregates.

- Monomer: The main therapeutic protein peak.

- LMWS (Low Molecular Weight Species): Elute after the monomer peak.
  Includes fragments, clips, and degradation products.

**Calculation Method:** \$\$\\ HMWS = \frac{A\_{HMWS}}{A\_{total}}
\times 100\$\$ \$\$\\ Monomer = \frac{A\_{monomer}}{A\_{total}} \times
100\$\$ \$\$\\ LMWS = \frac{A\_{LMWS}}{A\_{total}} \times 100\$\$

where A represents integrated peak areas.

**Regulatory Importance:**

- ICH Q6B requires aggregate content specification

- USP \<129\> provides guidance on aggregate testing

- Typical acceptance: HMWS \< 5%, Monomer \> 95%

## Note

For accurate results:

- Baseline correct the chromatogram first

- Ensure proper column resolution (especially for dimer separation)

- Use UV detection at 280 nm (or 214 nm for higher sensitivity)

## See also

Other sec-protein: [`step_sec_oligomer()`](step_sec_oligomer.md),
[`step_sec_protein()`](step_sec_protein.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(recipes)
library(measure)

# Analyze mAb aggregate content
rec <- recipe(~., data = mab_data) |>
  step_measure_input_long(uv280, location = vars(elution_time), col_name = "uv") |>
  step_sec_baseline() |>
  step_sec_aggregates(
    monomer_start = 8.5,
    monomer_end = 10.5
  ) |>
  prep()

# Automatic peak detection
rec <- recipe(~., data = mab_data) |>
  step_measure_input_long(uv280, location = vars(elution_time), col_name = "uv") |>
  step_sec_baseline() |>
  step_sec_aggregates(method = "tallest") |>
  prep()
} # }
```
