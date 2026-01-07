# SEC/GPC Region Exclusion

`step_sec_exclude_regions()` creates a *specification* of a recipe step
that marks regions for exclusion from analysis. Excluded regions can be
solvent peaks, artifacts, flow markers, or other features that should
not be included in baseline fitting or molecular weight integration.

## Usage

``` r
step_sec_exclude_regions(
  recipe,
  measures = NULL,
  regions = NULL,
  purpose = c("both", "baseline", "integration"),
  role = NA,
  trained = FALSE,
  skip = FALSE,
  id = recipes::rand_id("sec_exclude_regions")
)
```

## Arguments

- recipe:

  A recipe object. The step will be added to the sequence of operations
  for this recipe.

- measures:

  An optional character vector of measure column names to process. If
  `NULL` (the default), all measure columns will be processed.

- regions:

  A data frame or tibble specifying regions to exclude. Must contain
  columns `start` and `end` (x-axis values defining each region).
  Optional columns: `reason` (character describing why excluded),
  `sample_id` (for sample-specific exclusions).

- purpose:

  Character. What the exclusions apply to:

  - `"baseline"`: Exclude from baseline fitting only

  - `"integration"`: Exclude from MW integration only

  - `"both"` (default): Exclude from both baseline fitting and
    integration

- role:

  Not used by this step since no new variables are created.

- trained:

  A logical to indicate if the quantities for preprocessing have been
  estimated.

- skip:

  A logical. Should the step be skipped when the recipe is baked?

- id:

  A character string that is unique to this step to identify it.

## Value

An updated version of `recipe` with the new step added to the sequence
of any existing operations. An `.excluded_regions` column will be added
containing exclusion information for each sample.

## Details

Excluded regions are stored in a list column `.excluded_regions` where
each element is a tibble with columns:

- `start`: Start of excluded region (x-axis value)

- `end`: End of excluded region (x-axis value)

- `purpose`: What this exclusion applies to

- `reason`: Optional description of why this region is excluded

**Sample-specific exclusions:** If the `regions` data frame includes a
`sample_id` column, exclusions can be applied only to specific samples.
Rows without `sample_id` (or with `sample_id = NA`) are applied globally
to all samples.

**Common uses:**

- Solvent peaks that elute at the end of the chromatogram

- System peaks or artifacts at specific retention times

- Flow marker peaks used for retention time correction

- Air bubbles or other transient artifacts

## Tidying

When you
[`tidy()`](https://recipes.tidymodels.org/reference/tidy.recipe.html)
this step, a tibble with columns `terms`, `n_regions`, `purpose`, and
`id` is returned.

## See also

[`step_sec_baseline()`](https://jameshwade.github.io/measure-sec/reference/step_sec_baseline.md)
which can use exclusion information,
[`step_sec_mw_averages()`](https://jameshwade.github.io/measure-sec/reference/step_sec_mw_averages.md)
which respects integration exclusions.

Other sec-chromatography:
[`step_sec_band_broadening()`](https://jameshwade.github.io/measure-sec/reference/step_sec_band_broadening.md),
[`step_sec_baseline()`](https://jameshwade.github.io/measure-sec/reference/step_sec_baseline.md),
[`step_sec_detector_delay()`](https://jameshwade.github.io/measure-sec/reference/step_sec_detector_delay.md),
[`step_sec_integration_window()`](https://jameshwade.github.io/measure-sec/reference/step_sec_integration_window.md),
[`step_sec_mw_averages()`](https://jameshwade.github.io/measure-sec/reference/step_sec_mw_averages.md),
[`step_sec_mw_distribution()`](https://jameshwade.github.io/measure-sec/reference/step_sec_mw_distribution.md),
[`step_sec_mw_fractions()`](https://jameshwade.github.io/measure-sec/reference/step_sec_mw_fractions.md),
[`step_sec_peaks_deconvolve()`](https://jameshwade.github.io/measure-sec/reference/step_sec_peaks_deconvolve.md),
[`step_sec_peaks_detect()`](https://jameshwade.github.io/measure-sec/reference/step_sec_peaks_detect.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(recipes)
library(measure)

# Exclude a single region (solvent peak at end)
exclusions <- tibble::tibble(
  start = 18.5,
  end = 20.0,
  reason = "Solvent peak"
)

rec <- recipe(~., data = sec_data) |>
  step_measure_input_long(ri_signal, location = vars(elution_volume)) |>
  step_sec_exclude_regions(regions = exclusions) |>
  prep()

# Exclude multiple regions
exclusions <- tibble::tibble(
  start = c(8.0, 18.5),
  end = c(9.0, 20.0),
  reason = c("Void volume", "Solvent peak")
)

rec <- recipe(~., data = sec_data) |>
  step_measure_input_long(ri_signal, location = vars(elution_volume)) |>
  step_sec_exclude_regions(regions = exclusions, purpose = "integration") |>
  prep()

# Sample-specific exclusions
exclusions <- tibble::tibble(
  start = c(18.5, 15.0),
  end = c(20.0, 16.0),
  sample_id = c(NA, "sample_2"),  # NA = global, specific sample_id = per-sample
  reason = c("Solvent peak", "Artifact in sample_2")
)

rec <- recipe(~., data = sec_data) |>
  step_measure_input_long(ri_signal, location = vars(elution_volume)) |>
  step_sec_exclude_regions(regions = exclusions) |>
  prep()
} # }
```
