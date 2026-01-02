# Plot SEC Composition Distribution

Creates a composition plot showing copolymer or blend composition as a
function of molecular weight or elution time.

## Usage

``` r
plot_sec_composition(
  data,
  composition_col = "composition_a",
  x_axis = c("mw", "retention"),
  mw_col = "mw",
  sample_id = NULL,
  show_average = TRUE,
  component_names = NULL,
  show_distribution = TRUE,
  show_points = FALSE,
  y_limits = c(0, 1),
  log_mw = TRUE,
  ...
)
```

## Arguments

- data:

  A data frame containing SEC results with composition data, typically
  from
  [`step_sec_composition()`](https://jameshwade.github.io/measure-sec/reference/step_sec_composition.md).

- composition_col:

  Name of the composition measure column. Default is `"composition_a"`.
  Should contain values from 0 to 1.

- x_axis:

  Variable for x-axis: `"mw"` (molecular weight, default) or
  `"retention"` (elution time).

- mw_col:

  Name of molecular weight column. Default is `"mw"`.

- sample_id:

  Column name containing sample identifiers.

- show_average:

  Logical. Show average composition line? Default is `TRUE`.

- component_names:

  Named character vector with component labels. E.g.,
  `c(a = "Styrene", b = "Acrylate")`. Default shows "Component A" and
  "Component B".

- show_distribution:

  Logical. Show composition distribution as lines? Default is `TRUE`.

- show_points:

  Logical. Show individual data points? Default is `FALSE`.

- y_limits:

  Numeric vector of length 2 for y-axis limits. Default is `c(0, 1)` for
  fraction scale.

- log_mw:

  Logical. Use log scale for MW on x-axis? Default is `TRUE` when
  `x_axis = "mw"`.

- ...:

  Additional arguments passed to
  [`ggplot2::geom_line()`](https://ggplot2.tidyverse.org/reference/geom_path.html).

## Value

A ggplot2 object.

## Details

Composition plots are essential for characterizing copolymers and
blends:

**Uniform copolymers:** Horizontal line across all MW values indicates
consistent composition throughout the distribution.

**Compositional drift:** Slope indicates composition varies with chain
length, common in batch polymerization where monomer ratios change over
time.

**Blend separation:** Multiple distinct compositions indicate blend
separation or block copolymer structure.

The plot works with output from
[`step_sec_composition()`](https://jameshwade.github.io/measure-sec/reference/step_sec_composition.md)
which calculates weight fraction from UV/RI detector signals and known
response factors.

## See also

Other sec-visualization:
[`autoplot.sec_results()`](https://jameshwade.github.io/measure-sec/reference/autoplot.sec_results.md),
[`plot_sec()`](https://jameshwade.github.io/measure-sec/reference/plot_sec.md),
[`plot_sec_calibration()`](https://jameshwade.github.io/measure-sec/reference/plot_sec_calibration.md),
[`plot_sec_chromatogram()`](https://jameshwade.github.io/measure-sec/reference/plot_sec_chromatogram.md),
[`plot_sec_conformation()`](https://jameshwade.github.io/measure-sec/reference/plot_sec_conformation.md),
[`plot_sec_multidetector()`](https://jameshwade.github.io/measure-sec/reference/plot_sec_multidetector.md),
[`plot_sec_mwd()`](https://jameshwade.github.io/measure-sec/reference/plot_sec_mwd.md),
[`sec_results()`](https://jameshwade.github.io/measure-sec/reference/sec_results.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(recipes)
library(measure)
library(ggplot2)

# Process copolymer data with composition calculation
processed <- recipe(~., data = sec_copolymer) |>
  step_measure_input_long(ri_signal, location = vars(elution_time), col_name = "ri") |>
  step_measure_input_long(uv_254_signal, location = vars(elution_time), col_name = "uv") |>
  step_sec_baseline() |>
  step_sec_composition(
    uv_col = "uv",
    ri_col = "ri",
    component_a_uv = 1.0,
    component_a_ri = 0.185,
    component_b_uv = 0.01,
    component_b_ri = 0.084,
    output_col = "styrene_frac"
  ) |>
  step_sec_conventional_cal(standards = ps_standards) |>
  prep() |>
  bake(new_data = NULL)

# Basic composition plot vs MW
plot_sec_composition(
  processed,
  composition_col = "styrene_frac"
)

# Composition vs retention time with custom component names
plot_sec_composition(
  processed,
  composition_col = "styrene_frac",
  x_axis = "retention",
  component_names = c(a = "Styrene", b = "Acrylate")
)

# Customize appearance
plot_sec_composition(processed, composition_col = "styrene_frac") +
  theme_bw() +
  labs(title = "Styrene-Acrylate Copolymer Composition")
} # }
```
