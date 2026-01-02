# Plot SEC Conformation Data

Creates conformation plots showing structure-MW relationships from
multi-detector SEC data (Mark-Houwink plots, Rg-MW scaling).

## Usage

``` r
plot_sec_conformation(
  data,
  type = c("rg_mw", "eta_mw", "rh_mw"),
  mw_col = "mw",
  y_col = NULL,
  sample_id = NULL,
  show_fit = TRUE,
  show_exponent = TRUE,
  compare_linear = NULL,
  ...
)
```

## Arguments

- data:

  A data frame containing SEC results with MW and conformation data (Rg,
  intrinsic viscosity).

- type:

  Type of conformation plot:

  - `"rg_mw"`: Radius of gyration vs molecular weight (default)

  - `"eta_mw"`: Intrinsic viscosity vs molecular weight

  - `"rh_mw"`: Hydrodynamic radius vs molecular weight

- mw_col:

  Name of molecular weight column. Default is `"mw"`.

- y_col:

  Name of y-axis column. If `NULL`, auto-detected based on type.

- sample_id:

  Column name containing sample identifiers.

- show_fit:

  Logical. Show power-law fit line? Default is `TRUE`.

- show_exponent:

  Logical. Annotate slope/exponent on plot? Default is `TRUE`.

- compare_linear:

  Data frame with linear reference polymer for branching comparison.
  Should have same structure as main data.

- ...:

  Additional arguments passed to
  [`ggplot2::geom_point()`](https://ggplot2.tidyverse.org/reference/geom_point.html).

## Value

A ggplot2 object.

## Details

Conformation plots reveal polymer architecture:

**Rg-MW (radius of gyration):** Log-log slope indicates conformation:

- 0.33: Compact sphere

- 0.5-0.6: Random coil (linear polymer in good solvent)

- 1.0: Rigid rod Branched polymers show reduced Rg at same MW compared
  to linear.

**eta-MW (Mark-Houwink):** The relationship `[eta] = K * M^a` where:

- K and a are Mark-Houwink parameters

- a is approximately 0.5-0.8 for typical polymers

- Branched polymers show lower `[eta]` at same MW

## See also

Other sec-visualization:
[`autoplot.sec_results()`](https://jameshwade.github.io/measure-sec/reference/autoplot.sec_results.md),
[`plot_sec()`](https://jameshwade.github.io/measure-sec/reference/plot_sec.md),
[`plot_sec_calibration()`](https://jameshwade.github.io/measure-sec/reference/plot_sec_calibration.md),
[`plot_sec_chromatogram()`](https://jameshwade.github.io/measure-sec/reference/plot_sec_chromatogram.md),
[`plot_sec_composition()`](https://jameshwade.github.io/measure-sec/reference/plot_sec_composition.md),
[`plot_sec_multidetector()`](https://jameshwade.github.io/measure-sec/reference/plot_sec_multidetector.md),
[`plot_sec_mwd()`](https://jameshwade.github.io/measure-sec/reference/plot_sec_mwd.md),
[`sec_results()`](https://jameshwade.github.io/measure-sec/reference/sec_results.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Rg-MW plot from MALS data
plot_sec_conformation(processed_data, type = "rg_mw")

# Mark-Houwink plot
plot_sec_conformation(processed_data, type = "eta_mw")

# Compare branched to linear reference
plot_sec_conformation(
  branched_data,
  compare_linear = linear_reference
)
} # }
```
