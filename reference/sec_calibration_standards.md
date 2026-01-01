# SEC Calibration Standards - Complete Dataset

Comprehensive calibration standard data for SEC/GPC analysis, containing
polystyrene and PMMA narrow standards with certificate values, retention
data, and Mark-Houwink parameters. Based on commercial standard kits.

## Usage

``` r
sec_calibration_standards
```

## Format

A tibble with 26 rows and 19 columns:

- standard_name:

  Character. Standard identifier (e.g., "PS-67500", "PMMA-30300")

- polymer_type:

  Character. Either "polystyrene" or "pmma"

- kit_name:

  Character. Commercial kit name for reference

- mp:

  Numeric. Peak molecular weight in Da (most commonly used for
  calibration)

- mn:

  Numeric. Number-average molecular weight in Da

- mw:

  Numeric. Weight-average molecular weight in Da

- dispersity:

  Numeric. Polydispersity index (Mw/Mn), typically 1.02-1.05

- mp_uncertainty:

  Numeric. Relative uncertainty in Mp (e.g., 0.05 = 5%)

- log_mp:

  Numeric. log10(Mp) for calibration curve fitting

- log_mw:

  Numeric. log10(Mw)

- log_mn:

  Numeric. log10(Mn)

- retention_time:

  Numeric. Peak retention time in minutes

- retention_volume:

  Numeric. Peak retention volume in mL

- k_value:

  Numeric. Mark-Houwink K constant in mL/g

- a_value:

  Numeric. Mark-Houwink exponent (alpha)

- intrinsic_viscosity:

  Numeric. Intrinsic viscosity in mL/g

- log_hydrodynamic_vol:

  Numeric. log10(M \* \[eta\]) for universal calibration

- dn_dc:

  Numeric. Refractive index increment in mL/g

- notes:

  Character. Special notes (e.g., near exclusion limit)

## Source

Synthetic data based on typical commercial narrow standards (e.g.,
Agilent EasiVial, PSS ReadyCal) with realistic retention times for PLgel
Mixed-C columns in THF at 1.0 mL/min.

## Details

This dataset enables several key calibration workflows:

**Conventional Calibration:** Use `retention_time` (or
`retention_volume`) and `log_mp` with
[`step_sec_conventional_cal`](step_sec_conventional_cal.md) to build a
calibration curve.

**Universal Calibration:** Use `log_hydrodynamic_vol` vs retention for
polymer-independent calibration. The Mark-Houwink parameters (`k_value`,
`a_value`) enable conversion between polymers.

**Quality Assessment:** The `mp_uncertainty` values (from typical
certificates) enable uncertainty propagation through the calibration.
Dispersity values confirm standards are suitably narrow.

**Polymer Types:**

- 16 polystyrene standards (162 Da to 3,150,000 Da)

- 10 PMMA standards (602 Da to 1,190,000 Da)

**Mark-Houwink Parameters (THF, 35°C):**

- Polystyrene: K = 0.000141 mL/g, a = 0.700

- PMMA: K = 0.000128 mL/g, a = 0.690

## See also

[`sec_ps_standards`](sec_ps_standards.md) for polystyrene-only subset
[`sec_pmma_standards`](sec_pmma_standards.md) for PMMA-only subset
[`step_sec_conventional_cal`](step_sec_conventional_cal.md) for
conventional calibration
[`step_sec_universal_cal`](step_sec_universal_cal.md) for universal
calibration

Other sec-data: [`sec_pmma_standards`](sec_pmma_standards.md),
[`sec_ps_standards`](sec_ps_standards.md),
[`sec_triple_detect`](sec_triple_detect.md)

Other sec-calibration: [`sec_pmma_standards`](sec_pmma_standards.md),
[`sec_ps_standards`](sec_ps_standards.md),
[`step_sec_conventional_cal()`](step_sec_conventional_cal.md),
[`step_sec_universal_cal()`](step_sec_universal_cal.md)

## Examples

``` r
library(dplyr)
data(sec_calibration_standards)

# View the polystyrene standards
sec_calibration_standards |>
  filter(polymer_type == "polystyrene") |>
  select(standard_name, mp, log_mp, retention_time)
#> # A tibble: 16 × 4
#>    standard_name      mp log_mp retention_time
#>    <chr>           <dbl>  <dbl>          <dbl>
#>  1 PS-3150000    3150000   6.50           11.2
#>  2 PS-1870000    1870000   6.27           11.6
#>  3 PS-1090000    1090000   6.04           12.1
#>  4 PS-630000      630000   5.80           12.6
#>  5 PS-430000      430000   5.63           13.2
#>  6 PS-216000      216000   5.33           13.8
#>  7 PS-120000      120000   5.08           14.3
#>  8 PS-67500        67500   4.83           15.0
#>  9 PS-33500        33500   4.53           15.5
#> 10 PS-19800        19800   4.30           16.1
#> 11 PS-9680          9680   3.99           16.7
#> 12 PS-5030          5030   3.70           17.4
#> 13 PS-2970          2970   3.47           17.9
#> 14 PS-1050          1050   3.02           18.9
#> 15 PS-580            580   2.76           19.4
#> 16 PS-162            162   2.21           20.8

# Compare PS and PMMA at similar MW - PMMA elutes later (smaller Rh)
sec_calibration_standards |>
  filter(mp > 60000 & mp < 80000) |>
  select(standard_name, polymer_type, mp, retention_time)
#> # A tibble: 2 × 4
#>   standard_name polymer_type    mp retention_time
#>   <chr>         <chr>        <dbl>          <dbl>
#> 1 PMMA-67700    pmma         67700           14.8
#> 2 PS-67500      polystyrene  67500           15.0

# Prepare standards for conventional calibration
ps_cal <- sec_calibration_standards |>
  filter(polymer_type == "polystyrene") |>
  select(retention = retention_time, log_mw = log_mp)
```
