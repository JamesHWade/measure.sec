# Calculate Branching Frequency from Zimm-Stockmayer Theory

Calculates the number of branch points per molecule using
Zimm-Stockmayer theory, which relates the branching index (g) to the
branching frequency.

## Usage

``` r
measure_branching_frequency(
  g,
  architecture = c("random", "star_3", "star_4", "star_f", "comb"),
  arms = NULL,
  mw = NULL
)
```

## Arguments

- g:

  Numeric vector of branching ratios g = Rg^2(branched) / Rg^2(linear).

- architecture:

  Branching architecture model:

  - `"random"`: Random trifunctional branching (default)

  - `"star_3"`: 3-arm star polymer

  - `"star_4"`: 4-arm star polymer

  - `"star_f"`: f-arm star (requires `arms` parameter)

  - `"comb"`: Comb-like branching

- arms:

  Number of arms for star polymers (required for `"star_f"`).

- mw:

  Optional molecular weight vector (same length as `g`) for calculating
  total branches in sample.

## Value

A data frame of class `branching_frequency` containing:

- g:

  Input branching ratio

- branches_per_molecule:

  Estimated branch points per molecule

- branch_density:

  Branches per unit MW (if `mw` provided)

## Details

**Zimm-Stockmayer Theory:**

For random trifunctional branching, the relationship between g and the
average number of branch points (n_b) is:

\$\$g = \left\[\left(1 + \frac{n_b}{7}\right)^{1/2} +
\frac{4n_b}{9\pi}\right\]^{-1/2}\$\$

For star polymers with f arms: \$\$g = \frac{3f - 2}{f^2}\$\$

**Interpretation:**

- g = 1.0: Linear polymer (no branching)

- g ~ 0.5-0.8: Lightly branched

- g ~ 0.3-0.5: Moderately branched

- g \< 0.3: Highly branched/hyperbranched

## References

Zimm, B.H. and Stockmayer, W.H. (1949). "The Dimensions of Chain
Molecules Containing Branches and Rings." J. Chem. Phys., 17, 1301-1314.

## See also

Other sec-polymer:
[`measure_branching_index()`](https://jameshwade.github.io/measure-sec/reference/measure_branching_index.md),
[`measure_branching_model_comparison()`](https://jameshwade.github.io/measure-sec/reference/measure_branching_model_comparison.md),
[`measure_conformation_data()`](https://jameshwade.github.io/measure-sec/reference/measure_conformation_data.md),
[`measure_mh_parameters()`](https://jameshwade.github.io/measure-sec/reference/measure_mh_parameters.md),
[`measure_rg_mw_scaling()`](https://jameshwade.github.io/measure-sec/reference/measure_rg_mw_scaling.md)

## Examples

``` r
# Calculate branching for random trifunctional polymer
g_values <- c(0.9, 0.7, 0.5, 0.3)
bf <- measure_branching_frequency(g_values)
print(bf)
#> Branching Frequency Analysis (Zimm-Stockmayer)
#> ================================================== 
#> 
#> Architecture: random
#> N samples: 4
#> 
#> g ratio range: 0.300 - 0.900
#> Branches/molecule: 1.12 - 57.14 (mean: 19.74)
#> 
#> # A tibble: 4 × 2
#>       g branches_per_molecule
#>   <dbl>                 <dbl>
#> 1   0.9                  1.12
#> 2   0.7                  5.12
#> 3   0.5                 15.6 
#> 4   0.3                 57.1 

# For star polymer with 4 arms
bf_star <- measure_branching_frequency(0.625, architecture = "star_4")
```
