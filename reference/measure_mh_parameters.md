# Estimate Mark-Houwink Parameters

Estimates the Mark-Houwink K and a (alpha) parameters from intrinsic
viscosity and molecular weight data.

## Usage

``` r
measure_mh_parameters(
  mw,
  intrinsic_visc,
  weights = NULL,
  mw_range = NULL,
  log_fit = TRUE
)
```

## Arguments

- mw:

  Numeric vector of molecular weights.

- intrinsic_visc:

  Numeric vector of intrinsic viscosities (same length as `mw`).

- weights:

  Optional numeric vector of weights for weighted regression. Use
  concentration or signal intensity as weights for SEC data.

- mw_range:

  Optional numeric vector of length 2 specifying the MW range to use for
  fitting. Data outside this range is excluded.

- log_fit:

  Logical. Perform fit in log-log space (recommended)? Default is
  `TRUE`.

## Value

A list of class `mh_parameters` containing:

- K:

  Mark-Houwink K parameter

- a:

  Mark-Houwink a (alpha) exponent

- r_squared:

  R-squared of the fit

- n_points:

  Number of data points used

- mw_range:

  MW range of the data

- fit:

  The fitted linear model object

## Details

The Mark-Houwink equation relates intrinsic viscosity to molecular
weight:

\$\$\[\eta\] = K \cdot M^a\$\$

In log form: \$\$\log(\[\eta\]) = \log(K) + a \cdot \log(M)\$\$

The parameters K and a depend on:

- Polymer-solvent system

- Temperature

- Polymer microstructure (tacticity, branching)

**Interpretation of 'a' exponent:**

- a ~ 0.5: Theta solvent (polymer coil collapsed)

- a ~ 0.5-0.8: Good solvent (typical range)

- a ~ 0.8: Rigid rod or extended chain

- a \< 0.5: Branched or compact structures

## See also

Other sec-polymer:
[`measure_branching_index()`](https://jameshwade.github.io/measure-sec/reference/measure_branching_index.md),
[`measure_conformation_data()`](https://jameshwade.github.io/measure-sec/reference/measure_conformation_data.md)

## Examples

``` r
# Estimate Mark-Houwink parameters from triple-detection data
mw <- c(10000, 25000, 50000, 100000, 250000)
iv <- c(0.15, 0.28, 0.45, 0.72, 1.2)

mh <- measure_mh_parameters(mw, iv)
print(mh)
#> Mark-Houwink Parameters
#> ======================================== 
#> 
#> K = 3.8112e-04
#> a = 0.652
#> 
#> R-squared: 0.9981
#> Data points: 5
#> MW range: 10000 - 250000
#> 
#> Equation: [eta] = K * M^a
# K = 0.000114, a = 0.716 (typical for PS in THF)
```
