# Fit Rg-MW Scaling Relationship

Fits the power law relationship between radius of gyration (Rg) and
molecular weight (MW) to determine polymer conformation.

## Usage

``` r
measure_rg_mw_scaling(mw, rg, weights = NULL, mw_range = NULL)
```

## Arguments

- mw:

  Numeric vector of molecular weights.

- rg:

  Numeric vector of radii of gyration (same length as `mw`).

- weights:

  Optional numeric vector of weights (e.g., concentration).

- mw_range:

  Optional numeric vector of length 2 specifying MW range for fitting.

## Value

A list of class `rg_mw_scaling` containing:

- nu:

  Scaling exponent (slope in log-log space)

- prefactor:

  Prefactor K in Rg = K \* M^nu

- r_squared:

  R-squared of the fit

- conformation:

  Interpreted polymer conformation

- n_points:

  Number of data points used

- fit:

  The lm fit object

## Details

The Rg-MW relationship follows a power law:

\$\$R_g = K \cdot M^{\nu}\$\$

where nu (the Flory exponent) indicates polymer conformation:

**Interpretation of nu:**

- nu ~ 0.33: Compact/spherical (collapsed globule)

- nu ~ 0.50: Theta solvent (ideal chain)

- nu ~ 0.588: Good solvent (swollen coil)

- nu ~ 1.0: Rigid rod

**Typical Values:**

- Flexible polymers in good solvent: nu = 0.55-0.60

- Branched polymers: nu = 0.40-0.50

- Proteins (globular): nu = 0.30-0.35

- DNA: nu = 0.58-0.60

## References

Flory, P.J. (1953). "Principles of Polymer Chemistry." Cornell
University Press.

## See also

Other sec-polymer:
[`measure_branching_frequency()`](https://jameshwade.github.io/measure-sec/reference/measure_branching_frequency.md),
[`measure_branching_index()`](https://jameshwade.github.io/measure-sec/reference/measure_branching_index.md),
[`measure_branching_model_comparison()`](https://jameshwade.github.io/measure-sec/reference/measure_branching_model_comparison.md),
[`measure_conformation_data()`](https://jameshwade.github.io/measure-sec/reference/measure_conformation_data.md),
[`measure_mh_parameters()`](https://jameshwade.github.io/measure-sec/reference/measure_mh_parameters.md)

## Examples

``` r
# Fit Rg-MW for a linear polymer
mw <- c(10000, 50000, 100000, 500000, 1000000)
rg <- c(4.5, 12, 18, 45, 70)  # nm

scaling <- measure_rg_mw_scaling(mw, rg)
print(scaling)
#> Rg-MW Scaling Analysis
#> ======================================== 
#> 
#> Scaling Law: Rg = K * M^nu
#> 
#> nu (Flory exponent): 0.591
#>   95% CI: [0.574, 0.609]
#> K (prefactor): 1.9661e-02
#> 
#> Conformation: good solvent (swollen coil)
#> 
#> R-squared: 0.9997
#> Data points: 5
#> MW range: 10000 - 1000000
```
