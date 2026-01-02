# Compare Branching Models

Compares experimental branching data against theoretical predictions
from different branching models to identify the most likely
architecture.

## Usage

``` r
measure_branching_model_comparison(
  g,
  mw,
  g_prime = NULL,
  models = c("random", "star", "comb", "hyperbranched"),
  branch_frequency = NULL
)
```

## Arguments

- g:

  Numeric vector of experimental g ratios (Rg^2 branched / Rg^2 linear).

- mw:

  Numeric vector of molecular weights (same length as `g`).

- g_prime:

  Optional numeric vector of g' ratios (IV branched / IV linear).

- models:

  Character vector of models to compare. Default compares all:

  - `"random"`: Random trifunctional branching (Zimm-Stockmayer)

  - `"star"`: Star polymers (variable arms)

  - `"comb"`: Comb architecture

  - `"hyperbranched"`: Hyperbranched/dendritic

- branch_frequency:

  For fitting: assumed constant branching frequency (branches per 1000
  Da). If NULL, estimated from data.

## Value

A list of class `branching_model_comparison` containing:

- model_fits:

  Data frame with fit statistics for each model

- best_model:

  Name of the best-fitting model

- predictions:

  Data frame with predicted g for each model

- experimental:

  Input experimental data

## Details

**Model Equations:**

**Random trifunctional (Zimm-Stockmayer):** \$\$g = \left\[\left(1 +
\frac{n_b}{7}\right)^{1/2} + \frac{4n_b}{9\pi}\right\]^{-1/2}\$\$

**Star polymer (f arms):** \$\$g = \frac{3f - 2}{f^2}\$\$

**Comb polymer (n_b branches):** \$\$g \approx \frac{1}{1 + 2n_b/3}\$\$

**Hyperbranched (degree of branching DB):** \$\$g \approx
\left(\frac{1}{1 + DB \cdot n/2}\right)^{0.5}\$\$

Model selection uses residual sum of squares and AIC-like criteria.

## References

Zimm, B.H. and Stockmayer, W.H. (1949). J. Chem. Phys., 17, 1301-1314.

Burchard, W. (1999). "Solution Properties of Branched Macromolecules."
Adv. Polym. Sci., 143, 113-194.

## See also

Other sec-polymer:
[`measure_branching_frequency()`](https://jameshwade.github.io/measure-sec/reference/measure_branching_frequency.md),
[`measure_branching_index()`](https://jameshwade.github.io/measure-sec/reference/measure_branching_index.md),
[`measure_conformation_data()`](https://jameshwade.github.io/measure-sec/reference/measure_conformation_data.md),
[`measure_mh_parameters()`](https://jameshwade.github.io/measure-sec/reference/measure_mh_parameters.md),
[`measure_rg_mw_scaling()`](https://jameshwade.github.io/measure-sec/reference/measure_rg_mw_scaling.md)

## Examples

``` r
# Compare models for branched polymer data
mw <- c(50000, 100000, 200000, 500000)
g_exp <- c(0.85, 0.72, 0.58, 0.42)

comparison <- measure_branching_model_comparison(g_exp, mw)
print(comparison)
#> Branching Model Comparison
#> ============================================================ 
#> 
#> Model Fit Summary:
#> ------------------------------------------------------------ 
#>          model  parameter r_squared       rmse       aic
#>         random 0.04648146 0.9887283 0.01699313 -17.24806
#>           star 1.11325520 0.9470407 0.03683411 -11.05914
#>           comb 0.00502421 0.9670907 0.02903611 -12.96221
#>  hyperbranched 0.09231394 0.9932768 0.01312403 -19.31498
#> ------------------------------------------------------------ 
#> 
#> Best Model: hyperbranched (lowest AIC)
#>   R-squared: 0.9933
#>   RMSE: 0.0131
#>   Parameter: 0.0923
#> 
#> Note: Lower AIC indicates better fit with penalty for complexity.
#> Consider physical plausibility when selecting the final model.
```
