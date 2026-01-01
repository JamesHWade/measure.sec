# Calculate Branching Index

Calculates the branching index (g or g') for branched polymers by
comparing their properties to linear reference polymers of the same
molecular weight.

## Usage

``` r
measure_branching_index(
  mw,
  rg = NULL,
  intrinsic_visc = NULL,
  reference = "linear",
  mh_linear = NULL,
  rg_linear_fit = NULL,
  method = c("g", "g_prime", "both")
)
```

## Arguments

- mw:

  Numeric vector of molecular weights for branched samples.

- rg:

  Numeric vector of radius of gyration values (for g ratio).

- intrinsic_visc:

  Numeric vector of intrinsic viscosity values (for g').

- reference:

  Type of reference: either `"linear"` for comparison to linear polymer
  theory, or a data frame containing linear reference data.

- mh_linear:

  Mark-Houwink parameters for linear reference polymer (list with K and
  a). Required if `reference = "linear"` and using g'.

- rg_linear_fit:

  Rg-MW relationship for linear polymer: either a model or a list with
  `slope` and `intercept` for log(Rg) = intercept + slope \* log(M).

- method:

  Branching index type:

  - `"g"`: Radius of gyration ratio: g = Rg^2(branched) / Rg^2(linear)

  - `"g_prime"`: Viscosity ratio: g' = \[eta\](branched) /
    \[eta\](linear)

  - `"both"`: Calculate both g and g'

## Value

A data frame with columns:

- mw:

  Molecular weight

- g:

  Branching index from Rg (if calculated)

- g_prime:

  Branching index from viscosity (if calculated)

- branches_per_molecule:

  Estimated branch points (Zimm-Stockmayer)

## Details

Branching reduces the hydrodynamic size of polymers compared to their
linear counterparts. This is quantified by the branching ratios:

**Rg-based branching index (g):** \$\$g =
\frac{R_g^2(branched)}{R_g^2(linear)}\$\$

**Viscosity-based branching index (g'):** \$\$g' =
\frac{\[\eta\](branched)}{\[\eta\](linear)}\$\$

The relationship between g and g' depends on polymer architecture:
\$\$g' = g^\epsilon\$\$

where epsilon ~ 0.5-1.5 depending on branching type.

**Zimm-Stockmayer Model (random branching):** \$\$g = \frac{6}{n_b}
\left\[ \frac{1}{2} + \frac{(2 + n_b)^{1/2} - 1 - n_b/2}{n_b}
\right\]\$\$

where n_b is the number of branch points per molecule.

## See also

Other sec-polymer:
[`measure_conformation_data()`](measure_conformation_data.md),
[`measure_mh_parameters()`](measure_mh_parameters.md)

## Examples

``` r
# Calculate branching index from Rg data
mw <- c(100000, 200000, 500000)
rg_branched <- c(12, 18, 32)
rg_linear <- c(15, 22, 40)  # Reference linear polymer

g <- measure_branching_index(
  mw = mw,
  rg = rg_branched,
  reference = data.frame(mw = mw, rg = rg_linear),
  method = "g"
)
```
