# Create Conformation Plot Data

Prepares data for Mark-Houwink or conformation plots (log(\[eta\]) vs
log(M) or log(Rg) vs log(M)).

## Usage

``` r
measure_conformation_data(mw, y, y_type = c("iv", "rg"), fit_line = TRUE)
```

## Arguments

- mw:

  Numeric vector of molecular weights.

- y:

  Numeric vector of intrinsic viscosity or radius of gyration.

- y_type:

  Type of y-axis data: `"iv"` for intrinsic viscosity or `"rg"` for
  radius of gyration.

- fit_line:

  Logical. Include fitted line? Default is `TRUE`.

## Value

A data frame suitable for plotting with columns:

- log_mw:

  log10(MW)

- log_y:

  log10(y)

- mw:

  Original MW values

- y:

  Original y values

## See also

Other sec-polymer:
[`measure_branching_index()`](measure_branching_index.md),
[`measure_mh_parameters()`](measure_mh_parameters.md)

## Examples

``` r
mw <- c(10000, 50000, 100000, 500000)
iv <- c(0.15, 0.35, 0.50, 0.95)

plot_data <- measure_conformation_data(mw, iv, y_type = "iv")

# Plot with ggplot2
# ggplot(plot_data, aes(log_mw, log_y)) +
#   geom_point() +
#   geom_smooth(method = "lm")
```
