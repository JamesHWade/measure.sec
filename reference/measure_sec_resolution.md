# Calculate Peak Resolution

Calculates the resolution between two chromatographic peaks using the
USP or EP formula.

## Usage

``` r
measure_sec_resolution(
  retention_1,
  retention_2,
  width_1,
  width_2,
  width_type = c("baseline", "half_height", "tangent"),
  method = c("usp", "ep")
)
```

## Arguments

- retention_1:

  Retention time of the first peak (earlier eluting).

- retention_2:

  Retention time of the second peak (later eluting).

- width_1:

  Peak width of the first peak. See `width_type` for units.

- width_2:

  Peak width of the second peak.

- width_type:

  Type of peak width measurement:

  - `"baseline"` (default): Width at baseline (Wb)

  - `"half_height"`: Width at half height (W0.5h)

  - `"tangent"`: Width from tangent lines at inflection points

- method:

  Resolution formula to use:

  - `"usp"` (default): Rs = 2(t2 - t1) / (w1 + w2)

  - `"ep"`: Rs = 1.18(t2 - t1) / (w1_0.5h + w2_0.5h)

## Value

Numeric resolution value. Rs \> 1.5 indicates baseline separation.

## Details

Resolution quantifies the degree of separation between adjacent peaks:

**USP Formula (baseline width):** \$\$R_s = \frac{2(t_2 -
t_1)}{W\_{b1} + W\_{b2}}\$\$

**EP Formula (half-height width):** \$\$R_s = \frac{1.18(t_2 -
t_1)}{W\_{0.5h,1} + W\_{0.5h,2}}\$\$

**Interpretation:**

- Rs \< 1.0: Peaks overlap significantly

- Rs = 1.0: ~94% separation (4 sigma)

- Rs = 1.5: Baseline separation (~99.7%)

- Rs \> 2.0: Complete separation with gap

## See also

Other sec-qc:
[`measure_sec_asymmetry()`](https://jameshwade.github.io/measure-sec/reference/measure_sec_asymmetry.md),
[`measure_sec_plate_count()`](https://jameshwade.github.io/measure-sec/reference/measure_sec_plate_count.md),
[`measure_sec_recovery()`](https://jameshwade.github.io/measure-sec/reference/measure_sec_recovery.md),
[`measure_sec_suitability()`](https://jameshwade.github.io/measure-sec/reference/measure_sec_suitability.md)

## Examples

``` r
# Calculate resolution between monomer and dimer
measure_sec_resolution(
  retention_1 = 8.2,   # dimer (elutes first in SEC)
  retention_2 = 9.5,   # monomer
  width_1 = 0.4,
  width_2 = 0.5
)
#> [1] 2.888889
```
