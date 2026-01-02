# Calculate Peak Asymmetry Factor

Calculates the asymmetry factor (As) or tailing factor (Tf) for a
chromatographic peak.

## Usage

``` r
measure_sec_asymmetry(leading, tailing, method = c("usp", "ep"))
```

## Arguments

- leading:

  Width of the leading (front) half of the peak at the measurement
  height.

- tailing:

  Width of the tailing (back) half of the peak at the measurement
  height.

- method:

  Asymmetry calculation method:

  - `"usp"` (default): Tailing factor at 5% height

  - `"ep"`: Asymmetry factor at 10% height

## Value

Numeric asymmetry value. Values \> 1 indicate tailing, \< 1 indicate
fronting. Ideal value is 1.0.

## Details

Peak asymmetry indicates deviation from ideal Gaussian peak shape:

**USP Tailing Factor (at 5% height):** \$\$T_f =
\frac{W\_{0.05}}{2f}\$\$

where W_0.05 is the width at 5% height and f is the leading half-width.

**EP Asymmetry Factor (at 10% height):** \$\$A_s = \frac{b}{a}\$\$

where b is the tailing half-width and a is the leading half-width.

**Interpretation:**

- As = 1.0: Symmetric (ideal)

- As \< 0.9 or \> 1.2: Slight asymmetry (acceptable)

- As \< 0.8 or \> 1.5: Significant asymmetry (investigate)

- As \> 2.0: Severe tailing (column/sample issue)

## See also

Other sec-qc:
[`measure_sec_column_performance()`](https://jameshwade.github.io/measure-sec/reference/measure_sec_column_performance.md),
[`measure_sec_plate_count()`](https://jameshwade.github.io/measure-sec/reference/measure_sec_plate_count.md),
[`measure_sec_recovery()`](https://jameshwade.github.io/measure-sec/reference/measure_sec_recovery.md),
[`measure_sec_resolution()`](https://jameshwade.github.io/measure-sec/reference/measure_sec_resolution.md),
[`measure_sec_suitability()`](https://jameshwade.github.io/measure-sec/reference/measure_sec_suitability.md)

## Examples

``` r
# Calculate USP tailing factor
measure_sec_asymmetry(
  leading = 0.12,
  tailing = 0.15,
  method = "usp"
)
#> [1] 1.125
```
