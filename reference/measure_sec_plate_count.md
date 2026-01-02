# Calculate Theoretical Plate Count

Calculates the number of theoretical plates (N) for a chromatographic
peak, a measure of column efficiency.

## Usage

``` r
measure_sec_plate_count(
  retention,
  width,
  width_type = c("half_height", "baseline", "inflection"),
  dead_time = NULL
)
```

## Arguments

- retention:

  Retention time of the peak.

- width:

  Peak width. See `width_type` for measurement method.

- width_type:

  Type of peak width measurement:

  - `"half_height"` (default): Width at 50% height (W0.5h)

  - `"baseline"`: Width at baseline from tangent lines (Wb)

  - `"inflection"`: Width at inflection points (Wi)

- dead_time:

  Column dead time (t0). If provided, calculates effective plates
  (N_eff) using adjusted retention time.

## Value

Numeric plate count. Higher values indicate better efficiency.

## Details

Theoretical plate count measures column efficiency:

**Half-height width (most common):** \$\$N = 5.54
\left(\frac{t_R}{W\_{0.5h}}\right)^2\$\$

**Baseline width:** \$\$N = 16 \left(\frac{t_R}{W_b}\right)^2\$\$

**With dead time correction (effective plates):** \$\$N\_{eff} = 5.54
\left(\frac{t_R - t_0}{W\_{0.5h}}\right)^2\$\$

**Typical SEC Performance:**

- Analytical SEC columns: 10,000-40,000 plates/meter

- Preparative columns: 5,000-15,000 plates/meter

- UHPLC SEC: 50,000+ plates/meter

## See also

Other sec-qc:
[`measure_sec_asymmetry()`](https://jameshwade.github.io/measure-sec/reference/measure_sec_asymmetry.md),
[`measure_sec_column_performance()`](https://jameshwade.github.io/measure-sec/reference/measure_sec_column_performance.md),
[`measure_sec_recovery()`](https://jameshwade.github.io/measure-sec/reference/measure_sec_recovery.md),
[`measure_sec_resolution()`](https://jameshwade.github.io/measure-sec/reference/measure_sec_resolution.md),
[`measure_sec_suitability()`](https://jameshwade.github.io/measure-sec/reference/measure_sec_suitability.md)

## Examples

``` r
# Calculate plate count for a monomer peak
measure_sec_plate_count(
  retention = 9.5,
  width = 0.25,
  width_type = "half_height"
)
#> [1] 7999.76

# With dead time for effective plates
measure_sec_plate_count(
  retention = 9.5,
  width = 0.25,
  dead_time = 3.0
)
#> [1] 3745.04
```
