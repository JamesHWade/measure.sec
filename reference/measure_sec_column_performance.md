# Column Performance Metrics for SEC

Calculates comprehensive column performance metrics including HETP,
separation range, and efficiency parameters.

## Usage

``` r
measure_sec_column_performance(
  calibration_data,
  column_length = 30,
  column_diameter = 7.8,
  flow_rate = NULL,
  particle_size = NULL,
  dead_volume = NULL
)
```

## Arguments

- calibration_data:

  A data frame containing calibration data with columns:

  - `retention` or `retention_time`: Retention time/volume

  - `mw` or `molecular_weight`: Molecular weight

  - `width`: Peak width (optional, for plate calculation)

- column_length:

  Column length in cm. Default is 30 cm.

- column_diameter:

  Column inner diameter in mm. Default is 7.8 mm.

- flow_rate:

  Flow rate in mL/min (optional, for linear velocity).

- particle_size:

  Particle size in micrometers (optional, for reduced HETP).

- dead_volume:

  Column dead volume in mL (optional).

## Value

A list of class `sec_column_performance` containing:

- separation_range:

  MW range (exclusion limit to total permeation)

- selectivity:

  Slope of log(MW) vs retention (mL^-1)

- hetp:

  Height equivalent to theoretical plate (mm)

- plates_per_meter:

  Theoretical plates per meter

- reduced_hetp:

  HETP/particle_size (if particle_size provided)

- peak_capacity:

  Estimated peak capacity in separation range

- resolution_factor:

  Resolution per decade of MW

## Details

**Key SEC Column Performance Metrics:**

**HETP (Height Equivalent to Theoretical Plate):** \$\$HETP =
\frac{L}{N}\$\$

**Reduced HETP (h):** \$\$h = \frac{HETP}{d_p}\$\$ Optimal h ~ 2-3 for
well-packed columns.

**Selectivity (D):** \$\$D = \frac{d \log M}{dV_R}\$\$ Higher D means
better MW resolution.

**Peak Capacity:** \$\$n_c = 1 + \frac{\sqrt{N}}{4}
\ln\left(\frac{V\_{exclusion}}{V\_{total}}\right)\$\$

**Typical Performance Guidelines:**

- HETP: \< 50 um for analytical columns

- Reduced HETP: 2-5 for good columns

- Plates/meter: \> 20,000 for HPLC SEC

## See also

Other sec-qc:
[`measure_sec_asymmetry()`](https://jameshwade.github.io/measure-sec/reference/measure_sec_asymmetry.md),
[`measure_sec_plate_count()`](https://jameshwade.github.io/measure-sec/reference/measure_sec_plate_count.md),
[`measure_sec_recovery()`](https://jameshwade.github.io/measure-sec/reference/measure_sec_recovery.md),
[`measure_sec_resolution()`](https://jameshwade.github.io/measure-sec/reference/measure_sec_resolution.md),
[`measure_sec_suitability()`](https://jameshwade.github.io/measure-sec/reference/measure_sec_suitability.md)

## Examples

``` r
# From calibration standards
cal_data <- data.frame(
  retention = c(5.2, 6.1, 7.0, 8.2, 9.5, 10.8),
  mw = c(1200000, 400000, 100000, 30000, 5000, 580),
  width = c(0.4, 0.35, 0.30, 0.28, 0.25, 0.30)
)

perf <- measure_sec_column_performance(
  cal_data,
  column_length = 30,
  particle_size = 5
)
print(perf)
#> SEC Column Performance
#> ================================================== 
#> 
#> Separation Range:
#>   Exclusion limit: 1200000 Da
#>   Total permeation: 580 Da
#>   Log MW range: 3.32 decades
#> 
#> Calibration:
#>   Selectivity: 0.5801 log(MW)/unit
#>   R-squared: 0.9956
#> 
#> Column Efficiency:
#>   HETP: 0.070 mm (70.4 um)
#>   Plates/meter: 14203
#>   Reduced HETP (h): 14.08
#>   Average plates (N): 4261
#> 
#> Resolution:
#>   Peak capacity: 12.9
#>   Resolution/decade: 15.99
#> 
#> Column: 300 x 7.8 mm
#> Standards used: 6
```
