# Calculate Mass Recovery

Calculates the mass recovery (percent of injected mass detected) for SEC
analysis.

## Usage

``` r
measure_sec_recovery(detected_mass, injected_mass, units = "mg")
```

## Arguments

- detected_mass:

  Mass detected by integration of the chromatogram.

- injected_mass:

  Mass injected onto the column.

- units:

  Units for mass values. Both must be in the same units.

## Value

Numeric recovery percentage (0-100+).

## Details

Mass recovery verifies that the analytical system is detecting all of
the injected sample:

\$\$\\ Recovery = \frac{m\_{detected}}{m\_{injected}} \times 100\$\$

**Interpretation:**

- 95-105%: Excellent recovery (typical acceptance)

- 90-95% or 105-110%: Acceptable (investigate if persistent)

- \< 90%: Low recovery - possible column adsorption, precipitation

- \> 110%: High recovery - calibration issue, interference

**Common Causes of Low Recovery:**

- Sample adsorption to column packing

- Sample precipitation or aggregation on-column

- Detector calibration drift

- Integration baseline errors

- Sample degradation

## See also

Other sec-qc: [`measure_sec_asymmetry()`](measure_sec_asymmetry.md),
[`measure_sec_plate_count()`](measure_sec_plate_count.md),
[`measure_sec_resolution()`](measure_sec_resolution.md),
[`measure_sec_suitability()`](measure_sec_suitability.md)

## Examples

``` r
# Calculate recovery
measure_sec_recovery(
  detected_mass = 0.195,
  injected_mass = 0.200
)
#> [1] 97.5
# Returns 97.5%
```
