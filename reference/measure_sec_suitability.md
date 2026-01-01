# System Suitability Test for SEC

Performs comprehensive system suitability testing for SEC analysis,
evaluating resolution, plate count, asymmetry, and other quality
metrics.

## Usage

``` r
measure_sec_suitability(
  data = NULL,
  peaks,
  reference_peaks = NULL,
  injected_mass = NULL,
  criteria = NULL,
  column_length = NULL
)
```

## Arguments

- data:

  A data frame or tibble containing chromatogram data with at minimum
  retention times and peak parameters.

- peaks:

  A data frame with peak information. Must contain columns:

  - `retention`: Peak retention time

  - `width`: Peak width (at half height unless specified)

  - `area`: Peak area (for recovery calculation)

  - `height`: Peak height (optional, for asymmetry from raw data

- reference_peaks:

  Character vector of peak names to use for resolution calculation
  (e.g., c("dimer", "monomer")).

- injected_mass:

  Injected mass for recovery calculation (optional).

- criteria:

  A list of acceptance criteria. Default uses common biopharmaceutical
  criteria.

- column_length:

  Column length in cm (for plates per meter calculation).

## Value

A list of class `sec_suitability` containing:

- results:

  Data frame of calculated metrics and pass/fail status

- passed:

  Logical indicating if all criteria passed

- summary:

  Character summary of results

- criteria:

  Criteria used for evaluation

## Details

System suitability testing (SST) verifies that the chromatographic
system is performing adequately before, during, and after sample
analysis.

**Standard SEC SST Parameters:**

- Resolution: Rs \>= 1.5 between critical pair

- Plate count: N \>= specified minimum

- Tailing factor: 0.8 \<= Tf \<= 1.5

- Mass recovery: 95-105%

- Retention time RSD: \<= 1.0%

- Peak area RSD: \<= 2.0%

**Regulatory References:**

- USP \<621\> Chromatography

- ICH Q2(R1) Validation

- EP 2.2.46 Chromatographic Separation Techniques

## See also

Other sec-qc: [`measure_sec_asymmetry()`](measure_sec_asymmetry.md),
[`measure_sec_plate_count()`](measure_sec_plate_count.md),
[`measure_sec_recovery()`](measure_sec_recovery.md),
[`measure_sec_resolution()`](measure_sec_resolution.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Define peaks from integration results
peaks <- data.frame(
  name = c("aggregate", "dimer", "monomer", "fragment"),
  retention = c(7.2, 8.5, 9.8, 11.5),
  width = c(0.3, 0.25, 0.28, 0.35),
  area = c(2.1, 5.3, 89.2, 3.4)
)

# Run system suitability
sst <- measure_sec_suitability(
  peaks = peaks,
  reference_peaks = c("dimer", "monomer"),
  injected_mass = 0.200,
  column_length = 30
)

print(sst)
} # }
```
