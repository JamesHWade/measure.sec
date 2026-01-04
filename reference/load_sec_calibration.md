# Load SEC Calibration Parameters

Loads previously saved SEC calibration parameters for use in
[`step_sec_conventional_cal()`](https://jameshwade.github.io/measure-sec/reference/step_sec_conventional_cal.md).
This enables reusing established calibrations without refitting.

## Usage

``` r
load_sec_calibration(file)
```

## Arguments

- file:

  Path to the calibration file (`.rds` or `.yaml`/`.yml`).

## Value

A `sec_calibration` object that can be passed to the `calibration`
argument of
[`step_sec_conventional_cal()`](https://jameshwade.github.io/measure-sec/reference/step_sec_conventional_cal.md).

## Details

The loaded calibration can be used in two ways:

1.  **Direct use**: Pass to
    `step_sec_conventional_cal(calibration = cal)` to skip fitting and
    use the pre-established calibration.

2.  **Inspection**: Use [`print()`](https://rdrr.io/r/base/print.html)
    or [`summary()`](https://rdrr.io/r/base/summary.html) to view
    calibration details.

## See also

[`save_sec_calibration()`](https://jameshwade.github.io/measure-sec/reference/save_sec_calibration.md),
[`step_sec_conventional_cal()`](https://jameshwade.github.io/measure-sec/reference/step_sec_conventional_cal.md)

Other sec-calibration:
[`save_sec_calibration()`](https://jameshwade.github.io/measure-sec/reference/save_sec_calibration.md),
[`sec_calibration_standards`](https://jameshwade.github.io/measure-sec/reference/sec_calibration_standards.md),
[`sec_pmma_standards`](https://jameshwade.github.io/measure-sec/reference/sec_pmma_standards.md),
[`sec_ps_standards`](https://jameshwade.github.io/measure-sec/reference/sec_ps_standards.md),
[`step_sec_broad_standard()`](https://jameshwade.github.io/measure-sec/reference/step_sec_broad_standard.md),
[`step_sec_conventional_cal()`](https://jameshwade.github.io/measure-sec/reference/step_sec_conventional_cal.md),
[`step_sec_universal_cal()`](https://jameshwade.github.io/measure-sec/reference/step_sec_universal_cal.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Load a saved calibration
cal <- load_sec_calibration("ps_calibration.rds")
print(cal)

# Use in a recipe (no fitting needed)
rec <- recipe(~., data = new_samples) |>
  step_sec_conventional_cal(calibration = cal) |>
  prep()
} # }
```
