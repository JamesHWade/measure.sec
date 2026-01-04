# Save SEC Calibration Parameters

Extracts calibration parameters from a prepped recipe and saves them to
a file for later reuse. This allows calibrations to be established once
and applied to future analyses without re-fitting.

## Usage

``` r
save_sec_calibration(
  prepped_recipe,
  file,
  step_number = NULL,
  metadata = NULL,
  overwrite = FALSE
)
```

## Arguments

- prepped_recipe:

  A prepped recipe containing a trained calibration step (e.g.,
  `step_sec_conventional_cal`).

- file:

  Path to save the calibration file. File extension determines format:
  `.rds` for RDS format (default), `.yaml` or `.yml` for YAML format.

- step_number:

  Integer. Which step to extract calibration from. If NULL (default),
  finds the first calibration step automatically.

- metadata:

  Optional named list of additional metadata to include (e.g.,
  `instrument`, `column`, `analyst`, `date`).

- overwrite:

  Logical. If TRUE, overwrite existing file. Default FALSE.

## Value

Invisibly returns the calibration object that was saved.

## Details

The saved calibration includes:

- Calibration fit coefficients and polynomial degree

- Calibration range (valid elution range)

- Fit diagnostics (R², RMSE, per-standard results)

- Original standards data

- Settings (fit_type, extrapolation, log_output)

- Timestamp and version information

- Optional user-provided metadata

RDS format preserves all R objects exactly and is recommended for most
uses. YAML format is human-readable and useful for documentation or
version control.

## See also

[`load_sec_calibration()`](https://jameshwade.github.io/measure-sec/reference/load_sec_calibration.md),
[`step_sec_conventional_cal()`](https://jameshwade.github.io/measure-sec/reference/step_sec_conventional_cal.md)

Other sec-calibration:
[`load_sec_calibration()`](https://jameshwade.github.io/measure-sec/reference/load_sec_calibration.md),
[`sec_calibration_standards`](https://jameshwade.github.io/measure-sec/reference/sec_calibration_standards.md),
[`sec_pmma_standards`](https://jameshwade.github.io/measure-sec/reference/sec_pmma_standards.md),
[`sec_ps_standards`](https://jameshwade.github.io/measure-sec/reference/sec_ps_standards.md),
[`step_sec_broad_standard()`](https://jameshwade.github.io/measure-sec/reference/step_sec_broad_standard.md),
[`step_sec_conventional_cal()`](https://jameshwade.github.io/measure-sec/reference/step_sec_conventional_cal.md),
[`step_sec_universal_cal()`](https://jameshwade.github.io/measure-sec/reference/step_sec_universal_cal.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Create and prep a calibration recipe
ps_standards <- data.frame(
  retention = c(12.5, 13.2, 14.1, 15.0, 16.2, 17.5),
  log_mw = c(6.0, 5.5, 5.0, 4.5, 4.0, 3.5)
)

rec <- recipe(~., data = sample_data) |>
  step_sec_conventional_cal(standards = ps_standards, fit_type = "cubic") |>
  prep()

# Save calibration for later use
save_sec_calibration(
  rec,
  "ps_calibration.rds",
  metadata = list(
    column = "PLgel 5um Mixed-C",
    instrument = "Agilent 1260",
    analyst = "JW"
  )
)
} # }
```
