# Convert Detector Signal to Concentration

`step_sec_concentration()` creates a *specification* of a recipe step
that converts detector signals to absolute concentration values using
calibration factors and injection parameters.

## Usage

``` r
step_sec_concentration(
  recipe,
  measures = NULL,
  detector = c("ri", "uv", "auto"),
  injection_volume = NULL,
  injection_mass = NULL,
  sample_concentration = NULL,
  flow_rate = 1,
  concentration_units = "mg/mL",
  normalize_to_mass = TRUE,
  role = NA,
  trained = FALSE,
  skip = FALSE,
  id = recipes::rand_id("sec_concentration")
)
```

## Arguments

- recipe:

  A recipe object.

- measures:

  Character vector of detector column names to convert. If `NULL`, will
  convert all measure columns.

- detector:

  Type of detector signal being converted:

  - `"ri"`: Refractive index detector (assumes dn/dc normalized)

  - `"uv"`: UV detector (assumes extinction coefficient normalized)

  - `"auto"`: Attempt to detect from column names

- injection_volume:

  Injection volume in uL. Required for absolute concentration
  calculation.

- injection_mass:

  Injected mass in mg. Alternative to using `sample_concentration` and
  `injection_volume`.

- sample_concentration:

  Sample concentration in mg/mL. Used with `injection_volume` to
  calculate injected mass.

- flow_rate:

  Flow rate in mL/min for peak area calculations.

- concentration_units:

  Output concentration units. Default is `"mg/mL"`.

- normalize_to_mass:

  Logical. If `TRUE`, normalize the chromatogram so that the total area
  equals the injected mass. Default is `TRUE`.

- role:

  Role for generated columns.

- trained:

  Logical indicating if the step has been trained.

- skip:

  Logical. Should the step be skipped when baking?

- id:

  Unique step identifier.

## Value

An updated recipe with the new step added.

## Details

This step converts detector response (after dn/dc or extinction
coefficient normalization) to absolute concentration. The conversion
uses the known injected mass to normalize the chromatogram area:

\$\$c(t) = \frac{S(t) \times m\_{inj}}{\int S(t) \times F \times dt}\$\$

where:

- S(t) is the normalized detector signal

- m_inj is the injected mass

- F is the flow rate

**Workflow for concentration determination:**

1.  Baseline correct the chromatogram

2.  Apply detector-specific normalization (dn/dc or extinction
    coefficient)

3.  Apply this step with known injection parameters

4.  Result: concentration at each elution point

## See also

Other sec-detectors:
[`step_sec_dad()`](https://jameshwade.github.io/measure-sec/reference/step_sec_dad.md),
[`step_sec_dls()`](https://jameshwade.github.io/measure-sec/reference/step_sec_dls.md),
[`step_sec_intrinsic_visc()`](https://jameshwade.github.io/measure-sec/reference/step_sec_intrinsic_visc.md),
[`step_sec_lals()`](https://jameshwade.github.io/measure-sec/reference/step_sec_lals.md),
[`step_sec_mals()`](https://jameshwade.github.io/measure-sec/reference/step_sec_mals.md),
[`step_sec_rals()`](https://jameshwade.github.io/measure-sec/reference/step_sec_rals.md),
[`step_sec_ri()`](https://jameshwade.github.io/measure-sec/reference/step_sec_ri.md),
[`step_sec_uv()`](https://jameshwade.github.io/measure-sec/reference/step_sec_uv.md),
[`step_sec_viscometer()`](https://jameshwade.github.io/measure-sec/reference/step_sec_viscometer.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(recipes)
library(measure)

# Convert RI signal to concentration
rec <- recipe(~., data = sec_triple_detect) |>
  step_measure_input_long(ri_signal, location = vars(elution_time), col_name = "ri") |>
  step_sec_baseline() |>
  step_sec_ri(dn_dc = 0.185) |>
  step_sec_concentration(
    detector = "ri",
    injection_volume = 100,        # uL
    sample_concentration = 2.0,    # mg/mL
    flow_rate = 1.0                # mL/min
  ) |>
  prep()
} # }
```
