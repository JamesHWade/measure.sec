# RI Detector Processing for SEC

`step_sec_ri()` creates a *specification* of a recipe step that
processes refractive index (RI) detector signals for SEC analysis,
including application of dn/dc (refractive index increment) values.

## Usage

``` r
step_sec_ri(
  recipe,
  measures = NULL,
  dn_dc = NULL,
  dn_dc_column = NULL,
  instrument_constant = 1,
  output_col = NULL,
  role = NA,
  trained = FALSE,
  skip = FALSE,
  id = recipes::rand_id("sec_ri"),
  concentration_units = "g/mL"
)
```

## Arguments

- recipe:

  A recipe object.

- measures:

  Character vector of RI detector columns to process. If `NULL`, the
  step looks for columns containing "ri" in the name. Multiple columns
  can be processed in place; a named output requires exactly one.

- dn_dc:

  Refractive index increment (mL/g). Can be:

  - A single finite, nonzero numeric value applied to all samples

  - `NULL` to skip dn/dc normalization (signal remains in detector
    units)

- dn_dc_column:

  Character name of a column containing sample-specific finite, nonzero
  dn/dc values. Values may be positive or negative. Overrides `dn_dc` if
  provided.

- instrument_constant:

  RI detector instrument constant. Default is 1.0 (no adjustment). This
  positive conversion gain multiplies raw detector response to obtain
  refractive-index difference.

- output_col:

  Name for the output column. `NULL` modifies the selected RI measure in
  place. A different name creates a converted column and keeps the
  source measure unchanged.

- role:

  Role for generated columns.

- trained:

  Logical indicating if the step has been trained.

- skip:

  Logical. Should the step be skipped when baking?

- id:

  Unique step identifier.

- concentration_units:

  Units for the signal after dn/dc normalization: g/mL (default, the
  natural unit from RI divided by dn/dc) or mg/mL. This is ignored when
  dn/dc normalization is skipped.

## Value

An updated recipe with the new step added.

## Details

The RI detector response is proportional to concentration and dn/dc. Let
`S_raw` be the stored detector response and `K` be the positive
conversion gain supplied as `instrument_constant`. The calibrated
refractive-index difference and concentration are related by:

\$\$\Delta n = K \times S\_{raw} = c \times (dn/dc)\$\$

Thus the step calculates concentration in g/mL as
`S_raw * instrument_constant / dn_dc`, when dn/dc is in mL/g. A negative
dn/dc is physically valid for some solute/solvent combinations; the
[Brookhaven BI-DNDC
manual](https://www.brookhaveninstruments.com/wp-content/uploads/2023/08/BI-DNDC-Manual-Static-Mode.pdf)
describes negative values. Wyatt notes that dn/dc depends on solute,
solvent, wavelength, and temperature, so supply a value measured or
qualified for the relevant conditions ([application
note](https://www.wyatt.com/library/application-notes/determination-of-dn-dc.html)).

Set `concentration_units` to `"mg/mL"` to scale the result by 1000. The
default preserves the existing numeric behavior and can be read directly
by scattering steps that use a tagged concentration curve. With
`dn_dc = NULL` and no `dn_dc_column`, no concentration normalization is
performed. Using the default `instrument_constant = 1.0` leaves signal
values unchanged and adds no concentration-unit tag.

## See also

Other sec-detectors:
[`step_sec_concentration()`](https://jameshwade.github.io/measure.sec/reference/step_sec_concentration.md),
[`step_sec_dad()`](https://jameshwade.github.io/measure.sec/reference/step_sec_dad.md),
[`step_sec_dls()`](https://jameshwade.github.io/measure.sec/reference/step_sec_dls.md),
[`step_sec_intrinsic_visc()`](https://jameshwade.github.io/measure.sec/reference/step_sec_intrinsic_visc.md),
[`step_sec_lals()`](https://jameshwade.github.io/measure.sec/reference/step_sec_lals.md),
[`step_sec_mals()`](https://jameshwade.github.io/measure.sec/reference/step_sec_mals.md),
[`step_sec_rals()`](https://jameshwade.github.io/measure.sec/reference/step_sec_rals.md),
[`step_sec_uv()`](https://jameshwade.github.io/measure.sec/reference/step_sec_uv.md),
[`step_sec_viscometer()`](https://jameshwade.github.io/measure.sec/reference/step_sec_viscometer.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(recipes)
library(measure)

# Apply fixed dn/dc value
rec <- recipe(~., data = sec_triple_detect) |>
  step_measure_input_long(ri_signal, location = vars(elution_time), col_name = "ri") |>
  step_sec_ri(dn_dc = 0.185) |>
  prep()

# Use sample-specific dn/dc from a column
rec <- recipe(~., data = sec_triple_detect) |>
  step_measure_input_long(ri_signal, location = vars(elution_time), col_name = "ri") |>
  step_sec_ri(dn_dc_column = "dn_dc") |>
  prep()
} # }
```
