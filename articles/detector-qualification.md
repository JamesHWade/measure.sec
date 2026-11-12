# Check detector units, signed RI and clipped signal

This guide checks numerical contracts using fixed synthetic traces:
converting calibrated RI signals to concentration, conserving declared
injected mass, and inspecting negative-sample clipping. It does not
establish detector calibration, recovery, a suitable baseline, or
instrument performance.

## Convert a signed RI response

`instrument_constant` is a positive conversion gain: multiplying the
recorded signal by it gives the calibrated refractive-index difference.
Dividing that difference by `dn_dc` in mL/g gives concentration in g/mL.
A negative increment can therefore turn a negative RI response into
positive concentration without taking the absolute value of the detector
trace.

The [Brookhaven BI-DNDC
manual](https://www.brookhaveninstruments.com/wp-content/uploads/2023/08/BI-DNDC-Manual-Static-Mode.pdf)
describes negative refractive-index increments. The [Wyatt application
note](https://www.wyatt.com/library/application-notes/determination-of-dn-dc.html)
explains why the appropriate value depends on the solute, solvent,
temperature and wavelength. The values below are synthetic declarations,
not recommended values for a polymer or solvent.

``` r

time_min <- c(0, 1, 2, 4)
concentration_g_ml <- c(0, 0.001, 0.002, 0)
gain <- 2
flow_ml_min <- 0.5
injected_mass_mg <- 2
input <- tibble::tibble(
  sample_id = c("positive increment", "negative increment"),
  dn_dc = c(0.2, -0.2)
)
input$ri <- new_measure_list(lapply(input$dn_dc, function(increment) {
  new_measure_tbl(time_min, concentration_g_ml * increment / gain)
}))

ri_recipe <- recipe(~., data = input) |>
  step_sec_ri(
    measures = "ri", dn_dc_column = "dn_dc",
    instrument_constant = gain, output_col = "concentration"
  )
ri_result <- bake(prep(ri_recipe), new_data = NULL)
stopifnot(identical(ri_result$ri, input$ri),
          identical(attr(ri_result$concentration, "sec_concentration_unit"),
                    "g/mL"))
for (row in seq_len(nrow(input))) {
  stopifnot(max(abs(ri_result$concentration[[row]]$value -
                    concentration_g_ml)) < 1e-14)
}

recipe_for <- function(unit) {
  ri_recipe |>
    step_sec_concentration(
      measures = "concentration", injection_mass = injected_mass_mg,
      flow_rate = flow_ml_min, concentration_units = unit
    )
}
result_mg <- bake(prep(recipe_for("mg/mL")), new_data = NULL)
result_g <- bake(prep(recipe_for("g/mL")), new_data = NULL)
```

The RI-only check runs before injected-mass normalization, which would
otherwise cancel a constant conversion error. Both rows represent the
same concentration curve. The signed increment is read from each current
baked row; it is not replaced by the training rows’ values. Finite,
nonzero increments are required. Without `dn_dc` or a `dn_dc_column`,
the step applies only the declared instrument gain and does not create a
concentration-unit tag.

## Check the declared mass and units

The elution axis here is **minutes**, and flow is **mL/min**. Convert
other axis units before using this concentration step: monotonic numeric
coordinates do not establish their physical unit.

For the declared curve, trapezoidal concentration area is 0.004 g
min/mL. Multiplication by the flow gives 0.002 g, or 2 mg. The following
independent balance checks the baked results in both requested
concentration units.

``` r

recovered_mass_mg <- function(curve, unit) {
  values <- curve$value
  if (identical(unit, "g/mL")) values <- values * 1000
  widths <- abs(diff(curve$location))
  flow_ml_min * sum(widths *
    (utils::head(values, -1) + utils::tail(values, -1)) / 2)
}
checks <- dplyr::bind_rows(lapply(
  list(result_mg, result_g), function(result) {
    unit <- attr(result$concentration, "sec_concentration_unit", exact = TRUE)
    recovered <- vapply(result$concentration, recovered_mass_mg,
                        numeric(1), unit = unit)
    tibble::tibble(
      sample_id = result$sample_id,
      unit = unit,
      recovered_mass_mg = recovered
    )
  }
))
stopifnot(max(abs(checks$recovered_mass_mg - injected_mass_mg)) < 1e-12)
stopifnot(identical(result_mg$ri, input$ri),
          identical(result_g$ri, input$ri),
          is.null(attr(result_mg$ri, "sec_concentration_unit", exact = TRUE)))
for (row in seq_len(nrow(input))) {
  stopifnot(max(abs(result_g$concentration[[row]]$value - concentration_g_ml)) < 1e-14)
  stopifnot(max(abs(result_mg$concentration[[row]]$value / 1000 -
                    result_g$concentration[[row]]$value)) < 1e-14)
}
knitr::kable(checks, digits = 8)
```

| sample_id          | unit  | recovered_mass_mg |
|:-------------------|:------|------------------:|
| positive increment | mg/mL |                 2 |
| negative increment | mg/mL |                 2 |
| positive increment | g/mL  |                 2 |
| negative increment | g/mL  |                 2 |

The mass balance checks the supplied declaration, not an independently
measured injection or column recovery. Here
`output_col = "concentration"` keeps the raw `ri` curves unchanged and
creates a separate concentration column. Its unit tag is checked
independently of the raw signal. Omit `output_col` to retain the
existing in-place behavior, or explicitly use the source column name. A
new output name requires one selected RI measure and must be unused at
that stage of the recipe.

## Inspect clipping on the actual sampling grid

Mass normalization sets negative sampled values to zero. The warning
includes a structured `diagnostics` tibble with one row for every curve
in that current bake. Capture this condition when you need to retain the
diagnostics; the trained recipe does not keep a mutable record of its
latest bake.

These three curves have identical sampled values and one negative point
each. Their first sampling interval differs, changing the quadrature
weight carried by that point.

``` r

first_interval_min <- c(0.1, 1, 9)
clipping_input <- tibble::tibble(first_interval_min = first_interval_min)
clipping_input$ri <- new_measure_list(lapply(first_interval_min, function(first) {
  new_measure_tbl(c(0, first, 10), c(-1, 1, 1))
}))
clipping_recipe <- recipe(~., data = clipping_input) |>
  step_sec_concentration(measures = "ri", injection_mass = 1, flow_rate = 1)

captured <- new.env(parent = emptyenv())
trained <- withCallingHandlers(
  prep(clipping_recipe),
  measure_sec_negative_signal_clipped = function(condition) {
    captured$diagnostics <- condition$diagnostics
    invokeRestart("muffleWarning")
  }
)
clipped_result <- bake(trained, new_data = NULL)
diagnostics <- captured$diagnostics
diagnostics$first_interval_min <- first_interval_min[diagnostics$row]
stopifnot(max(abs(diagnostics$excluded_fraction - c(0.005, 0.05, 0.45))) < 1e-14)
knitr::kable(diagnostics, digits = 4)
```

| row | negative_points | points | positive_area | negative_area | excluded_fraction | first_interval_min |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 1 | 3 | 9.95 | 0.05 | 0.005 | 0.1 |
| 2 | 1 | 3 | 9.50 | 0.50 | 0.050 | 1.0 |
| 3 | 1 | 3 | 5.50 | 4.50 | 0.450 | 9.0 |

`positive_area` is the trapezoidal area of the samplewise positive part,
and `negative_area` is the same calculation on the magnitudes of
negative samples. Their units are incoming signal units times minutes.
`excluded_fraction` is
`negative_area / (positive_area + negative_area)`: the negative share of
**absolute sampled signal area**. It is not an estimate of lost or
recovered mass, nor an exact integration of zero crossings between the
original samples.

``` r

ggplot2::ggplot(diagnostics, ggplot2::aes(
  x = factor(first_interval_min, levels = first_interval_min),
  y = 100 * excluded_fraction
)) +
  ggplot2::geom_col(fill = "#4A90A4", width = 0.65) +
  ggplot2::geom_hline(yintercept = 100 / 3, linetype = "dashed", colour = "grey45") +
  ggplot2::scale_y_continuous(limits = c(0, 50), expand = c(0, 0)) +
  ggplot2::labs(
    x = "First sampling interval (min)",
    y = "Negative share of absolute sampled signal area (%)",
    caption = "Dashed line: negative-point percentage, 33.3% for all three curves."
  ) +
  ggplot2::theme_minimal(base_size = 11)
```

![Three curves each have one negative point out of three. Their negative
shares of absolute sampled signal area are 0.5, 5 and 45 percent as the
first sampling interval grows from 0.1 to 1 to 9 minutes. The dashed
line marks the common negative-point percentage of 33.3
percent.](detector-qualification_files/figure-html/clipping-figure-1.png)

The warning keeps its point count and reports the largest area share
with its current row position. Fractions are calculated per curve rather
than pooling raw signal areas across samples. Reordering curves changes
their diagnostic row positions; associate rows with identifiers from the
same baked dataset.
[`prep()`](https://recipes.tidymodels.org/reference/prep.html) processes
training data and can emit this warning as well.

A high excluded share calls for investigation of the input and baseline;
this guide defines no universal acceptance threshold. Choose any
acceptance rule from appropriate blanks, controls and method
requirements. The clipping share is not a complete uncertainty budget.
The detector-unit metadata on these curves also does not imply that
every generic preprocessing transform preserves it; check units before
downstream detector calculations.
