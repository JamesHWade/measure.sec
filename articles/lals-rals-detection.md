# Single-Angle Light Scattering: LALS and RALS

## Overview

[`step_sec_lals()`](https://jameshwade.github.io/measure.sec/reference/step_sec_lals.md)
and
[`step_sec_rals()`](https://jameshwade.github.io/measure.sec/reference/step_sec_rals.md)
combine a single-angle scattering signal with concentration to calculate
a molecular-weight estimate. An instrument calibration constant is
required for an absolute value. Neither step estimates radius of
gyration (Rg). This guide describes the angle assumptions and the inputs
these steps use.

## Setup

``` r

library(measure)
#> Loading required package: recipes
#> Loading required package: dplyr
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
#> 
#> Attaching package: 'recipes'
#> The following object is masked from 'package:stats':
#> 
#>     step
library(measure.sec)
library(recipes)
library(dplyr)
library(ggplot2)
```

## What a single angle tells you

A particle’s scattering response depends on the angle of observation. At
low angles the particle form factor approaches one. The steps take the
signal from one detector angle and calculate molecular weight using

``` math
M_w = \frac{S(\theta) \cdot C_{cal}}{K \cdot c_{g/mL}}.
```

Here $`S(\theta)`$ is the detector signal at the observed angle,
$`C_{cal}`$ is its calibration constant, $`K`$ is the optical constant,
and $`c_{g/mL}`$ is the concentration at the same elution time. Omitting
`calibration_constant` gives a relative estimate. These single-angle
steps do not estimate radius of gyration.

## From refractive index to a molecular-weight estimate

A light-scattering calculation needs the concentration at each elution
time. Peak height can serve as a scale for comparing chromatogram
shapes, but dividing by it removes the concentration scale. Keep a
calibrated channel for the molecular-weight calculation.

Here is a five-point example with a known answer. The concentration
peaks at 2 mg/mL. We use it to construct the change in refractive index,
Δn, and a scattering signal for a sample with constant molecular weight.
The values and calibration constant are chosen to check the calculation.

``` r

time_min <- 8:12
conc_mg_ml <- c(0.4, 1.0, 2.0, 1.0, 0.4)
conc_g_ml <- conc_mg_ml / 1000
dn_dc <- 0.185 # mL/g
target_mw <- 66430 # g/mol
lambda_cm <- 670e-7
optical_constant <- 4 * pi^2 * 1.333^2 * dn_dc^2 /
  (6.022e23 * lambda_cm^4)

channel <- function(values) {
  new_measure_list(list(new_measure_tbl(time_min, values)))
}

run <- tibble(
  ri = channel(dn_dc * conc_g_ml),
  lals = channel(optical_constant * conc_g_ml * target_mw)
)
```

[`step_sec_ri()`](https://jameshwade.github.io/measure.sec/reference/step_sec_ri.md)
divides Δn by dn/dc and expresses the result in the requested unit. Here
we smooth both channels with the same three-point moving average.
Because the constructed channels are proportional, their ratio should
remain constant even where smoothing changes the peak height.

``` r

smoothed <- recipe(~ ., data = run) |>
  step_sec_ri(
    measures = "ri", dn_dc = dn_dc,
    output_col = "concentration", concentration_units = "mg/mL"
  ) |>
  step_measure_smooth_ma(measures = c("concentration", "lals"), window = 3)

estimate <- smoothed |>
  step_sec_lals(
    measures = "lals", concentration_col = "concentration",
    optical_constant = optical_constant, calibration_constant = 1
  ) |>
  prep() |>
  bake(new_data = run)

measure_value_metadata(estimate$concentration)
#> $quantity
#> [1] "mass_concentration"
#> 
#> $unit
#> [1] "mg/mL"
```

The LALS step reads mg/mL from the smoothed concentration column and
converts to g/mL for the calculation. The dashed line in the right panel
marks the 66,430 g/mol value used to construct the input.

``` r

estimated_mw <- estimate$mw_lals[[1]]$value
stopifnot(isTRUE(all.equal(estimated_mw, rep(target_mw, 5))))

curves <- bind_rows(
  tibble(time_min, panel = "Concentration (mg/mL)", series = "Input",
         value = conc_mg_ml),
  tibble(time_min, panel = "Concentration (mg/mL)", series = "Smoothed",
         value = estimate$concentration[[1]]$value),
  tibble(time_min, panel = "Molecular weight (g/mol)", series = "Calculated",
         value = estimated_mw)
)

ggplot(curves, aes(time_min, value, colour = series, linetype = series)) +
  geom_blank(
    data = tibble(panel = "Molecular weight (g/mol)", bound = c(0, 70000)),
    aes(x = 8, y = bound), inherit.aes = FALSE
  ) +
  geom_hline(
    data = tibble(panel = "Molecular weight (g/mol)", target = target_mw),
    aes(yintercept = target), linetype = "dashed", colour = "grey40"
  ) +
  geom_line(data = filter(curves, series != "Calculated"), linewidth = 0.8) +
  geom_point() +
  facet_wrap(~panel, scales = "free_y") +
  scale_colour_manual(values = c(
    Input = "#555555", Smoothed = "#0072B2", Calculated = "#A54500"
  )) +
  labs(x = "Elution time (min)", y = NULL, colour = NULL, linetype = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")
```

![Two plots cover 8 to 12 minutes. The input concentration peaks at 2
mg/mL and its three-point moving average peaks near 1.33 mg/mL. Five
molecular-weight estimates lie on a horizontal reference at 66,430
g/mol.](lals-rals-detection_files/figure-html/lals-software-example-1.png)

Moving-average smoothing lowers and broadens the concentration peak.
Smoothing the proportional scattering channel the same way leaves the
calculated molecular weight at the constructed target across all five
points.

### Comparing peak shapes

For a shape comparison, branch from the smoothed recipe and divide each
concentration curve by its maximum. The peak now has height 1. Its
values are fractions of the peak height. Passing this column as the
concentration in
[`step_sec_lals()`](https://jameshwade.github.io/measure.sec/reference/step_sec_lals.md)
gives an error, even if you supply a concentration unit.

``` r

shape <- smoothed |>
  step_measure_normalize_max(measures = "concentration") |>
  prep() |>
  bake(new_data = run)

measure_value_metadata(shape$concentration)
#> $quantity
#> [1] "normalized_signal"
#> 
#> $unit
#> [1] "1"
```

### Reusing a saved concentration curve

Older SEC recipes stored the unit in an attribute called
`sec_concentration_unit`. SEC detector steps still read it. To smooth
such a curve with a core measure step, copy the recorded unit into the
declaration first. For a saved concentration column `x`, read
`unit <- attr(x, "sec_concentration_unit", exact = TRUE)`. If the
recorded unit is `"mg/mL"` or `"g/mL"`, use
`x <- set_measure_value_metadata(x, "mass_concentration", unit)`. This
labels the existing numbers. If the unit is missing, consult the
calibration record before declaring one.

## RALS: right-angle light scattering

RALS records at 90 degrees by default. Larger particles can scatter less
light at this angle than at a low angle. Use the instrument’s RALS
channel and its calibration constant with
[`step_sec_rals()`](https://jameshwade.github.io/measure.sec/reference/step_sec_rals.md).
Supply the smoothed concentration in the same way as for LALS.

## Choosing a step

| Step | Angle | Current calculation | Absolute values |
|----|----|----|----|
| [`step_sec_lals()`](https://jameshwade.github.io/measure.sec/reference/step_sec_lals.md) | Low angle; defaults to 7 degrees | Single-angle estimate | Requires an instrument calibration constant |
| [`step_sec_rals()`](https://jameshwade.github.io/measure.sec/reference/step_sec_rals.md) | 90 degrees by default | Single-angle estimate | Requires an instrument calibration constant |
| [`step_sec_mals()`](https://jameshwade.github.io/measure.sec/reference/step_sec_mals.md) | One supplied angle | Single-angle approximation | Requires an instrument calibration constant; no multi-angle fit or Rg output |

Choose the step that matches the detector channel and calibration
available. The package examples do not establish a size range or
comparative accuracy for the instruments.

## Combining LALS/RALS with MALS

Some systems record at several scattering angles. Multiple channels can
provide cross-checks, but
[`step_sec_mals()`](https://jameshwade.github.io/measure.sec/reference/step_sec_mals.md)
currently uses one angle and does not fit a multi-angle model. The LALS
and RALS steps each calculate from their own single-angle channel.

## Calibration constants

Both
[`step_sec_lals()`](https://jameshwade.github.io/measure.sec/reference/step_sec_lals.md)
and
[`step_sec_rals()`](https://jameshwade.github.io/measure.sec/reference/step_sec_rals.md)
need a `calibration_constant` for absolute estimates. Without it,
results are relative and carry no molecular weight unit.

### Determining calibration constants

1.  **Use a well-characterized standard** (e.g., BSA, narrow PS
    standard)
2.  **Run the standard** under your experimental conditions
3.  **Calculate the constant** that yields the known MW

``` r

# Example: Calibrating with BSA (MW = 66,430 g/mol)
# 1. Measure BSA with known concentration
# 2. Convert c to g/mL, then use the same equation as step_sec_lals/rals:
#    MW = raw_signal * cal_constant / (K * c)

bsa_mw <- 66430  # g/mol (numerically equal to Da)
bsa_conc_mg_ml <- 1.0
bsa_conc_g_ml <- bsa_conc_mg_ml / 1000
dn_dc_bsa <- 0.185
solvent_ri <- 1.333
wavelength <- 670

# Calculate K as in the current step implementation. Wavelength is converted
# from nm to cm, and dn/dc is in mL/g (= cm^3/g).
lambda_cm <- wavelength * 1e-7
K <- 4 * pi^2 * solvent_ri^2 * dn_dc_bsa^2 /
  (6.022e23 * lambda_cm^4)

# For a measured raw_signal of 0.5, solve the API equation for the constant.
raw_signal <- 0.5
calibration_constant <- bsa_mw * K * bsa_conc_g_ml / raw_signal
```

## Troubleshooting

### Common issues

| Problem              | Possible Cause    | Solution                       |
|----------------------|-------------------|--------------------------------|
| MW too high          | Aggregates, dust  | Filter samples, check baseline |
| MW too low           | Poor dn/dc value  | Measure dn/dc accurately       |
| Noisy signal         | Low concentration | Increase injection mass        |
| Inconsistent results | Calibration drift | Recalibrate with standard      |

### Signal-to-noise considerations

Check detector response, concentration, baseline, and flow stability
over the support used for the calculation. The current MALS step also
uses one angle; it does not average measurements across angles.

## See also

- [MALS
  Detection](https://jameshwade.github.io/measure.sec/articles/mals-detection.md) -
  Current single-angle implementation and multi-angle limits
- [Triple
  Detection](https://jameshwade.github.io/measure.sec/articles/triple-detection.md) -
  Combined RI + Viscometer + LS workflows
- [Getting
  Started](https://jameshwade.github.io/measure.sec/articles/getting-started.md) -
  Basic SEC workflow and concepts
- [Calibration
  Management](https://jameshwade.github.io/measure.sec/articles/calibration-management.md) -
  Save and reuse calibrations

## Session info

``` r

sessionInfo()
#> R version 4.6.2 (2026-10-29)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 26.04.1 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.32.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] ggplot2_4.0.4          measure.sec_0.0.0.9000 measure_0.0.1.9004    
#> [4] recipes_1.4.0          dplyr_1.2.1           
#> 
#> loaded via a namespace (and not attached):
#>  [1] gtable_0.3.6        xfun_0.62           bslib_0.12.0       
#>  [4] lattice_0.23-1      vctrs_0.7.3         tools_4.6.2        
#>  [7] generics_0.1.4      parallel_4.6.2      tibble_3.3.1       
#> [10] pkgconfig_2.0.3     Matrix_1.7-6        data.table_1.18.6.1
#> [13] RColorBrewer_1.1-3  S7_0.2.2            desc_1.4.3         
#> [16] lifecycle_1.0.5     compiler_4.6.2      farver_2.1.2       
#> [19] textshaping_1.0.5   codetools_0.2-20    htmltools_0.5.9    
#> [22] class_7.3-24        sass_0.4.10         yaml_2.3.12        
#> [25] prodlim_2026.03.11  tidyr_1.3.2         pillar_1.11.1      
#> [28] pkgdown_2.2.1       jquerylib_0.1.4     MASS_7.3-66        
#> [31] cachem_1.1.0        gower_1.0.2         rpart_4.1.27       
#> [34] parallelly_1.50.0   lava_1.11           tidyselect_1.2.1   
#> [37] digest_0.6.39       future_1.82.0       purrr_1.2.2        
#> [40] listenv_1.1.0       labeling_0.4.3      splines_4.6.2      
#> [43] fastmap_1.2.0       grid_4.6.2          cli_3.6.6          
#> [46] magrittr_2.0.5      survival_3.8-12     future.apply_1.21.0
#> [49] withr_3.0.3         scales_1.4.0        lubridate_1.9.5    
#> [52] timechange_0.4.0    rmarkdown_2.32      globals_0.20.0     
#> [55] otel_0.3.0          nnet_7.3-21         timeDate_4052.112  
#> [58] ragg_1.5.2          evaluate_1.0.5      knitr_1.52         
#> [61] hardhat_1.4.3       rlang_1.3.0         Rcpp_1.1.2         
#> [64] glue_1.8.1          ipred_0.9-16        jsonlite_2.0.0     
#> [67] R6_2.6.1            systemfonts_1.3.3   fs_2.1.0
```
