# Tutorial: Triple Detection from Raw Data

This tutorial uses the bundled raw multi-detector SEC example to find
detector offsets, align a delay standard, and process sample signals. It
also shows what the single-angle MALS step does and which calibration
inputs an absolute estimate needs.

## Prerequisites

- Completed [Getting
  Started](https://jameshwade.github.io/measure.sec/articles/getting-started.md)
  tutorial
- Familiarity with conventional SEC calibration
- Basic understanding of light scattering (optional but helpful)

## Why triple detection?

Conventional calibration maps elution against molecular-weight
standards. Triple detection adds scattering alongside
concentration-related RI and UV signals. Absolute molecular-weight
estimates require a suitable instrument calibration constant and
concentration with declared units. The “triple” refers to:

1.  **RI detector** → concentration
2.  **UV detector** → chromophore content (optional)
3.  MALS, for a single-angle estimate; absolute values require
    calibration

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
library(tidyr)
library(ggplot2)
```

## The raw data challenge

The bundled raw example records detector signals at different elution
coordinates because the detector cells are separated in the flow path.

``` r

# Load raw multi-detector data
data(sec_raw_multidetector)

# View samples
sec_raw_multidetector |>
  distinct(sample_id, description, mw, dn_dc) |>
  print()
#> # A tibble: 4 × 4
#>   sample_id   description                           mw dn_dc
#>   <chr>       <chr>                              <dbl> <dbl>
#> 1 PS-DelayStd Narrow PS for delay determination 100000 0.185
#> 2 Sample-1    PS sample with UV absorption       75000 0.185
#> 3 Sample-2    PMMA sample (weak UV)             150000 0.085
#> 4 Sample-3    Copolymer sample                   45000 0.15
```

The delay standard shows the offset between the three channels:

``` r

# Extract delay standard
delay_std <- sec_raw_multidetector |>
  filter(sample_id == "PS-DelayStd")

# Normalize each signal for comparison
delay_plot <- delay_std |>
  mutate(
    ri_norm = ri_mv / max(ri_mv),
    uv_norm = (uv_au * 1000) / max(uv_au * 1000),  # Scale AU to visible range
    mals_norm = mals_mv / max(mals_mv)
  ) |>
  select(time_min, RI = ri_norm, UV = uv_norm, MALS = mals_norm) |>
  pivot_longer(-time_min, names_to = "detector", values_to = "signal")

peak_time_for_plot <- function(signal) {
  smoothed <- stats::filter(signal, rep(1 / 5, 5), sides = 2)
  smoothed[is.na(smoothed)] <- signal[is.na(smoothed)]
  delay_std$time_min[which.max(smoothed)]
}
peak_markers <- tibble(
  detector = c("RI", "UV", "MALS"),
  peak_time = c(
    peak_time_for_plot(delay_std$ri_mv),
    peak_time_for_plot(delay_std$uv_au),
    peak_time_for_plot(delay_std$mals_mv)
  )
)

ggplot(delay_plot, aes(time_min, signal, color = detector)) +
  geom_line(linewidth = 0.8) +
  geom_vline(
    data = peak_markers,
    aes(xintercept = peak_time, colour = detector),
    linetype = "dashed", show.legend = FALSE
  ) +
  scale_color_manual(values = c("RI" = "#2E86AB", "UV" = "#A23B72", "MALS" = "#F18F01")) +
  labs(
    x = "Time (min)",
    y = "Normalized signal",
    color = "Detector"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")
```

![Line chart of independently normalized RI, UV, and MALS signals
against time in minutes for PS-DelayStd in the synthetic
\`sec_raw_multidetector\` dataset. Dashed vertical markers show the
five-point-smoothed maximum of each channel; the maxima occur at
different times and need not equal the geometric detector delays. Each
signal is scaled to its own maximum, so response heights are not
comparable.](tutorial-triple-detection_files/figure-html/plot-misaligned-1.png)

Independently normalized detector responses for PS-DelayStd in the
synthetic `sec_raw_multidetector` dataset, before correction. Dashed
markers show maxima from the five-point moving average used below; those
signal maxima are not the dataset’s declared geometric delay volumes.

The vertical markers show signal maxima, not physical detector offsets.
`sec_raw_multidetector` stores the declared geometric offsets separately
in `delay_uv_ml` and `delay_mals_ml`. A UV or scattering response can
peak at a different time from RI because the channels weight composition
and molecular weight differently.

## Step 1: compare peak-time offsets with declared delays

The example smooths each signal with a five-point moving average, finds
its maximum, then calculates target-minus-reference peak-time offsets.
These maxima are signal-dependent and do not by themselves establish the
geometric detector delay.

``` r

# Find peak time for each detector
find_peak_time <- function(time, signal) {
  # Simple peak finding: maximum after smoothing
  smoothed <- stats::filter(signal, rep(1/5, 5), sides = 2)
  smoothed[is.na(smoothed)] <- signal[is.na(smoothed)]
  time[which.max(smoothed)]
}

# Calculate peak times
ri_peak <- find_peak_time(delay_std$time_min, delay_std$ri_mv)
uv_peak <- find_peak_time(delay_std$time_min, delay_std$uv_au)
mals_peak <- find_peak_time(delay_std$time_min, delay_std$mals_mv)

cat("Peak retention times:\n")
#> Peak retention times:
cat("  RI:   ", round(ri_peak, 3), "min\n")
#>   RI:    15.552 min
cat("  UV:   ", round(uv_peak, 3), "min\n")
#>   UV:    15.132 min
cat("  MALS: ", round(mals_peak, 3), "min\n")
#>   MALS:  14.757 min

# Calculate peak-time offsets relative to RI (the reference)
peak_offset_uv_min <- uv_peak - ri_peak
peak_offset_mals_min <- mals_peak - ri_peak

cat("\nPeak-time offsets relative to RI:\n")
#> 
#> Peak-time offsets relative to RI:
cat("  UV:   ", round(peak_offset_uv_min, 4), "min\n")
#>   UV:    -0.42 min
cat("  MALS: ", round(peak_offset_mals_min, 4), "min\n")
#>   MALS:  -0.795 min

# Compare peak-time offsets with the dataset's declared geometric delays.
flow_rate <- 1.0
peak_offset_comparison <- tibble(
  detector = c("UV", "MALS"),
  peak_offset_min = c(peak_offset_uv_min, peak_offset_mals_min),
  peak_offset_equivalent_mL = c(peak_offset_uv_min, peak_offset_mals_min) *
    flow_rate,
  declared_geometric_delay_mL = c(
    unique(delay_std$delay_uv_ml),
    unique(delay_std$delay_mals_ml)
  )
)
knitr::kable(peak_offset_comparison, digits = 3)
```

| detector | peak_offset_min | peak_offset_equivalent_mL | declared_geometric_delay_mL |
|:---|---:|---:|---:|
| UV | -0.420 | -0.420 | -0.08 |
| MALS | -0.795 | -0.795 | 0.18 |

The table compares peak-time offsets with the geometric offsets declared
in the bundled dataset. They differ in this noisy example. UV and
scattering responses can weight composition and molecular weight
differently, so a global signal maximum is not by itself a
geometric-delay estimate. The recipe below uses the declared
delay-volume fields; use values established for the instrument method in
an analysis.

## Step 2: apply delay correction

Now we apply the dataset’s declared geometric delay volumes to align the
channels to the RI reference. The stored time coordinates are rounded,
with adjacent widths from 0.0016 to 0.0017 minutes, so use
`method = "interpolate"` for this irregular grid. A positive delay means
a later target and correction samples at a larger target coordinate,
shifting its corrected signal earlier on the fixed output axis. A
negative delay samples at a smaller target coordinate and shifts the
corrected signal later.

``` r

# Process delay standard with correction
# Use the dataset's declared geometric delay volumes and minute coordinate axis.
rec_aligned <- recipe(
  ri_mv + uv_au + mals_mv + time_min + dn_dc ~ sample_id,
  data = delay_std
) |>
  update_role(sample_id, new_role = "id") |>
  # Convert each detector to measure format
  step_measure_input_long(ri_mv, location = vars(time_min), col_name = "ri") |>
  step_measure_input_long(uv_au, location = vars(time_min), col_name = "uv") |>
  step_measure_input_long(mals_mv, location = vars(time_min), col_name = "mals") |>
  # Apply baseline correction
  step_sec_baseline(measures = c("ri", "uv", "mals")) |>
  # Correct inter-detector delays
  # A positive delay means the detector sees the sample after the reference.
  step_sec_detector_delay(
    reference = "ri",
    targets = c("uv", "mals"),
    delay_volumes = c(
      uv = unique(delay_std$delay_uv_ml),
      mals = unique(delay_std$delay_mals_ml)
    ),
    flow_rate = flow_rate,
    elution_units = "min",
    method = "interpolate"
  )

aligned <- prep(rec_aligned) |> bake(new_data = NULL)
```

``` r

# Extract and normalize aligned signals
ri_chrom <- aligned$ri[[1]]
uv_chrom <- aligned$uv[[1]]
mals_chrom <- aligned$mals[[1]]

# Normalize each signal
ri_norm <- ri_chrom$value / max(ri_chrom$value, na.rm = TRUE)
uv_norm <- uv_chrom$value / max(abs(uv_chrom$value), na.rm = TRUE)
mals_norm <- mals_chrom$value / max(mals_chrom$value, na.rm = TRUE)

# Build plot data
aligned_plot <- tibble(
  location = rep(ri_chrom$location, 3),
  signal = c(ri_norm, uv_norm, mals_norm),
  detector = rep(c("RI", "UV", "MALS"), each = length(ri_norm))
)

ggplot(aligned_plot, aes(location, signal, color = detector)) +
  geom_line(linewidth = 0.8, na.rm = TRUE) +
  scale_color_manual(values = c("RI" = "#2E86AB", "UV" = "#A23B72", "MALS" = "#F18F01")) +
  labs(
    x = "Time (min)",
    y = "Normalized signal",
    color = "Detector"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")
```

![Line chart of independently normalized RI, UV, and MALS signals
against time in minutes for PS-DelayStd in the synthetic
\`sec_raw_multidetector\` dataset, after interpolation using the
declared UV and MALS delay volumes. The response maxima do not all
coincide after correction. Missing values at unsupported edges leave
gaps rather than extending the
curves.](tutorial-triple-detection_files/figure-html/plot-aligned-1.png)

Independently normalized PS-DelayStd responses after applying the
declared geometric delay volumes from the bundled dataset. Values
outside interpolation support remain missing; signal maxima can remain
offset because the detector responses differ.

The recipe applies the declared path delays, so response maxima need not
coincide. This plot shows the corrected coordinate relationship and
missing edge support; peak alignment alone is not a validation of the
molecular-weight calculation.

## Step 3: process sample signals

The recipe below applies the same declared geometric delay volumes to
Sample-1, then processes its RI, UV, and MALS signals.

``` r

# Select a sample
sample_1 <- sec_raw_multidetector |>
  filter(sample_id == "Sample-1")

# Full triple detection recipe
rec_triple <- recipe(
  ri_mv + uv_au + mals_mv + time_min + dn_dc ~ sample_id,
  data = sample_1
) |>
  update_role(sample_id, new_role = "id") |>
  # Convert to measure format
  step_measure_input_long(ri_mv, location = vars(time_min), col_name = "ri") |>
  step_measure_input_long(uv_au, location = vars(time_min), col_name = "uv") |>
  step_measure_input_long(mals_mv, location = vars(time_min), col_name = "mals") |>
  # Baseline correction
  step_sec_baseline(measures = c("ri", "uv", "mals")) |>
  # Delay correction
  step_sec_detector_delay(
    reference = "ri",
    targets = c("uv", "mals"),
    delay_volumes = c(
      uv = unique(delay_std$delay_uv_ml),
      mals = unique(delay_std$delay_mals_ml)
    ),
    flow_rate = flow_rate,
    elution_units = "min",
    method = "interpolate"
  ) |>
  # Process RI for concentration
  step_sec_ri(measures = "ri", dn_dc_column = "dn_dc")

result_sample <- prep(rec_triple) |> bake(new_data = NULL)

# View the processed data
result_sample |>
  select(sample_id, ri, uv, mals, dn_dc)
#> # A tibble: 1 × 5
#>   sample_id           ri           uv         mals dn_dc
#>   <chr>           <meas>       <meas>       <meas> <dbl>
#> 1 Sample-1  [12,001 × 2] [12,001 × 2] [12,001 × 2] 0.185
```

## Step 4: interpret a MALS estimate

A MALS estimate combines the scattering response with concentration and
optical inputs. The basic relationship is:

``` math
M_w = \frac{R_\theta}{K^* \cdot c \cdot P(\theta)}
```

Where: - $`R_\theta`$ = excess Rayleigh ratio (from MALS signal) -
$`K^*`$ = optical constant (depends on dn/dc, wavelength) - $`c`$ =
concentration (from RI) - $`P(\theta)`$ = form factor (angular
dependence)

[`step_sec_mals()`](https://jameshwade.github.io/measure.sec/reference/step_sec_mals.md)
currently uses one angle and assumes the form factor is approximately
one. It requires a concentration curve with a known unit and an
instrument calibration constant for absolute estimates. Without
calibration, the output is relative and uncalibrated; the step does not
calculate Rg. The [MALS
guide](https://jameshwade.github.io/measure.sec/articles/mals-detection.md)
plots the current API output and documents these limits.

## Complete multi-sample analysis

The code below records each sample’s RI peak position and its supplied
molecular-weight metadata. It does not calculate or compare a
molecular-weight estimate.

``` r

analyze_sample <- function(sample_data) {
  sample_id <- unique(sample_data$sample_id)
  known_mw <- unique(sample_data$mw)

  rec <- recipe(
    ri_mv + uv_au + mals_mv + time_min + dn_dc ~ sample_id,
    data = sample_data
  ) |>
    update_role(sample_id, new_role = "id") |>
    step_measure_input_long(ri_mv, location = vars(time_min), col_name = "ri") |>
    step_measure_input_long(mals_mv, location = vars(time_min), col_name = "mals") |>
    step_sec_baseline(measures = c("ri", "mals")) |>
    step_sec_detector_delay(
      reference = "ri",
      targets = "mals",
      delay_volumes = c(mals = unique(sample_data$delay_mals_ml)),
      flow_rate = flow_rate,
      elution_units = "min",
      method = "interpolate"
    ) |>
    step_sec_ri(measures = "ri", dn_dc_column = "dn_dc")

  result <- prep(rec) |> bake(new_data = NULL)

  # Extract peak info for summary
  ri_chrom <- result$ri[[1]]
  peak_idx <- which.max(ri_chrom$value)
  peak_time <- ri_chrom$location[peak_idx]

  tibble(
    sample_id = sample_id,
    known_mw = known_mw,
    peak_time = peak_time,
    dn_dc = unique(sample_data$dn_dc)
  )
}

# Analyze all samples
all_results <- sec_raw_multidetector |>
  group_by(sample_id) |>
  group_split() |>
  purrr::map(analyze_sample) |>
  bind_rows()

print(all_results)
#> # A tibble: 4 × 4
#>   sample_id   known_mw peak_time dn_dc
#>   <chr>          <dbl>     <dbl> <dbl>
#> 1 PS-DelayStd   100000      15.6 0.185
#> 2 Sample-1       75000      14.6 0.185
#> 3 Sample-2      150000      14.1 0.085
#> 4 Sample-3       45000      16.0 0.15
```

## Overlay processed chromatograms

``` r

# Process all samples for plotting
process_for_plot <- function(sample_data) {
  sid <- unique(sample_data$sample_id)
  known_mw <- unique(sample_data$mw)

  rec <- recipe(
    ri_mv + time_min + dn_dc ~ sample_id,
    data = sample_data
  ) |>
    update_role(sample_id, new_role = "id") |>
    step_measure_input_long(ri_mv, location = vars(time_min), col_name = "ri") |>
    step_sec_baseline(measures = "ri")

  result <- prep(rec) |> bake(new_data = NULL)
  chrom <- result$ri[[1]]

  tibble(
    location = chrom$location,
    value = chrom$value,
    sample_id = sid,
    mw = known_mw,
    label = paste0(sid, "\n(", scales::comma(known_mw), " Da)")
  )
}

plot_data <- sec_raw_multidetector |>
  group_by(sample_id) |>
  group_split() |>
  purrr::map(process_for_plot) |>
  bind_rows()

ggplot(plot_data, aes(location, value, color = label)) +
  geom_line(linewidth = 0.8) +
  scale_color_viridis_d(option = "turbo") +
  labs(
    x = "Elution time (min)",
    y = "Baseline-corrected RI signal",
    color = "Sample"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "right")
```

![Line chart of baseline-corrected RI signal versus time in minutes for
the samples in the synthetic \`sec_raw_multidetector\` dataset. The
legend identifies each sample and its supplied molecular-weight
metadata. The higher-molecular-weight examples appear earlier on the
elution axis; this plot shows RI traces, not an independent
molecular-weight
estimate.](tutorial-triple-detection_files/figure-html/overlay-plot-1.png)

Baseline-corrected RI traces from the bundled raw example samples.
Labels report the dataset’s supplied molecular-weight metadata; the
curves are not molecular-weight calculations.

## When to use triple detection

| Scenario | Consider |
|----|----|
| Same polymer type as available standards | Conventional calibration |
| A standard-relative result is insufficient | Light scattering with a validated calibration |
| Intrinsic viscosity is needed | Add a viscometer signal and concentration calculation |
| Detector channels are offset | Measure and correct the delays before slice-level calculations |

## Troubleshooting

| Problem | Possible Cause | Solution |
|----|----|----|
| MALS signal very noisy | Low concentration or small molecules | Increase injection, use LALS/RALS |
| MW values too high | Aggregates or dust in sample | Filter sample, check for aggregates |
| MW values too low | Incorrect dn/dc | Measure dn/dc or use literature value |
| Peaks still misaligned | Wrong delay values | Re-run delay determination |
| Missing estimates near peak edges | Signal below the configured threshold | Check the threshold and select integration support based on the method |

## Summary

The example applies the dataset’s declared UV and MALS delay volumes,
converts them to time using the stated flow rate, and interpolates onto
the fixed time axis. The response maxima are diagnostics; they do not
estimate those path delays. Re-measure detector delays after changing
the flow path. The sample recipe processes RI, UV, and MALS signals; an
absolute MALS molecular-weight estimate also requires a calibrated
scattering channel and concentration with a declared unit.

## Next steps

| Topic | Resource |
|----|----|
| MALS theory and angular analysis | [MALS Detection Guide](https://jameshwade.github.io/measure.sec/articles/mals-detection.md) |
| Single-angle alternatives | [LALS/RALS Detection](https://jameshwade.github.io/measure.sec/articles/lals-rals-detection.md) |
| Adding viscometer for \[η\] | [Multi-Detector SEC](https://jameshwade.github.io/measure.sec/articles/triple-detection.md) |
| Conventional calibration | [Conventional Calibration Tutorial](https://jameshwade.github.io/measure.sec/articles/tutorial-conventional-calibration.md) |

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
#> [1] ggplot2_4.0.4          tidyr_1.3.2            measure.sec_0.0.0.9000
#> [4] measure_0.0.1.9003     recipes_1.4.0          dplyr_1.2.1           
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
#> [25] prodlim_2026.03.11  pillar_1.11.1       pkgdown_2.2.1      
#> [28] jquerylib_0.1.4     MASS_7.3-66         cachem_1.1.0       
#> [31] gower_1.0.2         rpart_4.1.27        parallelly_1.50.0  
#> [34] lava_1.11           tidyselect_1.2.1    digest_0.6.39      
#> [37] future_1.82.0       purrr_1.2.2         listenv_1.1.0      
#> [40] labeling_0.4.3      splines_4.6.2       fastmap_1.2.0      
#> [43] grid_4.6.2          cli_3.6.6           magrittr_2.0.5     
#> [46] utf8_1.2.6          survival_3.8-12     future.apply_1.21.0
#> [49] withr_3.0.3         scales_1.4.0        lubridate_1.9.5    
#> [52] timechange_0.4.0    rmarkdown_2.32      globals_0.20.0     
#> [55] otel_0.3.0          nnet_7.3-21         timeDate_4052.112  
#> [58] ragg_1.5.2          evaluate_1.0.5      knitr_1.52         
#> [61] hardhat_1.4.3       viridisLite_0.4.3   rlang_1.3.0        
#> [64] Rcpp_1.1.2          glue_1.8.1          ipred_0.9-16       
#> [67] jsonlite_2.0.0      R6_2.6.1            systemfonts_1.3.3  
#> [70] fs_2.1.0
```
