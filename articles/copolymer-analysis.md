# Copolymer Composition Analysis

## Overview

SEC with dual detection (UV and RI) can reveal compositional
heterogeneity in copolymers. By analyzing the UV/RI ratio across the
molecular weight distribution, you can:

- Detect compositional drift during polymerization
- Identify blocky vs random copolymer structures
- Quantify the composition at different molecular weights

This vignette covers:

1.  UV/RI ratio analysis
2.  Composition calculations
3.  Interpreting compositional heterogeneity

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

## Copolymer Detection Principles

### UV and RI Response

Different monomers have different detector responses:

| Monomer Type             | UV Response | RI Response    |
|--------------------------|-------------|----------------|
| Aromatic (styrene)       | Strong      | Moderate       |
| Aliphatic (acrylate)     | Weak/None   | Moderate       |
| UV-absorbing chromophore | Strong      | Variable       |
| Non-absorbing            | None        | Based on dn/dc |

### UV/RI Ratio Interpretation

![](copolymer-analysis_files/figure-html/uv-ri-concept-1.png)

## Example Dataset

``` r
data(sec_triple_detect, package = "measure.sec")

# Select copolymer samples
copolymers <- sec_triple_detect |>
  filter(polymer_type == "copolymer")

glimpse(copolymers)
#> Rows: 4,002
#> Columns: 11
#> $ sample_id        <chr> "Copoly-A", "Copoly-A", "Copoly-A", "Copoly-A", "Copo…
#> $ sample_type      <chr> "sample", "sample", "sample", "sample", "sample", "sa…
#> $ polymer_type     <chr> "copolymer", "copolymer", "copolymer", "copolymer", "…
#> $ elution_time     <dbl> 5.00, 5.01, 5.02, 5.03, 5.04, 5.05, 5.06, 5.07, 5.08,…
#> $ ri_signal        <dbl> 6.128782e-05, 5.685648e-04, 0.000000e+00, 4.859711e-0…
#> $ uv_signal        <dbl> 0.0000000000, 0.0004169298, 0.0000000000, 0.000000000…
#> $ mals_signal      <dbl> 4.816355e-06, 5.069294e-06, 3.309838e-06, 0.000000e+0…
#> $ known_mw         <dbl> 40000, 40000, 40000, 40000, 40000, 40000, 40000, 4000…
#> $ known_dispersity <dbl> 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5…
#> $ dn_dc            <dbl> 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15,…
#> $ extinction_coef  <dbl> 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6…
```

## UV/RI Ratio Analysis

### Basic Workflow

``` r
# UV/RI ratio analysis for compositional heterogeneity
rec_ratio <- recipe(
  ri_signal + uv_signal + elution_time ~ sample_id,
  data = copolymers
) |>
  update_role(sample_id, new_role = "id") |>
  # Convert signals to measure format
  step_measure_input_long(
    ri_signal,
    location = vars(elution_time),
    col_name = "ri"
  ) |>
  step_measure_input_long(
    uv_signal,
    location = vars(elution_time),
    col_name = "uv"
  ) |>
  # Baseline correction
  step_sec_baseline(measures = c("ri", "uv")) |>
  # Calculate UV/RI ratio
  step_sec_uv_ri_ratio(
    uv_col = "uv",
    ri_col = "ri",
    smooth = TRUE,          # Apply smoothing for noise reduction
    smooth_span = 0.1,      # Smoothing window (fraction of data)
    min_signal = 0.01       # Minimum signal threshold
  )

prepped_ratio <- prep(rec_ratio)
result_ratio <- bake(prepped_ratio, new_data = NULL)

# The result contains a uv_ri_ratio column with the ratio curve
```

### Extracting the Ratio Curve

``` r
# Get the ratio values at each elution time
ratio_data <- result_ratio |>
  select(sample_id, uv_ri_ratio) |>
  tidyr::unnest(uv_ri_ratio)

# Plot ratio vs elution time
ggplot(ratio_data, aes(location, value)) +
  geom_line() +
  facet_wrap(~sample_id) +
  labs(
    x = "Elution Time (min)",
    y = "UV/RI Ratio",
    title = "Compositional Profile"
  ) +
  theme_minimal()
```

## Composition Calculation

### Using Known Response Factors

When you know the response factors for each monomer, calculate actual
composition:

``` r
# Composition calculation with known response factors
rec_comp <- recipe(
  ri_signal + uv_signal + elution_time ~ sample_id,
  data = copolymers
) |>
  update_role(sample_id, new_role = "id") |>
  step_measure_input_long(ri_signal, location = vars(elution_time), col_name = "ri") |>
  step_measure_input_long(uv_signal, location = vars(elution_time), col_name = "uv") |>
  step_sec_baseline(measures = c("ri", "uv")) |>
  # Calculate composition from UV/RI
  step_sec_composition(
    uv_col = "uv",
    ri_col = "ri",
    # Response factors for component A (e.g., styrene)
    component_a_uv = 1.0,    # UV extinction coefficient
    component_a_ri = 0.185,  # RI response (dn/dc)
    # Response factors for component B (e.g., butadiene)
    component_b_uv = 0.1,    # UV extinction coefficient
    component_b_ri = 0.084   # RI response (dn/dc)
  )

prepped_comp <- prep(rec_comp)
result_comp <- bake(prepped_comp, new_data = NULL)

# Result contains:
# - composition_a: Weight fraction of component A at each point
# - composition_b: Weight fraction of component B at each point
```

### Response Factor Determination

Response factors should be measured experimentally:

1.  **dn/dc**: Measure with differential refractometer
2.  **UV extinction**: Measure with UV-Vis spectrophotometer

Typical values (THF, 25°C):

| Polymer       | dn/dc (mL/g) | ε₂₅₄ (mL/(g·cm)) |
|---------------|--------------|------------------|
| Polystyrene   | 0.185        | 1.0              |
| PMMA          | 0.084        | ~0.01            |
| Polybutadiene | 0.127        | 0.05             |
| PEG           | 0.069        | ~0               |

## Interpreting Results

### Uniform Composition

![](copolymer-analysis_files/figure-html/uniform-concept-1.png)

**Characteristics:** - Constant UV/RI ratio - Flat composition profile -
Indicates random or alternating copolymer

### Compositional Drift

![](copolymer-analysis_files/figure-html/drift-concept-1.png)

**Characteristics:** - Changing UV/RI ratio - Sloped composition
profile - Indicates gradient copolymer or reactivity differences

### Bimodal Composition

![](copolymer-analysis_files/figure-html/bimodal-concept-1.png)

**Characteristics:** - Abrupt ratio changes - Step-like composition
profile - Indicates blend or block copolymer

## Complete Workflow Example

``` r
# Full copolymer analysis workflow
rec_full <- recipe(
  ri_signal + uv_signal + elution_time + dn_dc + extinction_coef ~ sample_id,
  data = copolymers
) |>
  update_role(sample_id, new_role = "id") |>
  # Input signals
  step_measure_input_long(ri_signal, location = vars(elution_time), col_name = "ri") |>
  step_measure_input_long(uv_signal, location = vars(elution_time), col_name = "uv") |>

  # Preprocessing
  step_sec_detector_delay(reference = "ri", delay_volumes = c(uv = -0.05)) |>
  step_sec_baseline(measures = c("ri", "uv")) |>

  # Detector processing
  step_sec_ri(measures = "ri", dn_dc_column = "dn_dc") |>
  step_sec_uv(measures = "uv", extinction_column = "extinction_coef") |>

  # Composition analysis
  step_sec_uv_ri_ratio(uv_col = "uv", ri_col = "ri", smooth = TRUE) |>
  step_sec_composition(
    uv_col = "uv",
    ri_col = "ri",
    component_a_uv = 1.0,
    component_a_ri = 0.185,
    component_b_uv = 0.05,
    component_b_ri = 0.084
  ) |>

  # MW calculations
  step_sec_conventional_cal(
    standards = ps_standards,
    fit_type = "cubic"
  ) |>
  step_sec_mw_averages(measures = "log_mw")

prepped_full <- prep(rec_full)
result_full <- bake(prepped_full, new_data = NULL)

# Report composition statistics
result_full |>
  select(
    sample_id,
    mw_mw,
    composition_a_mean,
    composition_a_sd
  )
```

## Best Practices

### Detector Calibration

1.  **Verify response factors** with homopolymer standards
2.  **Account for detector delay** between UV and RI
3.  **Use same solvent** as for response factor determination

### Data Quality

1.  **Adequate signal-to-noise** on both detectors
2.  **Baseline correction** before ratio calculation
3.  **Signal threshold** to avoid noise in low-signal regions

### Interpretation Caveats

1.  **UV/RI ratio is relative**, not absolute composition
2.  **Assumes no specific interactions** between detector and polymer
3.  **Block copolymers** may show different behavior than random

## Session Info

``` r
sessionInfo()
#> R version 4.5.2 (2025-10-31)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.3 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
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
#> [1] ggplot2_4.0.1          measure.sec_0.0.0.9000 measure_0.0.1.9001    
#> [4] recipes_1.3.1          dplyr_1.1.4           
#> 
#> loaded via a namespace (and not attached):
#>  [1] gtable_0.3.6        xfun_0.55           bslib_0.9.0        
#>  [4] lattice_0.22-7      vctrs_0.6.5         tools_4.5.2        
#>  [7] generics_0.1.4      parallel_4.5.2      tibble_3.3.0       
#> [10] pkgconfig_2.0.3     Matrix_1.7-4        data.table_1.18.0  
#> [13] RColorBrewer_1.1-3  S7_0.2.1            desc_1.4.3         
#> [16] lifecycle_1.0.4     compiler_4.5.2      farver_2.1.2       
#> [19] textshaping_1.0.4   codetools_0.2-20    htmltools_0.5.9    
#> [22] class_7.3-23        sass_0.4.10         yaml_2.3.12        
#> [25] prodlim_2025.04.28  tidyr_1.3.2         pillar_1.11.1      
#> [28] pkgdown_2.2.0       jquerylib_0.1.4     MASS_7.3-65        
#> [31] cachem_1.1.0        gower_1.0.2         rpart_4.1.24       
#> [34] parallelly_1.46.0   lava_1.8.2          tidyselect_1.2.1   
#> [37] digest_0.6.39       future_1.68.0       purrr_1.2.0        
#> [40] listenv_0.10.0      labeling_0.4.3      splines_4.5.2      
#> [43] fastmap_1.2.0       grid_4.5.2          cli_3.6.5          
#> [46] magrittr_2.0.4      survival_3.8-3      future.apply_1.20.1
#> [49] withr_3.0.2         scales_1.4.0        lubridate_1.9.4    
#> [52] timechange_0.3.0    rmarkdown_2.30      globals_0.18.0     
#> [55] nnet_7.3-20         timeDate_4051.111   ragg_1.5.0         
#> [58] evaluate_1.0.5      knitr_1.51          hardhat_1.4.2      
#> [61] rlang_1.1.6         Rcpp_1.1.0          glue_1.8.0         
#> [64] ipred_0.9-15        jsonlite_2.0.0      R6_2.6.1           
#> [67] systemfonts_1.3.1   fs_1.6.6
```
