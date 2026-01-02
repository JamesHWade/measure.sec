# Getting Started with measure.sec

## Overview

**measure.sec** provides preprocessing steps for Size Exclusion
Chromatography (SEC) and Gel Permeation Chromatography (GPC) data
analysis. It extends the
[measure](https://github.com/JamesHWade/measure) package using the
[recipes](https://recipes.tidymodels.org/) framework.

This vignette covers:

1.  Installation and setup
2.  Understanding the data model
3.  A basic single-detector workflow
4.  Calculating molecular weight averages

## Installation

``` r
# Install from GitHub
# install.packages("pak")
pak::pak("JamesHWade/measure")
pak::pak("JamesHWade/measure-sec")
```

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

## The Data Model

SEC data in measure.sec uses the **measure** package’s nested tibble
structure:

- **`measure_tbl`**: A single chromatogram with `location` (elution
  time/volume) and `value` (detector response)
- **`measure_list`**: A list column containing multiple `measure_tbl`
  objects

This structure allows you to store complete chromatograms alongside
sample metadata.

## Example Dataset

The package includes `sec_triple_detect`, a synthetic multi-detector SEC
dataset:

``` r
data(sec_triple_detect, package = "measure.sec")

# Overview
glimpse(sec_triple_detect)
#> Rows: 24,012
#> Columns: 11
#> $ sample_id        <chr> "PS-1K", "PS-1K", "PS-1K", "PS-1K", "PS-1K", "PS-1K",…
#> $ sample_type      <chr> "standard", "standard", "standard", "standard", "stan…
#> $ polymer_type     <chr> "polystyrene", "polystyrene", "polystyrene", "polysty…
#> $ elution_time     <dbl> 5.00, 5.01, 5.02, 5.03, 5.04, 5.05, 5.06, 5.07, 5.08,…
#> $ ri_signal        <dbl> 6.926392e-04, 0.000000e+00, 3.199253e-04, 4.197175e-0…
#> $ uv_signal        <dbl> 0.0002034583, 0.0000000000, 0.0000000000, 0.000000000…
#> $ mals_signal      <dbl> 3.370385e-05, 3.483481e-05, 3.102092e-05, 3.261962e-0…
#> $ known_mw         <dbl> 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000,…
#> $ known_dispersity <dbl> 1.05, 1.05, 1.05, 1.05, 1.05, 1.05, 1.05, 1.05, 1.05,…
#> $ dn_dc            <dbl> 0.185, 0.185, 0.185, 0.185, 0.185, 0.185, 0.185, 0.18…
#> $ extinction_coef  <dbl> 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2…
```

The dataset contains:

- 12 polymer samples (polystyrene, PMMA, PEG, copolymers)
- RI, UV, and MALS detector signals
- Known molecular weights and dispersities
- Sample-specific optical constants (dn/dc, extinction coefficients)

``` r
# Sample types
sec_triple_detect |>
  distinct(sample_id, sample_type, polymer_type) |>
  print(n = 12)
#> # A tibble: 12 × 3
#>    sample_id sample_type polymer_type
#>    <chr>     <chr>       <chr>       
#>  1 PS-1K     standard    polystyrene 
#>  2 PS-10K    standard    polystyrene 
#>  3 PS-50K    standard    polystyrene 
#>  4 PS-100K   standard    polystyrene 
#>  5 PS-500K   standard    polystyrene 
#>  6 PMMA-Low  sample      pmma        
#>  7 PMMA-Med  sample      pmma        
#>  8 PMMA-High sample      pmma        
#>  9 PEG-5K    sample      peg         
#> 10 PEG-20K   sample      peg         
#> 11 Copoly-A  sample      copolymer   
#> 12 Copoly-B  sample      copolymer
```

## Basic Workflow: RI Detector Analysis

Let’s analyze a polystyrene sample using the RI detector:

``` r
# Select a single polystyrene standard
ps_sample <- sec_triple_detect |>
  filter(sample_id == "PS-50K")

# View the data structure
ps_sample |>
  select(sample_id, polymer_type, known_mw, ri_signal)
#> # A tibble: 2,001 × 4
#>    sample_id polymer_type known_mw ri_signal
#>    <chr>     <chr>           <dbl>     <dbl>
#>  1 PS-50K    polystyrene     50000 0        
#>  2 PS-50K    polystyrene     50000 0.000279 
#>  3 PS-50K    polystyrene     50000 0        
#>  4 PS-50K    polystyrene     50000 0        
#>  5 PS-50K    polystyrene     50000 0.000842 
#>  6 PS-50K    polystyrene     50000 0.000483 
#>  7 PS-50K    polystyrene     50000 0.0000325
#>  8 PS-50K    polystyrene     50000 0        
#>  9 PS-50K    polystyrene     50000 0        
#> 10 PS-50K    polystyrene     50000 0.000828 
#> # ℹ 1,991 more rows
```

### Step 1: Create a Recipe

Recipes define a sequence of preprocessing steps. Start by converting
raw signal columns to the measure format:

``` r
rec <- recipe(~., data = ps_sample) |>
  # Convert RI signal to measure format
  step_measure_input_long(
    ri_signal,
    location = vars(elution_time),
    col_name = "ri"
  )
```

### Step 2: Add Preprocessing Steps

Add baseline correction and RI processing:

``` r
rec <- recipe(~., data = ps_sample) |>
  step_measure_input_long(
    ri_signal,
    location = vars(elution_time),
    col_name = "ri"
  ) |>
  # Baseline correction
  step_sec_baseline(measures = "ri") |>
  # RI detector processing with dn/dc
  step_sec_ri(measures = "ri", dn_dc_column = "dn_dc")
```

### Step 3: Prep and Bake

[`prep()`](https://recipes.tidymodels.org/reference/prep.html) learns
parameters from the training data,
[`bake()`](https://recipes.tidymodels.org/reference/bake.html) applies
the transformations:

``` r
prepped <- prep(rec)
result <- bake(prepped, new_data = NULL)

# View the processed data
result |>
  select(sample_id, ri)
```

## Molecular Weight Averages

Calculate Mn, Mw, Mz, and dispersity using
[`step_sec_mw_averages()`](https://jameshwade.github.io/measure-sec/reference/step_sec_mw_averages.md):

``` r
rec <- recipe(~., data = ps_sample) |>
  step_measure_input_long(
    ri_signal,
    location = vars(elution_time),
    col_name = "ri"
  ) |>
  step_sec_baseline(measures = "ri") |>
  step_sec_ri(measures = "ri", dn_dc_column = "dn_dc") |>
  # Calculate MW averages using known MW values
  step_sec_mw_averages(mw_column = "known_mw")

prepped <- prep(rec)
result <- bake(prepped, new_data = NULL)

# View molecular weight results
result |>
  select(sample_id, mw_mn, mw_mw, mw_mz, mw_dispersity)
```

## Calibration Curves

For samples without absolute MW data, use calibration standards:

``` r
# Load polystyrene standards
data(sec_ps_standards, package = "measure.sec")

# View the standards
sec_ps_standards |>
  select(standard_name, mp, log_mp, retention_time) |>
  print(n = 8)
#> # A tibble: 16 × 4
#>   standard_name      mp log_mp retention_time
#>   <chr>           <dbl>  <dbl>          <dbl>
#> 1 PS-3150000    3150000   6.50           11.2
#> 2 PS-1870000    1870000   6.27           11.6
#> 3 PS-1090000    1090000   6.04           12.1
#> 4 PS-630000      630000   5.80           12.6
#> 5 PS-430000      430000   5.63           13.2
#> 6 PS-216000      216000   5.33           13.8
#> 7 PS-120000      120000   5.08           14.3
#> 8 PS-67500        67500   4.83           15.0
#> # ℹ 8 more rows

# Visualize calibration curve
ggplot(sec_ps_standards, aes(retention_time, log_mp)) +
  geom_point(size = 3, color = "#2E86AB") +
  geom_smooth(
    method = "lm",
    formula = y ~ poly(x, 3),
    se = TRUE,
    color = "#A23B72",
    fill = "#A23B72",
    alpha = 0.2
  ) +
  labs(
    x = "Retention Time (min)",
    y = expression(log[10](M[p])),
    title = "Polystyrene Calibration Curve"
  ) +
  theme_minimal()
```

![](getting-started_files/figure-html/calibration-1.png)

Apply conventional calibration with
[`step_sec_conventional_cal()`](https://jameshwade.github.io/measure-sec/reference/step_sec_conventional_cal.md):

``` r
# Prepare standards for calibration
ps_cal <- sec_ps_standards |>
  select(retention = retention_time, log_mw = log_mp)

rec <- recipe(~., data = ps_sample) |>
  step_measure_input_long(
    ri_signal,
    location = vars(elution_time),
    col_name = "ri"
  ) |>
  step_sec_baseline(measures = "ri") |>
  # Apply calibration
  step_sec_conventional_cal(
    standards = ps_cal,
    fit_type = "cubic"
  ) |>
  step_sec_mw_averages(measures = "log_mw")

prepped <- prep(rec)
result <- bake(prepped, new_data = NULL)
```

## Available Steps

The package provides steps for: \### Preprocessing -
[`step_sec_baseline()`](https://jameshwade.github.io/measure-sec/reference/step_sec_baseline.md):
SEC-optimized baseline correction -
[`step_sec_detector_delay()`](https://jameshwade.github.io/measure-sec/reference/step_sec_detector_delay.md):
Correct inter-detector delays

### Detector Processing

- [`step_sec_ri()`](https://jameshwade.github.io/measure-sec/reference/step_sec_ri.md):
  RI detector with dn/dc
- [`step_sec_uv()`](https://jameshwade.github.io/measure-sec/reference/step_sec_uv.md):
  UV detector with extinction coefficient
- [`step_sec_mals()`](https://jameshwade.github.io/measure-sec/reference/step_sec_mals.md),
  [`step_sec_lals()`](https://jameshwade.github.io/measure-sec/reference/step_sec_lals.md),
  [`step_sec_rals()`](https://jameshwade.github.io/measure-sec/reference/step_sec_rals.md):
  Light scattering
- [`step_sec_dls()`](https://jameshwade.github.io/measure-sec/reference/step_sec_dls.md):
  Dynamic light scattering
- [`step_sec_viscometer()`](https://jameshwade.github.io/measure-sec/reference/step_sec_viscometer.md):
  Differential viscometer

### Molecular Weight

- [`step_sec_mw_averages()`](https://jameshwade.github.io/measure-sec/reference/step_sec_mw_averages.md):
  Mn, Mw, Mz, dispersity
- [`step_sec_mw_fractions()`](https://jameshwade.github.io/measure-sec/reference/step_sec_mw_fractions.md):
  MW fractions above/below cutoffs
- [`step_sec_mw_distribution()`](https://jameshwade.github.io/measure-sec/reference/step_sec_mw_distribution.md):
  Differential/cumulative MWD
- [`step_sec_conventional_cal()`](https://jameshwade.github.io/measure-sec/reference/step_sec_conventional_cal.md):
  Narrow standard calibration
- [`step_sec_universal_cal()`](https://jameshwade.github.io/measure-sec/reference/step_sec_universal_cal.md):
  Universal calibration

### Composition & Protein

- [`step_sec_uv_ri_ratio()`](https://jameshwade.github.io/measure-sec/reference/step_sec_uv_ri_ratio.md):
  UV/RI ratio for heterogeneity
- [`step_sec_composition()`](https://jameshwade.github.io/measure-sec/reference/step_sec_composition.md):
  Copolymer composition
- [`step_sec_aggregates()`](https://jameshwade.github.io/measure-sec/reference/step_sec_aggregates.md):
  HMWS/monomer/LMWS quantitation
- [`step_sec_protein()`](https://jameshwade.github.io/measure-sec/reference/step_sec_protein.md):
  Complete protein SEC workflow

## Next Steps

For more detailed workflows, see:

- [`vignette("triple-detection")`](https://jameshwade.github.io/measure-sec/articles/triple-detection.md):
  Multi-detector SEC with MALS
- [`vignette("protein-sec")`](https://jameshwade.github.io/measure-sec/articles/protein-sec.md):
  Protein aggregate analysis
- [`vignette("copolymer-analysis")`](https://jameshwade.github.io/measure-sec/articles/copolymer-analysis.md):
  Composition analysis
- [`vignette("system-suitability")`](https://jameshwade.github.io/measure-sec/articles/system-suitability.md):
  QC and SST testing
- [`vignette("sec-analysis")`](https://jameshwade.github.io/measure-sec/articles/sec-analysis.md):
  Comprehensive reference

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
#> [34] nlme_3.1-168        parallelly_1.46.0   lava_1.8.2         
#> [37] tidyselect_1.2.1    digest_0.6.39       future_1.68.0      
#> [40] purrr_1.2.0         listenv_0.10.0      labeling_0.4.3     
#> [43] splines_4.5.2       fastmap_1.2.0       grid_4.5.2         
#> [46] cli_3.6.5           magrittr_2.0.4      utf8_1.2.6         
#> [49] survival_3.8-3      future.apply_1.20.1 withr_3.0.2        
#> [52] scales_1.4.0        lubridate_1.9.4     timechange_0.3.0   
#> [55] rmarkdown_2.30      globals_0.18.0      nnet_7.3-20        
#> [58] timeDate_4051.111   ragg_1.5.0          evaluate_1.0.5     
#> [61] knitr_1.51          hardhat_1.4.2       mgcv_1.9-3         
#> [64] rlang_1.1.6         Rcpp_1.1.0          glue_1.8.0         
#> [67] ipred_0.9-15        jsonlite_2.0.0      R6_2.6.1           
#> [70] systemfonts_1.3.1   fs_1.6.6
```
