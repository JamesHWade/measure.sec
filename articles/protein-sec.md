# Protein SEC Analysis

## Overview

Size Exclusion Chromatography (SEC) is a critical analytical technique
for biopharmaceutical characterization. It separates proteins based on
hydrodynamic size, enabling quantification of:

- **High Molecular Weight Species (HMWS)**: Aggregates, dimers,
  oligomers
- **Monomer**: The intended product
- **Low Molecular Weight Species (LMWS)**: Fragments, degradation
  products

This vignette covers:

1.  Protein SEC basics
2.  Aggregate quantitation workflows
3.  Detailed oligomer analysis
4.  Regulatory considerations

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

## Protein SEC Principles

### Elution Order

In SEC, proteins elute in order of decreasing size:

![](protein-sec_files/figure-html/elution-concept-1.png)

### Detection

Protein SEC typically uses UV detection at 280 nm (tryptophan/tyrosine
absorbance).

## Basic Aggregate Quantitation

### Using step_sec_aggregates()

For simple HMWS/monomer/LMWS quantitation:

``` r
# Load or create protein SEC data
protein_data <- sec_triple_detect |>
  filter(sample_type == "sample") |>
  head(1)

# Note: Aggregate analysis requires defined peak boundaries
# The monomer_start and monomer_end values depend on your chromatographic conditions
rec_agg <- recipe(
  uv_signal + elution_time ~ sample_id,
  data = protein_data
) |>
  update_role(sample_id, new_role = "id") |>
  # Convert UV signal to measure format
  step_measure_input_long(
    uv_signal,
    location = vars(elution_time),
    col_name = "uv"
  ) |>
  # Baseline correction
  step_sec_baseline(measures = "uv") |>
  # Aggregate quantitation
  step_sec_aggregates(
    measures = "uv",
    monomer_start = 9.5,   # Monomer peak start
    monomer_end = 11.0,    # Monomer peak end
    method = "manual"
  )

prepped_agg <- prep(rec_agg)
result_agg <- bake(prepped_agg, new_data = NULL)

# View aggregate results
result_agg |>
  select(sample_id, purity_hmws, purity_monomer, purity_lmws)
```

### Automatic Peak Detection

Let the algorithm find the tallest peak as the monomer:

``` r
step_sec_aggregates(
  measures = "uv",
  method = "tallest"  # Uses tallest peak as monomer reference
)
```

## Complete Protein Workflow

### Using step_sec_protein()

The
[`step_sec_protein()`](https://jameshwade.github.io/measure-sec/reference/step_sec_protein.md)
function provides a streamlined workflow combining baseline correction,
aggregate analysis, and optional oligomer detection:

``` r
# Comprehensive protein SEC workflow
# Note: step_sec_protein provides a streamlined workflow for protein analysis
rec_protein <- recipe(
  uv_signal + elution_time ~ sample_id,
  data = protein_data
) |>
  update_role(sample_id, new_role = "id") |>
  step_measure_input_long(
    uv_signal,
    location = vars(elution_time),
    col_name = "uv"
  ) |>
  step_sec_protein(
    type = "native",
    monomer_mw = 150000,        # mAb ~150 kDa
    extinction_coef = 1.4,      # E1% at 280 nm
    include_oligomer = TRUE,
    baseline_method = "linear"
  )

prepped_protein <- prep(rec_protein)
result_protein <- bake(prepped_protein, new_data = NULL)

# View comprehensive results
result_protein |>
  select(
    sample_id,
    protein_hmws_pct,
    protein_monomer_pct,
    protein_lmws_pct,
    protein_dimer_pct,
    protein_species_count
  )
```

### Output Columns

The
[`step_sec_protein()`](https://jameshwade.github.io/measure-sec/reference/step_sec_protein.md)
step creates:

| Column                | Description                     |
|-----------------------|---------------------------------|
| `protein_hmws_pct`    | % high molecular weight species |
| `protein_monomer_pct` | % main peak (monomer)           |
| `protein_lmws_pct`    | % low molecular weight species  |
| `protein_main_start`  | Start of monomer region         |
| `protein_main_end`    | End of monomer region           |

With `include_oligomer = TRUE`:

| Column                      | Description                      |
|-----------------------------|----------------------------------|
| `protein_monomer_oligo_pct` | % monomer from oligomer analysis |
| `protein_dimer_pct`         | % dimer                          |
| `protein_trimer_pct`        | % trimer                         |
| `protein_hmw_oligo_pct`     | % higher oligomers               |
| `protein_lmw_oligo_pct`     | % fragments                      |
| `protein_species_count`     | Number of detected species       |

## Native vs Denaturing SEC

### Native SEC

Preserves non-covalent interactions:

``` r
step_sec_protein(
  type = "native",
  monomer_mw = 150000
)
```

Use for: - Detecting reversible aggregates - Oligomer state assessment -
Native quaternary structure

### Denaturing SEC (SDS-SEC)

Disrupts non-covalent interactions:

``` r
step_sec_protein(
  type = "denaturing",
  monomer_mw = 150000
)
```

Use for: - Covalent aggregate detection - Clipped/truncated species -
Heavy/light chain analysis

## Detailed Oligomer Analysis

### Using step_sec_oligomer()

For more control over oligomer detection:

``` r
# Detailed oligomer analysis with explicit control
rec_oligo <- recipe(
  uv_signal + elution_time ~ sample_id,
  data = protein_data
) |>
  update_role(sample_id, new_role = "id") |>
  step_measure_input_long(
    uv_signal,
    location = vars(elution_time),
    col_name = "uv"
  ) |>
  step_sec_baseline(measures = "uv") |>
  step_sec_oligomer(
    measures = "uv",
    monomer_mw = 150000,
    mw_tolerance = 0.15,     # 15% MW tolerance for species assignment
    min_peak_area = 0.001    # Minimum 0.1% to report
  )
```

### Species Assignment

The oligomer step identifies species based on expected MW:

| Species         | Expected MW      | Example (mAb) |
|-----------------|------------------|---------------|
| Fragment        | \< 0.7 × monomer | \< 105 kDa    |
| Monomer         | 1.0 × monomer    | 150 kDa       |
| Dimer           | 2.0 × monomer    | 300 kDa       |
| Trimer          | 3.0 × monomer    | 450 kDa       |
| Higher oligomer | \> 3 × monomer   | \> 450 kDa    |

## Regulatory Considerations

### ICH Guidelines

- **ICH Q6B**: Specifications for biotechnology products
  - Aggregate content is a critical quality attribute
  - Typical specification: HMWS \< 5%
- **ICH Q5E**: Comparability of biotechnology products
  - Aggregate profiles should be comparable

### USP Guidelines

- **USP \<129\>**: Analytical procedures for proteins
- **USP \<1032\>**: Design and development of biological assays

### Typical Specifications

| Parameter | Typical Limit | Notes                                 |
|-----------|---------------|---------------------------------------|
| HMWS      | \< 5%         | May be tighter for high-dose products |
| Monomer   | \> 95%        | Primary quality indicator             |
| LMWS      | \< 2%         | Product-specific                      |
| Dimer     | Report        | Not always specified separately       |

## Best Practices

### Sample Preparation

1.  **Filter samples** through 0.22 μm before injection
2.  **Minimize time** between prep and analysis
3.  **Use appropriate buffer** (match formulation buffer)
4.  **Control temperature** during autosampler storage

### Method Parameters

1.  **Column selection**: Appropriate exclusion limit for aggregates
2.  **Mobile phase**: Typically PBS or phosphate buffer
3.  **Flow rate**: 0.5-1.0 mL/min
4.  **Injection volume**: 10-100 μL (depending on concentration)
5.  **Detection**: 280 nm (standard) or 214 nm (higher sensitivity)

### Integration

1.  **Consistent baseline** selection across samples
2.  **Same integration limits** for all samples in a study
3.  **Report detection limit** for minor species

## Example: Stability Study Analysis

``` r
# Stability data structure example
# In practice, chrom1-4 and time would be your actual chromatogram data vectors
# For example:
#   time <- seq(5, 20, by = 0.05)
#   chrom1 <- your_uv_signal_at_T0
#   chrom2 <- your_uv_signal_at_T1M
#   etc.

stability_data <- tibble(
  sample_id = c("T0", "T1M", "T3M", "T6M"),
  timepoint = c(0, 1, 3, 6),
  uv_signal = list(chrom1, chrom2, chrom3, chrom4),  # Your chromatogram vectors
  elution_time = list(time, time, time, time)        # Shared time axis
)

# Analyze all timepoints
rec <- recipe(~., data = stability_data) |>
  step_measure_input_long(
    uv_signal,
    location = vars(elution_time),
    col_name = "uv"
  ) |>
  step_sec_protein(
    monomer_mw = 150000,
    include_oligomer = TRUE
  )

prepped <- prep(rec)
result <- bake(prepped, new_data = NULL)

# Plot stability trend
ggplot(result, aes(timepoint, protein_monomer_pct)) +
  geom_line(linewidth = 1, color = "#2E86AB") +
  geom_point(size = 3, color = "#2E86AB") +
  geom_hline(yintercept = 95, linetype = "dashed", color = "red") +
  annotate("text", x = 5, y = 94.5, label = "Specification: > 95%", color = "red") +
  labs(
    x = "Time (months)",
    y = "% Monomer",
    title = "Stability Profile: Monomer Content"
  ) +
  ylim(90, 100) +
  theme_minimal()
```

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
