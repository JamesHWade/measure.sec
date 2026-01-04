# Exporting Results

## Overview

After processing SEC data with measure.sec, you’ll need to export
results for reports, further analysis, or regulatory submissions. This
guide covers:

1.  Creating summary tables with MW averages
2.  Extracting slice-by-slice data for detailed analysis
3.  Comparing multiple samples side-by-side
4.  Generating automated reports
5.  Exporting to Excel, CSV, and other formats

## Quick Reference: Export Functions

| Function                                                                                                         | Purpose                 | Output                              |
|------------------------------------------------------------------------------------------------------------------|-------------------------|-------------------------------------|
| [`measure_sec_summary_table()`](https://jameshwade.github.io/measure-sec/reference/measure_sec_summary_table.md) | Per-sample MW averages  | Tibble with Mn, Mw, Mz, dispersity  |
| [`measure_sec_slice_table()`](https://jameshwade.github.io/measure-sec/reference/measure_sec_slice_table.md)     | Point-by-point data     | Long or wide format tibble          |
| [`measure_sec_compare()`](https://jameshwade.github.io/measure-sec/reference/measure_sec_compare.md)             | Multi-sample comparison | Summary, differences, optional plot |
| [`measure_sec_report()`](https://jameshwade.github.io/measure-sec/reference/measure_sec_report.md)               | Automated reports       | HTML, PDF, or Word document         |

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

## Creating Processed Data for Export

Before exporting, we need processed SEC data. Let’s process a single
sample to demonstrate the export functions:

``` r
# Load multi-detector SEC data
data(sec_triple_detect)

# Filter to a single sample for processing
ps_sample <- sec_triple_detect |>
  filter(sample_id == "PS-100K")

# View what we're working with
ps_sample |>
  select(sample_id, polymer_type, known_mw, known_dispersity) |>
  distinct()
#> # A tibble: 1 × 4
#>   sample_id polymer_type known_mw known_dispersity
#>   <chr>     <chr>           <dbl>            <dbl>
#> 1 PS-100K   polystyrene    100000             1.01
```

``` r
# Create and prep a recipe for one sample
# Note: the formula tells recipes how to group the data
rec <- recipe(
  ri_signal + elution_time + known_mw ~ sample_id,
  data = ps_sample
) |>
  update_role(sample_id, new_role = "id") |>
  step_measure_input_long(
    ri_signal,
    location = vars(elution_time),
    col_name = "ri"
  ) |>
  step_sec_baseline(measures = "ri")

prepped <- prep(rec)
processed <- bake(prepped, new_data = NULL)

# View the structure - now has nested measure_list column
processed
#> # A tibble: 1 × 4
#>   sample_id known_mw          ri elution_time 
#>   <chr>        <dbl>      <meas> <list>       
#> 1 PS-100K     100000 [2,001 × 2] <dbl [2,001]>
```

For demonstrating multi-sample comparisons later, let’s process a few
more samples:

``` r
# Helper function to process one sample
process_sample <- function(sample_id_filter) {
  data <- sec_triple_detect |>
    filter(sample_id == sample_id_filter)

  recipe(
    ri_signal + elution_time + known_mw ~ sample_id,
    data = data
  ) |>
    update_role(sample_id, new_role = "id") |>
    step_measure_input_long(ri_signal, location = vars(elution_time), col_name = "ri") |>
    step_sec_baseline(measures = "ri") |>
    prep() |>
    bake(new_data = NULL)
}

# Process three samples
ps100k <- process_sample("PS-100K")
ps500k <- process_sample("PS-500K")
pmma_high <- process_sample("PMMA-High")
```

## Summary Tables

### Basic Summary Table

[`measure_sec_summary_table()`](https://jameshwade.github.io/measure-sec/reference/measure_sec_summary_table.md)
creates a one-row-per-sample table with key metrics:

``` r
# Create summary table
summary_tbl <- measure_sec_summary_table(
  processed,
  sample_id = "sample_id"
)

print(summary_tbl)
#> SEC Analysis Summary
#> ============================================================ 
#> Samples: 1 
#> 
#> # A tibble: 1 × 1
#>   sample_id
#>   <chr>    
#> 1 PS-100K
```

### Including Additional Columns

You can include any numeric column from your data:

``` r
# Include known MW for validation
summary_with_known <- measure_sec_summary_table(
  processed,
  sample_id = "sample_id",
  additional_cols = c("known_mw")
)

print(summary_with_known)
#> SEC Analysis Summary
#> ============================================================ 
#> Samples: 1 
#> 
#> # A tibble: 1 × 2
#>   sample_id known_mw
#>   <chr>        <dbl>
#> 1 PS-100K     100000
```

### Controlling Decimal Places

For regulatory submissions that require specific precision:

``` r
summary_precise <- measure_sec_summary_table(
  processed,
  sample_id = "sample_id",
  digits = 4
)
```

## Slice Tables: Point-by-Point Data

### Understanding Slice Data

SEC analysis works on a point-by-point basis across the chromatogram.
Each “slice” represents one data point with its elution time and
corresponding signal values. Use
[`measure_sec_slice_table()`](https://jameshwade.github.io/measure-sec/reference/measure_sec_slice_table.md)
to extract this detailed data.

### Long Format (Default)

Long format is best for plotting with ggplot2:

``` r
# Extract slice data in long format
slices_long <- measure_sec_slice_table(
  processed,
  measures = "ri",
  sample_id = "sample_id"
)

# View structure
head(slices_long, 10)
#> # A tibble: 10 × 5
#>    sample_id slice location measure   value
#>    <chr>     <int>    <dbl> <chr>     <dbl>
#>  1 PS-100K       1     5    ri      0.00132
#>  2 PS-100K       2     5.01 ri      0.00189
#>  3 PS-100K       3     5.02 ri      0.00156
#>  4 PS-100K       4     5.03 ri      0.00123
#>  5 PS-100K       5     5.04 ri      0.00157
#>  6 PS-100K       6     5.05 ri      0.00116
#>  7 PS-100K       7     5.06 ri      0.00188
#>  8 PS-100K       8     5.07 ri      0.00110
#>  9 PS-100K       9     5.08 ri      0.00135
#> 10 PS-100K      10     5.09 ri      0.00126
```

``` r
# Plot using the slice table
ggplot(slices_long, aes(x = location, y = value, color = sample_id)) +
  geom_line() +
  labs(
    x = "Elution Time (min)",
    y = "RI Signal",
    title = "Chromatogram Overlay from Slice Data"
  ) +
  theme_minimal()
```

![](exporting-results_files/figure-html/plot-slices-1.png)

### Wide Format

Wide format is better for spreadsheet export or correlating multiple
measures:

``` r
# If you had multiple measures (e.g., ri and mw), wide format would look like:
slices_wide <- measure_sec_slice_table(
  processed,
  measures = c("ri"),
  sample_id = "sample_id",
  pivot = TRUE
)

head(slices_wide)
```

### Exporting Slice Data to CSV

``` r
# Export for external analysis
write.csv(slices_long, "sec_slice_data.csv", row.names = FALSE)

# Or use readr for consistent formatting
readr::write_csv(slices_long, "sec_slice_data.csv")
```

## Comparing Multiple Samples

### Basic Comparison

[`measure_sec_compare()`](https://jameshwade.github.io/measure-sec/reference/measure_sec_compare.md)
provides side-by-side comparison with differences from a reference:

``` r
# Compare the samples we processed earlier
comparison <- measure_sec_compare(
  ps100k, ps500k, pmma_high,
  samples = c("PS 100K", "PS 500K", "PMMA High"),
  metrics = "mw_averages",
  plot = FALSE
)

print(comparison)
#> SEC Multi-Sample Comparison
#> ============================================================ 
#> Samples: 3 
#> Reference: PS 100K 
#> 
#> Summary:
#> # A tibble: 3 × 1
#>   sample   
#>   <chr>    
#> 1 PS 100K  
#> 2 PS 500K  
#> 3 PMMA High
```

### Understanding the Comparison Output

The comparison object contains:

- **`$summary`**: All samples with their metrics
- **`$differences`**: Absolute and percent differences from reference
- **`$plot`**: MWD overlay (if requested and available)
- **`$reference`**: Which sample is the reference

``` r
# Access individual components
comparison$summary      # Metrics for all samples
comparison$differences  # Differences from reference
comparison$reference    # Reference sample name
```

### Setting a Different Reference

By default, the first sample is the reference. You can change this:

``` r
# Use PS 500K as reference
comparison_ref <- measure_sec_compare(
  ps100k, ps500k, pmma,
  samples = c("PS 100K", "PS 500K", "PMMA High"),
  reference = "PS 500K"
)
```

### Batch-to-Batch Comparison

A common use case is comparing production batches to a reference lot:

``` r
# Compare production batches
batch_comparison <- measure_sec_compare(
  reference_lot,
  batch_001,
  batch_002,
  batch_003,
  samples = c("Reference", "Batch 001", "Batch 002", "Batch 003"),
  reference = "Reference"
)

# Check for significant deviations
batch_comparison$differences |>
  filter(abs(mw_mw_pct) > 5)  # Flag batches > 5% different
```

## Automated Reports

### Available Templates

``` r
# See available report templates
list_sec_templates()
#> # A tibble: 3 × 3
#>   template description                                         formats        
#>   <chr>    <chr>                                               <chr>          
#> 1 standard Summary table, chromatogram, and MWD plot           html, pdf, docx
#> 2 detailed All plots, multi-detector view, optional slice data html, pdf, docx
#> 3 qc       System suitability with pass/fail metrics           html, pdf, docx
```

### Standard Report

The standard template includes: - Summary table with MW averages -
Chromatogram overlay - Molecular weight distribution plot

``` r
# Generate HTML report
measure_sec_report(
  processed,
  template = "standard",
  output_format = "html",
  title = "Polymer SEC Analysis",
  author = "Lab Analyst"
)
```

### Detailed Report

For more comprehensive documentation:

``` r
measure_sec_report(
  processed,
  template = "detailed",
  output_format = "pdf",
  title = "Comprehensive SEC Report",
  sample_id = "sample_id",
  include_slice_table = TRUE  # Append raw data
)
```

### QC Report

For system suitability testing:

``` r
measure_sec_report(
  sst_data,
  template = "qc",
  output_format = "html",
  specs = list(
    plate_count_min = 10000,
    asymmetry_min = 0.8,
    asymmetry_max = 1.5,
    resolution_min = 1.5
  )
)
```

### Saving Reports to Specific Locations

``` r
# Save to a specific path
measure_sec_report(
  processed,
  template = "standard",
  output_file = "reports/polymer_analysis_2024-01-15.html",
  open = FALSE  # Don't open automatically
)
```

## Integration with Other Tools

### Export to Excel

For spreadsheet-based workflows:

``` r
# Requires writexl package
library(writexl)

# Create a workbook with multiple sheets
write_xlsx(
  list(
    "Summary" = summary_tbl,
    "Slice Data" = slices_long
  ),
  "sec_results.xlsx"
)
```

### Export to Database

For LIMS integration:

``` r
# Using DBI for database connection
library(DBI)

con <- dbConnect(RSQLite::SQLite(), "lims.db")

# Write summary to database
dbWriteTable(con, "sec_results", summary_tbl, append = TRUE)

dbDisconnect(con)
```

### Export for GraphPad Prism

Prism prefers wide format with specific column arrangements:

``` r
# Create Prism-friendly format
prism_data <- slices_long |>
  select(location, sample_id, value) |>
  tidyr::pivot_wider(names_from = sample_id, values_from = value)

write.csv(prism_data, "sec_for_prism.csv", row.names = FALSE)
```

### Export for Python Analysis

For interoperability with Python workflows:

``` r
# Save as feather for fast Python loading
arrow::write_feather(slices_long, "sec_data.feather")

# Or as parquet for columnar storage
arrow::write_parquet(slices_long, "sec_data.parquet")
```

## Best Practices

### Traceability

Always include metadata for regulatory compliance:

``` r
# Add audit trail information
summary_with_audit <- summary_tbl |>
  mutate(
    analysis_date = Sys.Date(),
    analyst = "JW",
    instrument = "Agilent 1260",
    column = "PLgel Mixed-C",
    software_version = packageVersion("measure.sec")
  )
```

### File Naming Conventions

Use descriptive, consistent file names:

``` r
# Good: includes key information
"sec_summary_batch123_2024-01-15.xlsx"
"mwd_comparison_stability_t12m.pdf"

# Avoid: ambiguous names
"results.xlsx"
"output.csv"
```

### Version Control for Reports

Save report parameters for reproducibility:

``` r
# Save analysis parameters alongside results
analysis_params <- list(
  date = Sys.time(),
  calibration = "ps_2024-01.rds",
  baseline_method = "linear",
  integration_limits = c(8, 18),
  package_version = packageVersion("measure.sec")
)

saveRDS(analysis_params, "sec_analysis_params.rds")
```

## Troubleshooting

### Missing Columns in Summary

If MW columns are missing from your summary table, check that you’ve run
the calibration and MW averaging steps:

``` r
# Ensure you have MW data
names(processed)  # Should include Mn, Mw, Mz, or mw_mn, mw_mw, mw_mz
```

### Empty Slice Tables

If
[`measure_sec_slice_table()`](https://jameshwade.github.io/measure-sec/reference/measure_sec_slice_table.md)
returns no data:

``` r
# Check that measure columns exist
find_measure_cols <- function(data) {
  names(data)[vapply(data, inherits, logical(1), "measure_list")]
}
find_measure_cols(processed)
```

### Report Generation Fails

If
[`measure_sec_report()`](https://jameshwade.github.io/measure-sec/reference/measure_sec_report.md)
fails:

1.  Verify Quarto is installed:
    [`quarto::quarto_version()`](https://quarto-dev.github.io/quarto-r/reference/quarto_version.html)
2.  For PDF output, ensure LaTeX is installed
3.  Check data has required columns

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
#> [46] magrittr_2.0.4      utf8_1.2.6          survival_3.8-3     
#> [49] future.apply_1.20.1 withr_3.0.2         scales_1.4.0       
#> [52] lubridate_1.9.4     timechange_0.3.0    rmarkdown_2.30     
#> [55] globals_0.18.0      nnet_7.3-20         timeDate_4051.111  
#> [58] ragg_1.5.0          evaluate_1.0.5      knitr_1.51         
#> [61] hardhat_1.4.2       rlang_1.1.6         Rcpp_1.1.0         
#> [64] glue_1.8.0          ipred_0.9-15        jsonlite_2.0.0     
#> [67] R6_2.6.1            systemfonts_1.3.1   fs_1.6.6
```
