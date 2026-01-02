# Generate SEC Summary Table

Creates a summary table of SEC analysis results with key metrics for
each sample.

## Usage

``` r
measure_sec_summary_table(
  data,
  mw_col = NULL,
  include_mw = TRUE,
  include_fractions = TRUE,
  include_purity = TRUE,
  sample_id = NULL,
  additional_cols = NULL,
  digits = 2
)
```

## Arguments

- data:

  A data frame containing SEC results.

- mw_col:

  Column name containing molecular weight averages (list column with Mn,
  Mw, Mz, dispersity).

- include_mw:

  Logical. Include molecular weight averages? Default is `TRUE`.

- include_fractions:

  Logical. Include MW fractions if available? Default is `TRUE`.

- include_purity:

  Logical. Include purity metrics (HMWS, monomer, LMWS) if available?
  Default is `TRUE`.

- sample_id:

  Column name for sample identifiers.

- additional_cols:

  Character vector of additional columns to include in the summary.

- digits:

  Number of decimal places for numeric columns. Default is 2.

## Value

A tibble with one row per sample containing:

- sample_id:

  Sample identifier

- Mn:

  Number-average molecular weight

- Mw:

  Weight-average molecular weight

- Mz:

  Z-average molecular weight

- dispersity:

  Polydispersity index (Mw/Mn)

- purity_hmws:

  Percent high MW species (if available)

- purity_monomer:

  Percent monomer (if available)

- purity_lmws:

  Percent low MW species (if available)

## Details

This function creates a publication-ready summary table of SEC results.
It automatically detects and includes available metrics.

**Typical Summary Metrics:**

- Molecular weight averages: Mn, Mw, Mz

- Dispersity (PDI): Mw/Mn

- Purity metrics: %HMWS, %Monomer, %LMWS

- MW fractions: % above/below cutoffs

- Recovery: % mass balance

## See also

Other sec-export:
[`measure_sec_compare()`](https://jameshwade.github.io/measure-sec/reference/measure_sec_compare.md),
[`measure_sec_slice_table()`](https://jameshwade.github.io/measure-sec/reference/measure_sec_slice_table.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Generate summary after SEC processing
summary_tbl <- measure_sec_summary_table(
  processed_data,
  sample_id = "sample_name"
)

# Print formatted table
print(summary_tbl)

# Export to Excel
writexl::write_xlsx(summary_tbl, "sec_summary.xlsx")
} # }
```
