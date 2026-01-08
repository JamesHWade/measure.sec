# Match Sample Names to Recognized Standards

Attempts to match sample names to recognized calibration standards using
pattern matching on normalized names.

## Usage

``` r
match_standards(sample_names, standards = NULL, return_all = FALSE)
```

## Arguments

- sample_names:

  Character vector of sample names to match.

- standards:

  Optional named list of standards with regex patterns. If `NULL`
  (default), uses
  [`get_recognized_standards()`](https://jameshwade.github.io/measure-sec/reference/get_recognized_standards.md).

- return_all:

  Logical. If `TRUE`, returns all matches including non-matches as NA.
  If `FALSE` (default), returns only matched entries.

## Value

A tibble with columns:

- `sample_name`: Original sample name

- `normalized_name`: Cleaned/normalized name

- `matched_standard`: Canonical standard name (or NA if no match)

- `standard_type`: Type of standard (PS, PMMA, etc.)

## See also

Other standards:
[`auto_detect_standards()`](https://jameshwade.github.io/measure-sec/reference/auto_detect_standards.md),
[`get_recognized_standards()`](https://jameshwade.github.io/measure-sec/reference/get_recognized_standards.md),
[`get_standard_names()`](https://jameshwade.github.io/measure-sec/reference/get_standard_names.md),
[`get_standard_type()`](https://jameshwade.github.io/measure-sec/reference/get_standard_type.md),
[`is_standard()`](https://jameshwade.github.io/measure-sec/reference/is_standard.md),
[`normalize_standard_name()`](https://jameshwade.github.io/measure-sec/reference/normalize_standard_name.md)

## Examples

``` r
samples <- c("20231215_PS_A", "Unknown Sample", "Std B", "PMMA-C")
match_standards(samples)
#> # A tibble: 2 × 4
#>   sample_name   normalized_name matched_standard standard_type
#>   <chr>         <chr>           <chr>            <chr>        
#> 1 20231215_PS_A ps a            PS A             PS           
#> 2 PMMA-C        pmma c          PMMA C           PMMA         

# Return all samples including non-matches
match_standards(samples, return_all = TRUE)
#> # A tibble: 4 × 4
#>   sample_name    normalized_name matched_standard standard_type
#>   <chr>          <chr>           <chr>            <chr>        
#> 1 20231215_PS_A  ps a            PS A             PS           
#> 2 Unknown Sample unknown sample  NA               NA           
#> 3 Std B          standard b      NA               NA           
#> 4 PMMA-C         pmma c          PMMA C           PMMA         
```
