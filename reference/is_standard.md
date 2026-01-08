# Check if Sample Names Contain Standards

Quick check to determine if any sample names match recognized standards.

## Usage

``` r
is_standard(sample_names)
```

## Arguments

- sample_names:

  Character vector of sample names.

## Value

Logical vector indicating whether each name matches a standard.

## See also

Other standards:
[`auto_detect_standards()`](https://jameshwade.github.io/measure-sec/reference/auto_detect_standards.md),
[`get_recognized_standards()`](https://jameshwade.github.io/measure-sec/reference/get_recognized_standards.md),
[`get_standard_names()`](https://jameshwade.github.io/measure-sec/reference/get_standard_names.md),
[`get_standard_type()`](https://jameshwade.github.io/measure-sec/reference/get_standard_type.md),
[`match_standards()`](https://jameshwade.github.io/measure-sec/reference/match_standards.md),
[`normalize_standard_name()`](https://jameshwade.github.io/measure-sec/reference/normalize_standard_name.md)

## Examples

``` r
is_standard(c("PS A", "Unknown", "PMMA-B"))
#> [1]  TRUE FALSE  TRUE
```
