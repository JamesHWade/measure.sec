# Normalize Standard Names for Matching

Normalizes sample and standard names to improve matching by:

- Converting to lowercase

- Replacing separators (-, \_) with spaces

- Removing date prefixes (YYYYMMDD, YYYY-MM-DD, etc.)

- Removing common lab/vendor identifiers

- Normalizing "std" variations to "standard"

- Cleaning up whitespace

## Usage

``` r
normalize_standard_name(name)
```

## Arguments

- name:

  Character vector of names to normalize.

## Value

Character vector of normalized names.

## See also

Other standards:
[`auto_detect_standards()`](https://jameshwade.github.io/measure-sec/reference/auto_detect_standards.md),
[`get_recognized_standards()`](https://jameshwade.github.io/measure-sec/reference/get_recognized_standards.md),
[`get_standard_names()`](https://jameshwade.github.io/measure-sec/reference/get_standard_names.md),
[`get_standard_type()`](https://jameshwade.github.io/measure-sec/reference/get_standard_type.md),
[`is_standard()`](https://jameshwade.github.io/measure-sec/reference/is_standard.md),
[`match_standards()`](https://jameshwade.github.io/measure-sec/reference/match_standards.md)

## Examples

``` r
normalize_standard_name("20231215_PS_A")
#> [1] "ps a"
normalize_standard_name("Std B")
#> [1] "standard b"
normalize_standard_name("Phillips PS-C")
#> [1] "ps c"
normalize_standard_name(c("DOW_PMMA_A", "2023-01-15 StdD"))
#> [1] "pmma a"     "standard d"
```
