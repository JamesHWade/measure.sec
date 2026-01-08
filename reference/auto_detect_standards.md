# Auto-Detect Standards in a Data Frame

Scans a data frame for columns that might contain standard names and
attempts to match them to recognized standards.

## Usage

``` r
auto_detect_standards(data, name_col = NULL, add_columns = TRUE)
```

## Arguments

- data:

  A data frame containing sample information.

- name_col:

  Name of the column containing sample names. If `NULL`, attempts to
  auto-detect from common column names.

- add_columns:

  Logical. If `TRUE`, adds `matched_standard` and `standard_type`
  columns to the data frame. If `FALSE`, returns only the matching
  results.

## Value

If `add_columns = TRUE`, the original data frame with added columns. If
`FALSE`, a tibble with matching results.

## See also

Other standards:
[`get_recognized_standards()`](https://jameshwade.github.io/measure-sec/reference/get_recognized_standards.md),
[`get_standard_names()`](https://jameshwade.github.io/measure-sec/reference/get_standard_names.md),
[`get_standard_type()`](https://jameshwade.github.io/measure-sec/reference/get_standard_type.md),
[`is_standard()`](https://jameshwade.github.io/measure-sec/reference/is_standard.md),
[`match_standards()`](https://jameshwade.github.io/measure-sec/reference/match_standards.md),
[`normalize_standard_name()`](https://jameshwade.github.io/measure-sec/reference/normalize_standard_name.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Auto-detect and add standard columns
data <- data.frame(
  sample_name = c("PS A", "Unknown", "PMMA B"),
  signal = c(100, 200, 150)
)
auto_detect_standards(data)
} # }
```
