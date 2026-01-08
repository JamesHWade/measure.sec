# Get Standard Type

Classifies a canonical standard name into its polymer type.

## Usage

``` r
get_standard_type(standard_name)
```

## Arguments

- standard_name:

  Character string of the standard name (canonical form).

## Value

Character string of the standard type: "PS", "PMMA", "PEG", "Pullulan",
"Generic", or "Other".

## See also

Other standards:
[`auto_detect_standards()`](https://jameshwade.github.io/measure-sec/reference/auto_detect_standards.md),
[`get_recognized_standards()`](https://jameshwade.github.io/measure-sec/reference/get_recognized_standards.md),
[`get_standard_names()`](https://jameshwade.github.io/measure-sec/reference/get_standard_names.md),
[`is_standard()`](https://jameshwade.github.io/measure-sec/reference/is_standard.md),
[`match_standards()`](https://jameshwade.github.io/measure-sec/reference/match_standards.md),
[`normalize_standard_name()`](https://jameshwade.github.io/measure-sec/reference/normalize_standard_name.md)

## Examples

``` r
get_standard_type("PS A")
#> [1] "PS"
get_standard_type("Standard 2")
#> [1] "Generic"
get_standard_type("PMMA C")
#> [1] "PMMA"
```
