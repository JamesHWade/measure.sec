# Get Standard Names

Returns a character vector of all recognized standard names.

## Usage

``` r
get_standard_names()
```

## Value

Character vector of canonical standard names.

## See also

Other standards:
[`auto_detect_standards()`](https://jameshwade.github.io/measure-sec/reference/auto_detect_standards.md),
[`get_recognized_standards()`](https://jameshwade.github.io/measure-sec/reference/get_recognized_standards.md),
[`get_standard_type()`](https://jameshwade.github.io/measure-sec/reference/get_standard_type.md),
[`is_standard()`](https://jameshwade.github.io/measure-sec/reference/is_standard.md),
[`match_standards()`](https://jameshwade.github.io/measure-sec/reference/match_standards.md),
[`normalize_standard_name()`](https://jameshwade.github.io/measure-sec/reference/normalize_standard_name.md)

## Examples

``` r
get_standard_names()
#>  [1] "PS 1683"    "PS A"       "PS B"       "PS C"       "PS D"      
#>  [6] "PMMA A"     "PMMA B"     "PMMA C"     "PMMA D"     "PEG A"     
#> [11] "PEG B"      "PEG C"      "PEG D"      "Pullulan A" "Pullulan B"
#> [16] "Pullulan C" "Pullulan D" "Standard 1" "Standard 2" "Standard 3"
#> [21] "Standard 4"
```
