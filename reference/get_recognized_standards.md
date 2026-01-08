# Get Recognized SEC Calibration Standards

Returns a list of recognized calibration standards with their matching
patterns. This is the central registry for standard name recognition.

## Usage

``` r
get_recognized_standards()
```

## Value

A named list where names are canonical standard names and values are
regex patterns for matching normalized sample names.

## Details

Standard types supported:

- **PS (Polystyrene)**: PS A, PS B, PS C, PS D, PS 1683

- **PMMA (Polymethyl methacrylate)**: PMMA A, PMMA B, PMMA C, PMMA D

- **PEG/PEO (Polyethylene glycol/oxide)**: PEG A, PEG B, PEG C, PEG D

- **Pullulan**: Pullulan A, Pullulan B, Pullulan C, Pullulan D

- **Generic**: Standard 1-4 (site-specific, no default MW values)

Patterns assume input has been normalized via
[`normalize_standard_name()`](https://jameshwade.github.io/measure-sec/reference/normalize_standard_name.md)
(lowercase, "std" converted to "standard", separators normalized).

## See also

Other standards:
[`auto_detect_standards()`](https://jameshwade.github.io/measure-sec/reference/auto_detect_standards.md),
[`get_standard_names()`](https://jameshwade.github.io/measure-sec/reference/get_standard_names.md),
[`get_standard_type()`](https://jameshwade.github.io/measure-sec/reference/get_standard_type.md),
[`is_standard()`](https://jameshwade.github.io/measure-sec/reference/is_standard.md),
[`match_standards()`](https://jameshwade.github.io/measure-sec/reference/match_standards.md),
[`normalize_standard_name()`](https://jameshwade.github.io/measure-sec/reference/normalize_standard_name.md)

## Examples

``` r
standards <- get_recognized_standards()
names(standards)
#>  [1] "PS 1683"    "PS A"       "PS B"       "PS C"       "PS D"      
#>  [6] "PMMA A"     "PMMA B"     "PMMA C"     "PMMA D"     "PEG A"     
#> [11] "PEG B"      "PEG C"      "PEG D"      "Pullulan A" "Pullulan B"
#> [16] "Pullulan C" "Pullulan D" "Standard 1" "Standard 2" "Standard 3"
#> [21] "Standard 4"
```
