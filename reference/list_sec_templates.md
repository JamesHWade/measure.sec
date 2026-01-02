# List Available SEC Report Templates

Lists all available report templates with descriptions.

## Usage

``` r
list_sec_templates()
```

## Value

A tibble with template names and descriptions.

## Examples

``` r
list_sec_templates()
#> # A tibble: 3 × 3
#>   template description                                         formats        
#>   <chr>    <chr>                                               <chr>          
#> 1 standard Summary table, chromatogram, and MWD plot           html, pdf, docx
#> 2 detailed All plots, multi-detector view, optional slice data html, pdf, docx
#> 3 qc       System suitability with pass/fail metrics           html, pdf, docx
```
