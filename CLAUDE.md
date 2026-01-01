# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Package Overview

`measure.sec` is an R package that extends the `measure` package with preprocessing steps specific to Size Exclusion Chromatography (SEC) and Gel Permeation Chromatography (GPC). It provides a recipes-style interface for SEC/GPC data processing, calibration, and molecular weight calculations.

This is a "technique pack" that registers with the `measure` package registry system.

## Development Commands

```r
# Generate documentation from roxygen2 comments
devtools::document()

# Run all tests
devtools::test()

# Run a single test file (filter matches test-*.R files)
devtools::test(filter = "conventional")  # runs test-step-conventional-cal.R

# Full CRAN-like package check
devtools::check()

# Test coverage
covr::package_coverage()

# Rebuild README.md from README.Rmd
devtools::build_readme()

# Build pkgdown site
pkgdown::build_site()
```

## Architecture

### Recipe Step Pattern

Every preprocessing step follows the recipes framework with three required methods:

```r
# Constructor function
step_sec_*() → calls recipes::add_step()

# Three S3 methods required:
prep.step_sec_*()   # Learn parameters from training data
bake.step_sec_*()   # Apply transformation to new data
tidy.step_sec_*()   # Return step parameters as tibble
```

Step implementations are in `R/` with naming convention: `step_sec_{operation}.R`

### Internal Data Format

The package uses the `measure` package's S3 class hierarchy:

- **`measure_tbl`**: Single measurement (tibble with `location` and `value` columns)
- **`measure_list`**: Collection of measurements stored as a list column

Key helper functions from `measure`:
- `find_measure_cols()`: Detect measure columns by class
- `check_for_measure()`: Validation utilities
- `new_measure_tbl()` / `new_measure_list()`: Constructors

### File Organization

**Detector Processing:**
- `step_sec_ri.R`: RI detector with dn/dc conversion
- `step_sec_uv.R`: UV detector with extinction coefficient
- `step_sec_mals.R`: Multi-angle light scattering
- `step_sec_viscometer.R`: Viscometer detector
- `step_sec_detector_delay.R`: Inter-detector delay correction

**Calibration:**
- `step_sec_conventional_cal.R`: Narrow standard calibration (polynomial fits)
- `step_sec_universal_cal.R`: Universal calibration with Mark-Houwink parameters

**Molecular Weight Calculations:**
- `step_sec_mw_averages.R`: Mn, Mw, Mz, dispersity
- `step_sec_mw_fractions.R`: MW fraction analysis
- `step_sec_mw_distribution.R`: MW distribution curves

**Signal Processing:**
- `step_sec_baseline.R`: SEC-optimized baseline correction
- `step_sec_concentration.R`: Signal to concentration conversion

**Analysis:**
- `step_sec_aggregates.R`: Aggregate/fragment analysis
- `step_sec_composition.R`: Copolymer composition
- `step_sec_uv_ri_ratio.R`: UV/RI ratio for composition
- `step_sec_intrinsic_visc.R`: Intrinsic viscosity calculation

**Supporting Files:**
- `qc-functions.R`: QC metrics (plate count, asymmetry, resolution)
- `polymer-analysis.R`: Polymer characterization (branching, M-H parameters)
- `export-functions.R`: Summary table and slice table exports
- `data.R`: Dataset documentation

## Testing

Tests use testthat edition 3. Test files mirror the source files:
- `tests/testthat/test-step-*.R`: Step-specific tests
- `tests/testthat/test-edge-cases.R`: Edge case coverage
- `tests/testthat/test-qc-functions.R`: QC function tests

## Code Style

- Tidyverse style guide
- Roxygen2 with Markdown syntax for documentation
- All exported functions need `@export` tag
- S3 methods should include `print()` and `tidy()` methods

## Linting & Formatting

The project uses two complementary code quality tools:

```bash
# Linter - catches logic/efficiency issues
jarl check .        # Check for issues
jarl check . --fix  # Auto-fix issues

# Formatter - enforces consistent code style
air format .        # Format all R files
```

**jarl** catches:
- `vector_logic`: Using `|` instead of `||` in `if()` statements
- Other potential bugs and inefficiencies

**air** enforces:
- 2-space indentation in function signatures
- One argument per line for long function calls
- Closing parenthesis on own line for multi-line constructs
- Consistent line length limits

Run both before committing to pass CI checks.

## PR Workflow

Before creating a pull request, run the following checks:

```bash
# 1. Format all R files (including tests)
air format .

# 2. Check for linting issues and auto-fix
jarl check . --fix

# 3. Generate/update documentation
R -e 'devtools::document()'

# 4. Run full package check (should pass with 0 errors, 0 warnings)
R -e 'devtools::check()'

# 5. Build pkgdown site to verify documentation renders correctly
R -e 'pkgdown::build_site()'
```

When adding new exported functions, ensure they are included in `_pkgdown.yml` under the appropriate reference section.

## Technique Pack Registration

The package registers with `measure` in `R/zzz.R`:

```r
.onLoad <- function(libname, pkgname) {

measure::register_measure_pack(
    pack_name = pkgname,
    technique = "SEC/GPC",
    description = "Size Exclusion / Gel Permeation Chromatography"
  )
  # ... register individual steps
}
```

Users can discover available steps with:
```r
measure::measure_steps(techniques = "SEC/GPC")
```

## Dependencies

- **measure**: Core package providing `measure_tbl`/`measure_list` classes and registry
- **recipes**: tidymodels preprocessing framework
- **tibble**, **dplyr**, **purrr**, **rlang**: tidyverse infrastructure
- **cli**: User-friendly messages and warnings

## Calibration Quality Metrics

The `step_sec_conventional_cal()` provides comprehensive diagnostics:

- **Per-standard**: residuals, % deviation, prediction intervals
- **Overall**: R², adjusted R², RMSE, residual standard error
- Access via `tidy(prepped_recipe, number = step_number)`
- Nested `standard_results` tibble contains per-standard diagnostics
