# ==============================================================================
# Dataset Documentation
# ==============================================================================

#' Multi-Detector SEC Data
#'
#' A synthetic dataset containing multi-detector Size Exclusion Chromatography
#' (SEC) data for 12 polymer samples with realistic signal characteristics.
#'
#' @format A tibble with 24,012 rows and 11 columns:
#' \describe{
#'   \item{sample_id}{Character. Unique sample identifier (e.g., "PS-10K", "PMMA-Low")}
#'   \item{sample_type}{Character. Either "standard" (narrow dispersity calibrants) or "sample"}
#'   \item{polymer_type}{Character. Polymer type: "polystyrene", "pmma", "peg", or "copolymer"}
#'   \item{elution_time}{Numeric. Elution time in minutes (5-25 min range)}
#'   \item{ri_signal}{Numeric. Refractive index detector signal (reference detector)}
#'   \item{uv_signal}{Numeric. UV detector signal at 280 nm}
#'   \item{mals_signal}{Numeric. Multi-angle light scattering detector signal}
#'   \item{known_mw}{Numeric. True weight-average molecular weight (Mw) in g/mol}
#'   \item{known_dispersity}{Numeric. True dispersity (Mw/Mn)}
#'   \item{dn_dc}{Numeric. Refractive index increment in mL/g}
#'   \item{extinction_coef}{Numeric. UV extinction coefficient in mL/(mg*cm)}
#' }
#'
#' @details
#' The dataset includes realistic features commonly encountered in SEC analysis:
#'
#' **Sample Composition:**
#' \itemize{
#'   \item 5 polystyrene standards (1K to 500K MW, narrow dispersity ~1.01-1.05)
#'   \item 3 PMMA samples (25K to 200K MW, broader dispersity 1.8-2.2)
#'   \item 2 PEG samples (5K and 20K MW, low dispersity ~1.1-1.15)
#'   \item 2 copolymer samples (40K and 80K MW, intermediate dispersity 1.5-1.7)
#' }
#'
#' **Multi-Detector Features:**
#' \itemize{
#'   \item Inter-detector volume delays: UV is 0.05 mL before RI, MALS is 0.15 mL after RI
#'   \item Different detector responses based on polymer chemistry
#'   \item PEG has no UV response (extinction_coef = 0)
#'   \item MALS signal scales with MW for absolute MW determination
#' }
#'
#' **Signal Characteristics:**
#' \itemize{
#'   \item Gaussian noise appropriate for each detector
#'   \item Slight baseline drift
#'   \item Log-normal peak shapes with tailing
#' }
#'
#' **Typical SEC Workflow:**
#' 1
#' . Convert to measure format with \code{\link[measure]{step_measure_input_long}}
#' 2. Correct inter-detector delays with \code{\link{step_sec_detector_delay}}
#' 3. Apply baseline correction with \code{\link{step_sec_baseline}}
#' 4. Process detectors with \code{\link{step_sec_ri}} or \code{\link{step_sec_uv}}
#' 5. Convert to concentration with \code{\link{step_sec_concentration}}
#' 6. Calculate MW averages with \code{\link{step_sec_mw_averages}}
#'
#' @family sec-data
#'
#' @examples
#' \dontrun{
#' library(recipes)
#' library(measure)
#' library(measure.sec)
#'
#' # Load the dataset
#' data(sec_triple_detect)
#'
#' # View sample distribution
#' table(sec_triple_detect$polymer_type)
#'
#' # Plot RI chromatograms for polystyrene standards
#' library(ggplot2)
#' sec_triple_detect |>
#'   dplyr::filter(polymer_type == "polystyrene") |>
#'   ggplot(aes(elution_time, ri_signal, color = sample_id)) +
#'   geom_line() +
#'   labs(x = "Elution Time (min)", y = "RI Signal", color = "Sample")
#'
#' # Process with SEC recipe
#' rec <- recipe(~ ., data = sec_triple_detect) |>
#'   step_measure_input_long(ri_signal, location = vars(elution_time), col_name = "ri") |>
#'   step_sec_baseline(measures = "ri") |>
#'   step_sec_ri(dn_dc_column = "dn_dc") |>
#'   prep()
#' }
#'
#' @source Synthetic data generated for package testing and examples.
"sec_triple_detect"


# ==============================================================================
# Calibration Standards Datasets
# ==============================================================================

#' SEC Calibration Standards - Complete Dataset
#'
#' Comprehensive calibration standard data for SEC/GPC analysis, containing
#' polystyrene and PMMA narrow standards with certificate values, retention
#' data, and Mark-Houwink parameters. Based on commercial standard kits.
#'
#' @format A tibble with 26 rows and 19 columns:
#' \describe{
#'   \item{standard_name}{Character. Standard identifier (e.g., "PS-67500", "PMMA-30300")}
#'   \item{polymer_type}{Character. Either "polystyrene" or "pmma"}
#'   \item{kit_name}{Character. Commercial kit name for reference}
#'   \item{mp}{Numeric. Peak molecular weight in Da (most commonly used for calibration)}
#'   \item{mn}{Numeric. Number-average molecular weight in Da}
#'   \item{mw}{Numeric. Weight-average molecular weight in Da}
#'   \item{dispersity}{Numeric. Polydispersity index (Mw/Mn), typically 1.02-1.05}
#'   \item{mp_uncertainty}{Numeric. Relative uncertainty in Mp (e.g., 0.05 = 5%)}
#'   \item{log_mp}{Numeric. log10(Mp) for calibration curve fitting}
#'   \item{log_mw}{Numeric. log10(Mw)}
#'   \item{log_mn}{Numeric. log10(Mn)}
#'   \item{retention_time}{Numeric. Peak retention time in minutes}
#'   \item{retention_volume}{Numeric. Peak retention volume in mL}
#'   \item{k_value}{Numeric. Mark-Houwink K constant in mL/g}
#'   \item{a_value}{Numeric. Mark-Houwink exponent (alpha)}
#'   \item{intrinsic_viscosity}{Numeric. Intrinsic viscosity in mL/g}
#'   \item{log_hydrodynamic_vol}{Numeric. log10(M * \[eta\]) for universal calibration}
#'   \item{dn_dc}{Numeric. Refractive index increment in mL/g}
#'   \item{notes}{Character. Special notes (e.g., near exclusion limit)}
#' }
#'
#' @details
#' This dataset enables several key calibration workflows:
#'
#' **Conventional Calibration:**
#' Use `retention_time` (or `retention_volume`) and `log_mp` with
#' \code{\link{step_sec_conventional_cal}} to build a calibration curve.
#'
#' **Universal Calibration:**
#' Use `log_hydrodynamic_vol` vs retention for polymer-independent calibration.
#' The Mark-Houwink parameters (`k_value`, `a_value`) enable conversion between
#' polymers.
#'
#' **Quality Assessment:**
#' The `mp_uncertainty` values (from typical certificates) enable uncertainty
#' propagation through the calibration. Dispersity values confirm standards
#' are suitably narrow.
#'
#' **Polymer Types:**
#' \itemize{
#'   \item 16 polystyrene standards (162 Da to 3,150,000 Da)
#'   \item 10 PMMA standards (602 Da to 1,190,000 Da)
#' }
#'
#' **Mark-Houwink Parameters (THF, 35°C):**
#' \itemize{
#'   \item Polystyrene: K = 0.000141 mL/g, a = 0.700
#'   \item PMMA: K = 0.000128 mL/g, a = 0.690
#' }
#'
#' @seealso
#' \code{\link{sec_ps_standards}} for polystyrene-only subset
#' \code{\link{sec_pmma_standards}} for PMMA-only subset
#' \code{\link{step_sec_conventional_cal}} for conventional calibration
#' \code{\link{step_sec_universal_cal}} for universal calibration
#'
#' @family sec-data
#' @family sec-calibration
#'
#' @examples
#' library(dplyr)
#' data(sec_calibration_standards)
#'
#' # View the polystyrene standards
#' sec_calibration_standards |>
#'   filter(polymer_type == "polystyrene") |>
#'   select(standard_name, mp, log_mp, retention_time)
#'
#' # Compare PS and PMMA at similar MW - PMMA elutes later (smaller Rh)
#' sec_calibration_standards |>
#'   filter(mp > 60000 & mp < 80000) |>
#'   select(standard_name, polymer_type, mp, retention_time)
#'
#' # Prepare standards for conventional calibration
#' ps_cal <- sec_calibration_standards |>
#'   filter(polymer_type == "polystyrene") |>
#'   select(retention = retention_time, log_mw = log_mp)
#'
#' @source Synthetic data based on typical commercial narrow standards
#'   (e.g., Agilent EasiVial, PSS ReadyCal) with realistic retention times
#'   for PLgel Mixed-C columns in THF at 1.0 mL/min.
"sec_calibration_standards"


#' SEC Polystyrene Calibration Standards
#'
#' Polystyrene narrow molecular weight standards for SEC/GPC conventional
#' calibration. A convenient subset of \code{\link{sec_calibration_standards}}
#' containing only polystyrene standards.
#'
#' @format A tibble with 16 rows and 12 columns:
#' \describe{
#'   \item{standard_name}{Character. Standard identifier (e.g., "PS-67500")}
#'   \item{mp}{Numeric. Peak molecular weight in Da}
#'   \item{log_mp}{Numeric. log10(Mp) for calibration curve fitting}
#'   \item{retention_time}{Numeric. Peak retention time in minutes}
#'   \item{retention_volume}{Numeric. Peak retention volume in mL}
#'   \item{mn}{Numeric. Number-average molecular weight in Da}
#'   \item{mw}{Numeric. Weight-average molecular weight in Da}
#'   \item{dispersity}{Numeric. Polydispersity index (Mw/Mn)}
#'   \item{mp_uncertainty}{Numeric. Relative uncertainty in Mp}
#'   \item{k_value}{Numeric. Mark-Houwink K constant (0.000141 mL/g)}
#'   \item{a_value}{Numeric. Mark-Houwink exponent (0.700)}
#'   \item{dn_dc}{Numeric. Refractive index increment (0.185 mL/g)}
#' }
#'
#' @details
#' Polystyrene is the most widely used SEC calibration standard due to its:
#' \itemize{
#'   \item Availability in narrow dispersity grades across wide MW range
#'   \item Well-characterized Mark-Houwink parameters in common solvents
#'   \item Strong UV absorption for dual detection
#'   \item Good solubility and stability
#' }
#'
#' The 16 standards span 162 Da to 3,150,000 Da, covering typical analytical
#' SEC columns. Standards are pre-sorted by descending molecular weight
#' (elution order).
#'
#' **Usage with step_sec_conventional_cal:**
#' \preformatted{
#' library(dplyr)
#' standards <- sec_ps_standards |>
#'   select(retention = retention_time, log_mw = log_mp)
#'
#' recipe(~., data = my_data) |>
#'   step_sec_conventional_cal(standards = standards, fit_type = "cubic")
#' }
#'
#' @seealso
#' \code{\link{sec_calibration_standards}} for full dataset with PMMA
#' \code{\link{sec_pmma_standards}} for PMMA standards
#' \code{\link{step_sec_conventional_cal}} for calibration step
#'
#' @family sec-data
#' @family sec-calibration
#'
#' @examples
#' data(sec_ps_standards)
#'
#' # Quick look at the calibration range
#' range(sec_ps_standards$mp)
#' range(sec_ps_standards$retention_time)
#'
#' # Plot calibration curve
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   library(ggplot2)
#'   ggplot(sec_ps_standards, aes(retention_time, log_mp)) +
#'     geom_point(size = 3) +
#'     geom_smooth(method = "lm", formula = y ~ poly(x, 3), se = FALSE) +
#'     labs(
#'       x = "Retention Time (min)",
#'       y = expression(log[10](M[p])),
#'       title = "PS Calibration Curve (Cubic Fit)"
#'     ) +
#'     theme_minimal()
#' }
#'
#' @source Synthetic data based on typical commercial narrow PS standards.
"sec_ps_standards"


#' SEC PMMA Calibration Standards
#'
#' Poly(methyl methacrylate) narrow molecular weight standards for SEC/GPC
#' calibration. A convenient subset of \code{\link{sec_calibration_standards}}
#' containing only PMMA standards.
#'
#' @format A tibble with 10 rows and 12 columns:
#' \describe{
#'   \item{standard_name}{Character. Standard identifier (e.g., "PMMA-67700")}
#'   \item{mp}{Numeric. Peak molecular weight in Da}
#'   \item{log_mp}{Numeric. log10(Mp) for calibration curve fitting}
#'   \item{retention_time}{Numeric. Peak retention time in minutes}
#'   \item{retention_volume}{Numeric. Peak retention volume in mL}
#'   \item{mn}{Numeric. Number-average molecular weight in Da}
#'   \item{mw}{Numeric. Weight-average molecular weight in Da}
#'   \item{dispersity}{Numeric. Polydispersity index (Mw/Mn)}
#'   \item{mp_uncertainty}{Numeric. Relative uncertainty in Mp}
#'   \item{k_value}{Numeric. Mark-Houwink K constant (0.000128 mL/g)}
#'   \item{a_value}{Numeric. Mark-Houwink exponent (0.690)}
#'   \item{dn_dc}{Numeric. Refractive index increment (0.084 mL/g)}
#' }
#'
#' @details
#' PMMA standards are useful for:
#' \itemize{
#'   \item Calibrating for PMMA or acrylate samples
#'   \item Validating universal calibration by comparing PS vs PMMA curves
#'   \item Demonstrating polymer-specific hydrodynamic volume differences
#' }
#'
#' At equivalent molecular weight, PMMA has a smaller hydrodynamic volume than
#' PS in THF, so PMMA standards elute later than PS of the same MW. This
#' demonstrates why conventional calibration is polymer-specific.
#'
#' **Universal Calibration Validation:**
#' When plotted as log(M * \[eta\]) vs retention time, PS and PMMA should fall
#' on the same curve, confirming universal calibration is valid for the column.
#'
#' @seealso
#' \code{\link{sec_calibration_standards}} for full dataset with PS
#' \code{\link{sec_ps_standards}} for polystyrene standards
#' \code{\link{step_sec_universal_cal}} for universal calibration
#'
#' @family sec-data
#' @family sec-calibration
#'
#' @examples
#' data(sec_pmma_standards)
#' data(sec_ps_standards)
#'
#' # Compare PS and PMMA calibration curves
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   library(ggplot2)
#'   library(dplyr)
#'
#'   bind_rows(
#'     sec_ps_standards |> mutate(polymer = "PS"),
#'     sec_pmma_standards |> mutate(polymer = "PMMA")
#'   ) |>
#'     ggplot(aes(retention_time, log_mp, color = polymer)) +
#'     geom_point(size = 3) +
#'     geom_smooth(method = "lm", formula = y ~ poly(x, 3), se = FALSE) +
#'     labs(
#'       x = "Retention Time (min)",
#'       y = expression(log[10](M[p])),
#'       title = "PS vs PMMA Calibration Curves",
#'       subtitle = "PMMA elutes later at same MW (smaller hydrodynamic volume)"
#'     ) +
#'     theme_minimal()
#' }
#'
#' @source Synthetic data based on typical commercial narrow PMMA standards.
"sec_pmma_standards"
