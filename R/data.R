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

# ==============================================================================
# Additional Example Datasets
# ==============================================================================

#' Copolymer SEC Data for Composition Analysis
#'
#' A synthetic dataset containing SEC chromatograms of styrene-acrylate
#' copolymers with varying compositions, designed for demonstrating UV/RI
#' ratio analysis for composition determination.
#'
#' @format A tibble with 4,206 rows and 8 columns:
#' \describe{
#'   \item{sample_id}{Character. Sample identifier (e.g., "Copoly-20S")}
#'   \item{elution_time}{Numeric. Elution time in minutes}
#'   \item{ri_signal}{Numeric. Refractive index detector signal}
#'   \item{uv_254_signal}{Numeric. UV detector signal at 254 nm}
#'   \item{styrene_fraction}{Numeric. Styrene content (0-1)}
#'   \item{mw}{Numeric. Weight-average molecular weight in Da}
#'   \item{dispersity}{Numeric. Polydispersity index (Mw/Mn)}
#'   \item{description}{Character. Sample description}
#' }
#'
#' @details
#' The dataset includes 6 samples spanning the full composition range:
#' \itemize{
#'   \item Pure polyacrylate (0% styrene) - no UV absorption
#'   \item 20%, 40%, 60%, 80% styrene copolymers
#'   \item Pure polystyrene (100% styrene) - strong UV absorption
#' }
#'
#' **UV/RI Ratio Analysis:**
#' The UV signal at 254 nm is selective for styrene units, while the RI signal
#' responds to total mass. The UV/RI ratio across the chromatogram reveals
#' composition as a function of molecular weight, enabling detection of
#' compositional drift.
#'
#' **Typical Workflow:**
#' 1. Load data and convert to measure format
#' 2. Apply \code{\link{step_sec_uv_ri_ratio}} to calculate ratios
#' 3. Calibrate ratio to composition using homopolymer standards
#' 4. Plot composition vs molecular weight
#'
#' @family sec-data
#'
#' @examples
#' data(sec_copolymer)
#'
#' # View composition range
#' unique(sec_copolymer[, c("sample_id", "styrene_fraction")])
#'
#' # Plot RI vs UV for different compositions
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   library(ggplot2)
#'   ggplot(sec_copolymer, aes(elution_time)) +
#'     geom_line(aes(y = ri_signal, color = "RI")) +
#'     geom_line(aes(y = uv_254_signal, color = "UV 254nm")) +
#'     facet_wrap(~sample_id) +
#'     labs(x = "Elution Time (min)", y = "Signal", color = "Detector") +
#'     theme_minimal()
#' }
#'
#' @source Synthetic data generated for package testing and examples.
"sec_copolymer"

#' Protein SEC Data with Aggregates
#'
#' A synthetic dataset containing SEC chromatograms of a monoclonal antibody
#' (mAb) under various stress conditions, showing monomer, dimer, higher-order
#' aggregates, and fragments.
#'
#' @format A tibble with 6,255 rows and 9 columns:
#' \describe{
#'   \item{sample_id}{Character. Sample identifier (e.g., "mAb-Reference")}
#'   \item{elution_time}{Numeric. Elution time in minutes}
#'   \item{uv_280_signal}{Numeric. UV detector signal at 280 nm}
#'   \item{uv_214_signal}{Numeric. UV detector signal at 214 nm}
#'   \item{description}{Character. Sample treatment description}
#'   \item{monomer_pct}{Numeric. Known monomer percentage}
#'   \item{dimer_pct}{Numeric. Known dimer percentage}
#'   \item{hmw_pct}{Numeric. Known high molecular weight aggregate percentage}
#'   \item{fragment_pct}{Numeric. Known fragment percentage}
#' }
#'
#' @details
#' The dataset simulates a typical mAb (~150 kDa) under various conditions:
#' \itemize{
#'   \item Reference standard (>98% monomer)
#'   \item Heat-stressed samples (40°C for 1-2 weeks)
#'   \item Aged sample (12-month stability)
#'   \item Freeze-thaw stressed sample
#' }
#'
#' **Species Present:**
#' \itemize{
#'   \item High molecular weight (HMW) aggregates (~600 kDa)
#'   \item Dimer (~300 kDa)
#'   \item Monomer (~150 kDa)
#'   \item Fragments (~50 kDa, Fab-like)
#' }
#'
#' **Typical Workflow:**
#' 1. Load data and apply baseline correction
#' 2. Use \code{\link{step_sec_aggregates}} to identify and quantify species
#' 3. Calculate percent area for each species
#' 4. Compare to acceptance criteria
#'
#' @family sec-data
#'
#' @examples
#' data(sec_protein)
#'
#' # View sample conditions
#' unique(sec_protein[, c("sample_id", "description", "monomer_pct")])
#'
#' # Plot stressed vs reference
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   library(ggplot2)
#'   library(dplyr)
#'
#'   sec_protein |>
#'     filter(sample_id %in% c("mAb-Reference", "mAb-Stressed-2")) |>
#'     ggplot(aes(elution_time, uv_280_signal, color = sample_id)) +
#'     geom_line() +
#'     labs(
#'       x = "Elution Time (min)",
#'       y = "UV 280 nm Signal",
#'       title = "Protein SEC: Reference vs Stressed"
#'     ) +
#'     theme_minimal()
#' }
#'
#' @source Synthetic data generated for package testing and examples.
"sec_protein"

#' Branched Polymer SEC Data
#'
#' A synthetic dataset containing SEC chromatograms of linear, branched, and
#' star polymers, designed for demonstrating branching analysis using
#' multi-detector SEC.
#'
#' @format A tibble with 5,608 rows and 10 columns:
#' \describe{
#'   \item{sample_id}{Character. Sample identifier (e.g., "Linear-50K")}
#'   \item{elution_time}{Numeric. Elution time in minutes}
#'   \item{ri_signal}{Numeric. Refractive index detector signal}
#'   \item{visc_signal}{Numeric. Viscometer detector signal}
#'   \item{mals_signal}{Numeric. Multi-angle light scattering signal}
#'   \item{topology}{Character. Polymer topology: "linear", "branched", or "star"}
#'   \item{mw}{Numeric. Weight-average molecular weight in Da}
#'   \item{branching_index}{Numeric. Branching index g' (1.0 for linear)}
#'   \item{intrinsic_visc}{Numeric. Intrinsic viscosity in mL/g}
#'   \item{rg}{Numeric. Radius of gyration in nm}
#' }
#'
#' @details
#' The dataset includes polyethylene-like samples with three topologies:
#' \itemize{
#'   \item Linear polymers (g' = 1.0) at 50K, 100K, 200K MW
#'   \item Branched polymers (g' = 0.55-0.75) at same MW range
#'   \item Star polymers (g' = 0.50-0.60) at 50K, 100K MW
#' }
#'
#' **Key Observation:**
#' At the same molecular weight, branched polymers have smaller hydrodynamic
#' volume and elute LATER than their linear counterparts. This is the basis
#' for branching analysis by comparing SEC retention to absolute MW from MALS.
#'
#' **Branching Index (g'):**
#' g' = \eqn{[\eta]_{branched}} / \eqn{[\eta]_{linear}} at same MW
#'
#' Values less than 1.0 indicate branching. Lower values mean more compact
#' (more branched) structures.
#'
#' **Typical Workflow:**
#' 1. Apply \code{\link{step_sec_mals}} for absolute MW
#' 2. Apply \code{\link{step_sec_viscometer}} for intrinsic viscosity
#' 3. Use \code{\link{measure_branching_index}} to calculate g'
#' 4. Plot Mark-Houwink relationship to compare topologies
#'
#' @family sec-data
#'
#' @examples
#' data(sec_branched)
#'
#' # Compare topologies
#' unique(sec_branched[, c("sample_id", "topology", "mw", "branching_index")])
#'
#' # Plot linear vs branched at same MW
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   library(ggplot2)
#'   library(dplyr)
#'
#'   sec_branched |>
#'     filter(mw == 100000) |>
#'     ggplot(aes(elution_time, ri_signal, color = topology)) +
#'     geom_line() +
#'     labs(
#'       x = "Elution Time (min)",
#'       y = "RI Signal",
#'       title = "100K MW: Linear vs Branched vs Star",
#'       subtitle = "Branched polymers elute later (smaller Vh)"
#'     ) +
#'     theme_minimal()
#' }
#'
#' @source Synthetic data generated for package testing and examples.
"sec_branched"

#' System Suitability Test Data
#'
#' A synthetic dataset containing SEC system suitability test (SST) runs
#' showing column degradation over time, useful for QC testing and
#' demonstrating column performance metrics.
#'
#' @format A tibble with 10,005 rows and 7 columns:
#' \describe{
#'   \item{run_id}{Character. Run identifier (e.g., "SST-Day0")}
#'   \item{elution_time}{Numeric. Elution time in minutes}
#'   \item{ri_signal}{Numeric. Refractive index detector signal}
#'   \item{column_age_days}{Numeric. Column age in days}
#'   \item{expected_plate_count}{Numeric. Expected theoretical plate count}
#'   \item{expected_asymmetry}{Numeric. Expected peak asymmetry factor}
#'   \item{expected_resolution}{Numeric. Expected resolution between peaks}
#' }
#'
#' @details
#' The dataset simulates SST runs using two polystyrene standards (50K and
#' 100K MW) measured at different column ages (0, 30, 60, 90, 120 days).
#'
#' **Column Degradation Effects:**
#' As the column ages, typical degradation patterns include:
#' \itemize{
#'   \item Decreased plate count (broader peaks)
#'   \item Increased asymmetry (more tailing)
#'   \item Decreased resolution between peaks
#' }
#'
#' **QC Metrics Available:**
#' Use the QC functions to calculate and verify:
#' \itemize{
#'   \item \code{\link{measure_sec_plate_count}} - Theoretical plates
#'   \item \code{\link{measure_sec_asymmetry}} - Peak asymmetry
#'   \item \code{\link{measure_sec_resolution}} - Peak resolution
#'   \item \code{\link{measure_sec_suitability}} - Combined SST check
#' }
#'
#' **Typical Workflow:**
#' 1. Load daily SST run
#' 2. Calculate QC metrics
#' 3. Compare to specifications
#' 4. Track trends over time
#'
#' @family sec-data
#'
#' @examples
#' data(sec_system_suitability)
#'
#' # View column age progression
#' unique(sec_system_suitability[,
#'   c("run_id", "column_age_days", "expected_plate_count")])
#'
#' # Plot degradation over time
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   library(ggplot2)
#'
#'   ggplot(sec_system_suitability, aes(elution_time, ri_signal, color = run_id)) +
#'     geom_line() +
#'     labs(
#'       x = "Elution Time (min)",
#'       y = "RI Signal",
#'       title = "System Suitability: Column Degradation Over Time"
#'     ) +
#'     theme_minimal()
#' }
#'
#' @source Synthetic data generated for package testing and examples.
"sec_system_suitability"

# ==============================================================================
# Raw Chromatogram Datasets for Tutorials
# ==============================================================================

#' Raw SEC Calibration Standards
#'
#' A dataset containing raw SEC chromatograms of polystyrene narrow standards
#' with realistic noise, baseline drift, and injection artifacts. Designed to
#' mimic data exported directly from SEC instruments for tutorial purposes.
#'
#' @format A tibble with approximately 130,000 rows and 6 columns:
#' \describe{
#'   \item{standard_name}{Character. Standard identifier (e.g., "PS-67500")}
#'   \item{mp}{Numeric. Peak molecular weight in Da from certificate}
#'   \item{log_mp}{Numeric. log10(Mp) for calibration curve fitting}
#'   \item{dispersity}{Numeric. Polydispersity index from certificate}
#'   \item{time_min}{Numeric. Elution time in minutes}
#'   \item{ri_mv}{Numeric. RI detector signal in millivolts (raw, unprocessed)}
#' }
#'
#' @details
#' This dataset represents "raw" data as it would come from an SEC instrument,
#' before any processing. Key characteristics include:
#'
#' **Realistic Signal Features:**
#' \itemize{
#'   \item Gaussian noise on detector signal
#'   \item Slow baseline drift from temperature fluctuations
#'   \item Injection artifacts at start of run
#'   \item Peak tailing typical of SEC columns
#' }
#'
#' **Standards Included:**
#' 12 polystyrene narrow standards spanning 580 Da to 930,000 Da, covering
#' the typical analytical SEC range. Standards are based on commercial kit
#' values.
#'
#' **Typical Tutorial Workflow:**
#' 1. Load raw data and inspect for quality
#' 2. Apply baseline correction with \code{\link{step_sec_baseline}}
#' 3. Identify peak retention times
#' 4. Build calibration curve with \code{\link{step_sec_conventional_cal}}
#'
#' @seealso
#' \code{\link{sec_raw_unknowns}} for unknown samples to analyze
#' \code{\link{sec_ps_standards}} for pre-processed calibration data
#' \code{\link{step_sec_conventional_cal}} for calibration step
#'
#' @family sec-data
#' @family sec-raw
#'
#' @examples
#' data(sec_raw_standards)
#'
#' # View available standards
#' unique(sec_raw_standards[, c("standard_name", "mp", "dispersity")])
#'
#' # Plot a single standard (shows noise and baseline)
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   library(ggplot2)
#'   library(dplyr)
#'
#'   sec_raw_standards |>
#'     filter(standard_name == "PS-67500") |>
#'     ggplot(aes(time_min, ri_mv)) +
#'     geom_line() +
#'     labs(
#'       x = "Time (min)",
#'       y = "RI Signal (mV)",
#'       title = "Raw PS-67500 Standard",
#'       subtitle = "Note baseline drift and noise"
#'     ) +
#'     theme_minimal()
#' }
#'
#' @source Synthetic data generated for package tutorials.
"sec_raw_standards"

#' Raw SEC Unknown Samples with Known Molecular Weights
#'
#' A dataset containing raw SEC chromatograms of polymer samples with known
#' true molecular weight values for validation. Includes challenging cases
#' like bimodal distributions and high molecular weight aggregates.
#'
#' @format A tibble with approximately 65,000 rows and 8 columns:
#' \describe{
#'   \item{sample_id}{Character. Sample identifier}
#'   \item{description}{Character. Sample description and type}
#'   \item{true_mw}{Numeric. Known weight-average MW in Da (NA for bimodal)}
#'   \item{true_mn}{Numeric. Known number-average MW in Da}
#'   \item{true_mz}{Numeric. Known z-average MW in Da}
#'   \item{true_dispersity}{Numeric. Known dispersity (Mw/Mn)}
#'   \item{time_min}{Numeric. Elution time in minutes}
#'   \item{ri_mv}{Numeric. RI detector signal in millivolts (raw)}
#' }
#'
#' @details
#' This dataset provides unknown samples with known "true" MW values, enabling
#' validation of the complete SEC workflow from raw data to final results.
#'
#' **Sample Types:**
#' \itemize{
#'   \item \strong{Unknown-A}: Broad distribution (dispersity ~2.0)
#'   \item \strong{Unknown-B}: Medium dispersity (~1.3)
#'   \item \strong{Unknown-C}: Narrow distribution (~1.1) for accuracy check
#'   \item \strong{Unknown-Bimodal}: Two-peak mixture (50K + 200K)
#'   \item \strong{Unknown-HMW}: Very high MW (~1.5M) with aggregate shoulder
#'   \item \strong{Unknown-LMW}: Low MW oligomers (~3.5K)
#' }
#'
#' **Educational Value:**
#' Students can compare their calculated MW values against the true values
#' to validate their analysis workflow and understand sources of error.
#'
#' @seealso
#' \code{\link{sec_raw_standards}} for calibration standards
#' \code{\link{step_sec_mw_averages}} for MW calculation
#'
#' @family sec-data
#' @family sec-raw
#'
#' @examples
#' data(sec_raw_unknowns)
#'
#' # View sample information
#' unique(sec_raw_unknowns[,
#'   c("sample_id", "description", "true_mw", "true_dispersity")])
#'
#' # Plot the bimodal sample
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   library(ggplot2)
#'   library(dplyr)
#'
#'   sec_raw_unknowns |>
#'     filter(sample_id == "Unknown-Bimodal") |>
#'     ggplot(aes(time_min, ri_mv)) +
#'     geom_line() +
#'     labs(
#'       x = "Time (min)",
#'       y = "RI Signal (mV)",
#'       title = "Bimodal Distribution",
#'       subtitle = "Mixture of 50K and 200K components"
#'     ) +
#'     theme_minimal()
#' }
#'
#' @source Synthetic data generated for package tutorials.
"sec_raw_unknowns"

#' Raw Multi-Detector SEC Data
#'
#' A dataset containing raw SEC chromatograms with RI, UV, and MALS detector
#' signals that have NOT been corrected for inter-detector delay. Designed
#' for teaching triple detection workflows from raw data.
#'
#' @format A tibble with approximately 48,000 rows and 12 columns:
#' \describe{
#'   \item{sample_id}{Character. Sample identifier}
#'   \item{description}{Character. Sample description}
#'   \item{mw}{Numeric. Known weight-average MW in Da}
#'   \item{dispersity}{Numeric. Known dispersity}
#'   \item{dn_dc}{Numeric. Refractive index increment in mL/g}
#'   \item{ext_coef}{Numeric. UV extinction coefficient}
#'   \item{time_min}{Numeric. Elution time in minutes}
#'   \item{ri_mv}{Numeric. RI detector signal in millivolts}
#'   \item{uv_au}{Numeric. UV detector signal in absorbance units}
#'   \item{mals_mv}{Numeric. MALS detector signal in millivolts}
#'   \item{delay_uv_ml}{Numeric. True UV detector delay in mL (negative = before RI)}
#'   \item{delay_mals_ml}{Numeric. True MALS detector delay in mL (positive = after RI)}
#' }
#'
#' @details
#' This dataset simulates raw multi-detector SEC data before inter-detector
#' delay correction. The detectors are physically separated in the flow path,
#' so peaks appear at different times in each detector.
#'
#' **Detector Configuration:**
#' \itemize{
#'   \item UV detector is 0.08 mL BEFORE RI (peaks appear earlier)
#'   \item MALS detector is 0.18 mL AFTER RI (peaks appear later)
#'   \item RI is the reference detector (delay = 0)
#' }
#'
#' **Samples Included:**
#' \itemize{
#'   \item \strong{PS-DelayStd}: Narrow PS standard for determining delays
#'   \item \strong{Sample-1}: PS sample with strong UV absorption
#'   \item \strong{Sample-2}: PMMA sample with weak UV
#'   \item \strong{Sample-3}: Copolymer sample
#' }
#'
#' **Tutorial Workflow:**
#' 1. Load raw multi-detector data
#' 2. Use delay standard to determine inter-detector offsets
#' 3. Apply \code{\link{step_sec_detector_delay}} to align signals
#' 4. Process with \code{\link{step_sec_mals}} for absolute MW
#'
#' @seealso
#' \code{\link{step_sec_detector_delay}} for delay correction
#' \code{\link{step_sec_mals}} for MALS processing
#' \code{\link{sec_triple_detect}} for pre-processed multi-detector data
#'
#' @family sec-data
#' @family sec-raw
#'
#' @examples
#' data(sec_raw_multidetector)
#'
#' # View sample information
#' unique(sec_raw_multidetector[, c("sample_id", "description", "mw")])
#'
#' # Plot delay standard showing detector offset (peaks not aligned)
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   library(ggplot2)
#'   library(dplyr)
#'   library(tidyr)
#'
#'   sec_raw_multidetector |>
#'     filter(sample_id == "PS-DelayStd") |>
#'     select(time_min, ri_mv, uv_au, mals_mv) |>
#'     # Normalize for comparison
#'     mutate(
#'       ri_norm = ri_mv / max(ri_mv),
#'       uv_norm = uv_au / max(uv_au),
#'       mals_norm = mals_mv / max(mals_mv)
#'     ) |>
#'     select(time_min, RI = ri_norm, UV = uv_norm, MALS = mals_norm) |>
#'     pivot_longer(-time_min, names_to = "detector", values_to = "signal") |>
#'     ggplot(aes(time_min, signal, color = detector)) +
#'     geom_line() +
#'     labs(
#'       x = "Time (min)",
#'       y = "Normalized Signal",
#'       title = "Raw Multi-Detector Data (Before Delay Correction)",
#'       subtitle = "Note: Peaks are NOT aligned - UV leads, MALS lags"
#'     ) +
#'     theme_minimal()
#' }
#'
#' @source Synthetic data generated for package tutorials.
"sec_raw_multidetector"
