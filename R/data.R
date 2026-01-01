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
