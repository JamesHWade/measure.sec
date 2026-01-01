# ==============================================================================
# Package Registration
#
# Registers measure.sec with the measure package on load.
# ==============================================================================

.onLoad <- function(libname, pkgname) {
  if (requireNamespace("measure", quietly = TRUE)) {
    # Register as a technique pack (use pkgname for portability)
    measure::register_measure_pack(
      pack_name = pkgname,
      technique = "SEC/GPC",
      description = "Size Exclusion / Gel Permeation Chromatography"
    )

    # Register all SEC steps
    .register_sec_steps(pkgname)
  }
}

#' Register SEC/GPC Steps
#'
#' Registers all steps from this package with the measure registry.
#'
#' @param pkgname Package name (passed from .onLoad).
#' @noRd
.register_sec_steps <- function(pkgname) {
  # Step registration data
  # Format: list(step_name, category, description)
  steps <- list(
    # Preprocessing
    list(
      "step_sec_detector_delay",
      "preprocessing",
      "Correct inter-detector volume delays"
    ),

    # Detector processing
    list("step_sec_ri", "detector", "RI detector with dn/dc handling"),
    list("step_sec_uv", "detector", "UV detector with extinction coefficient"),
    list("step_sec_mals", "detector", "Multi-angle light scattering for absolute MW"),
    list("step_sec_viscometer", "detector", "Differential viscometer processing"),
    list("step_sec_concentration", "calculation", "Convert signal to concentration"),
    list("step_sec_intrinsic_visc", "calculation", "Intrinsic viscosity from viscometer"),

    # Baseline correction
    list("step_sec_baseline", "baseline", "SEC-optimized baseline correction"),

    # Composition analysis
    list("step_sec_uv_ri_ratio", "composition", "UV/RI ratio for heterogeneity"),
    list("step_sec_composition", "composition", "Copolymer composition from UV/RI"),

    # Molecular weight calculations
    list("step_sec_mw_averages", "calculation", "Calculate Mn, Mw, Mz, dispersity"),
    list("step_sec_mw_fractions", "calculation", "Calculate MW fractions"),
    list("step_sec_mw_distribution", "processing", "Generate MW distribution curves"),

    # Protein SEC
    list("step_sec_aggregates", "protein", "Aggregate and fragment quantitation")

    # TODO: Add when implemented:
    # list("step_sec_calibrate", "calibration", "Calibration curve fitting")
  )

  for (s in steps) {
    measure::register_measure_step(s[[1]], pkgname, s[[2]], s[[3]])
  }
}
