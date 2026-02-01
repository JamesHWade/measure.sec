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

		# Register SEC-specific peak detection algorithms
		.register_sec_peak_algorithms(pkgname)
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
		list(
			"step_sec_dad",
			"detector",
			"Diode array detector processing (multi-wavelength UV)"
		),
		list(
			"step_sec_mals",
			"detector",
			"Multi-angle light scattering for absolute MW"
		),
		list("step_sec_lals", "detector", "Low-angle light scattering (LALS)"),
		list("step_sec_rals", "detector", "Right-angle light scattering (RALS)"),
		list("step_sec_dls", "detector", "Dynamic light scattering (DLS)"),
		list(
			"step_sec_viscometer",
			"detector",
			"Differential viscometer processing"
		),
		list(
			"step_sec_concentration",
			"calculation",
			"Convert signal to concentration"
		),
		list(
			"step_sec_intrinsic_visc",
			"calculation",
			"Intrinsic viscosity from viscometer"
		),

		# Baseline correction
		list("step_sec_baseline", "baseline", "SEC-optimized baseline correction"),

		# Band broadening correction
		list(
			"step_sec_band_broadening",
			"signal-processing",
			"Axial dispersion (band broadening) correction"
		),

		# Composition analysis
		list(
			"step_sec_uv_ri_ratio",
			"composition",
			"UV/RI ratio for heterogeneity"
		),
		list(
			"step_sec_composition",
			"composition",
			"Copolymer composition from UV/RI"
		),

		# Molecular weight calculations
		list(
			"step_sec_mw_averages",
			"calculation",
			"Calculate Mn, Mw, Mz, dispersity"
		),
		list("step_sec_mw_fractions", "calculation", "Calculate MW fractions"),
		list(
			"step_sec_mw_distribution",
			"processing",
			"Generate MW distribution curves"
		),

		# Protein SEC
		list(
			"step_sec_aggregates",
			"protein",
			"Aggregate and fragment quantitation"
		),
		list(
			"step_sec_oligomer",
			"protein",
			"Detailed oligomeric species analysis"
		),
		list(
			"step_sec_protein",
			"protein",
			"Protein SEC workflow (baseline + aggregates + oligomer)"
		),

		# Calibration
		list(
			"step_sec_conventional_cal",
			"calibration",
			"Conventional calibration using narrow standards"
		),
		list(
			"step_sec_universal_cal",
			"calibration",
			"Universal calibration using Mark-Houwink"
		),

		# Peak detection
		list(
			"step_sec_peaks_detect",
			"peaks",
			"SEC-optimized peak detection with finderskeepers algorithm"
		),

		# Integration
		list(
			"step_sec_integration_window",
			"integration",
			"Define integration window for MW calculations"
		),

		# Exclusion regions
		list(
			"step_sec_exclude_regions",
			"preprocessing",
			"Mark regions for exclusion from baseline/integration"
		),

		# Peak boundary refinement
		list(
			"step_sec_peaks_refine",
			"peaks",
			"Refine peak boundaries using height fraction method"
		),

		# Peak deconvolution
		list(
			"step_sec_peaks_deconvolve",
			"peaks",
			"SEC-optimized peak deconvolution with EMG model"
		)
	)

	for (s in steps) {
		measure::register_measure_step(s[[1]], pkgname, s[[2]], s[[3]])
	}
}

#' Register SEC-specific Peak Detection Algorithms
#'
#' Registers SEC-specific peak detection algorithms with the measure registry.
#' Currently registers the finderskeepers algorithm.
#'
#' @param pkgname Package name (passed from .onLoad).
#' @noRd
.register_sec_peak_algorithms <- function(pkgname) {
	# Note: The finderskeepers algorithm is implemented directly in
	# step_sec_peaks_detect.R rather than registered with measure's registry.
	# This is because the algorithm has SEC-specific parameters (loess_span,
	# ist_points, ist_nonlinearity) that don't fit the standard measure
	# algorithm interface. The step handles algorithm selection internally.
	#
	# If measure's registry adds support for technique-specific parameter
	# sets in the future, this function could register finderskeepers there.

	invisible(NULL)
}
