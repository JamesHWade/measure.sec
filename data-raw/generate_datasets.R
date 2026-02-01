# ==============================================================================
# Generate SEC Triple-Detector Dataset
#
# Creates realistic multi-detector SEC data for testing and examples.
# ==============================================================================

library(tibble)
library(dplyr)
library(purrr)

set.seed(42)

# -----------------------------------------------------------------------------
# Parameters
# -----------------------------------------------------------------------------

# Elution time range (minutes)
time_range <- c(5, 25)
time_points <- seq(time_range[1], time_range[2], by = 0.01)
n_points <- length(time_points)

# Flow rate (mL/min)
flow_rate <- 1.0

# Inter-detector delays (mL) - positive means detector sees sample later than RI
# UV is before RI (negative), MALS is after RI (positive)
delay_volumes <- c(ri = 0, uv = -0.05, mals = 0.15)

# Convert delay volumes to time offsets
delay_times <- delay_volumes / flow_rate # in minutes

# Noise levels (as fraction of peak height)
noise_levels <- c(ri = 0.002, uv = 0.003, mals = 0.005)

# Baseline drift amplitude
baseline_drift <- 0.001

# -----------------------------------------------------------------------------
# Sample Definitions
# -----------------------------------------------------------------------------

# Polymer properties
# dn/dc values in mL/g at 633 nm
# extinction_coef in mL/(mg*cm) at 280 nm (for UV-active polymers)

samples <- tribble(
	~sample_id,
	~sample_type,
	~polymer_type,
	~known_mw,
	~known_dispersity,
	~dn_dc,
	~extinction_coef,
	# Polystyrene standards (narrow dispersity, UV-active)
	"PS-1K",
	"standard",
	"polystyrene",
	1000,
	1.05,
	0.185,
	1.2,
	"PS-10K",
	"standard",
	"polystyrene",
	10000,
	1.03,
	0.185,
	1.2,
	"PS-50K",
	"standard",
	"polystyrene",
	50000,
	1.02,
	0.185,
	1.2,
	"PS-100K",
	"standard",
	"polystyrene",
	100000,
	1.01,
	0.185,
	1.2,
	"PS-500K",
	"standard",
	"polystyrene",
	500000,
	1.01,
	0.185,
	1.2,
	# PMMA samples (broader dispersity, weak UV)
	"PMMA-Low",
	"sample",
	"pmma",
	25000,
	1.8,
	0.084,
	0.1,
	"PMMA-Med",
	"sample",
	"pmma",
	75000,
	2.0,
	0.084,
	0.1,
	"PMMA-High",
	"sample",
	"pmma",
	200000,
	2.2,
	0.084,
	0.1,
	# PEG samples (no UV absorption)
	"PEG-5K",
	"sample",
	"peg",
	5000,
	1.1,
	0.135,
	0.0,
	"PEG-20K",
	"sample",
	"peg",
	20000,
	1.15,
	0.135,
	0.0,
	# Copolymers (variable UV)
	"Copoly-A",
	"sample",
	"copolymer",
	40000,
	1.5,
	0.150,
	0.6,
	"Copoly-B",
	"sample",
	"copolymer",
	80000,
	1.7,
	0.160,
	0.8
)

# -----------------------------------------------------------------------------
# Helper Functions
# -----------------------------------------------------------------------------

# Convert MW to elution time (log-linear calibration)
# Higher MW elutes earlier in SEC
mw_to_time <- function(mw, slope = -0.8, intercept = 20) {
	intercept + slope * log10(mw)
}

# Generate a chromatogram peak (log-normal shape for SEC)
# Peak width increases with dispersity
generate_peak <- function(time, peak_time, dispersity, intensity = 1) {
	# Log-normal parameters
	# sigma controls peak width, related to dispersity
	sigma <- 0.3 + 0.15 * (dispersity - 1)

	# Calculate log-normal shape
	if (peak_time <= 0) {
		return(rep(0, length(time)))
	}

	# Use shifted log-normal for realistic SEC peak shape
	mu <- log(peak_time)
	x_shifted <- time - (peak_time - exp(mu - sigma^2))
	x_shifted[x_shifted <= 0] <- 1e-10

	y <- dlnorm(x_shifted, meanlog = mu, sdlog = sigma)
	y <- y / max(y) * intensity

	# Add slight asymmetry (tailing)
	tail_factor <- 1 + 0.1 * (time - peak_time) / (max(time) - min(time))
	tail_factor[tail_factor < 0.9] <- 0.9
	y <- y * tail_factor

	y
}

# Add baseline drift
add_drift <- function(signal, time, amplitude) {
	# Slow sinusoidal drift
	drift <- amplitude * sin(2 * pi * time / (max(time) - min(time)))
	signal + drift
}

# Add Gaussian noise
add_noise <- function(signal, noise_level) {
	signal + rnorm(length(signal), 0, noise_level * max(abs(signal)))
}

# Apply detector delay (shift signal in time)
apply_delay <- function(signal, time, delay_time) {
	if (abs(delay_time) < 1e-6) {
		return(signal)
	}

	# Interpolate to get delayed signal
	approx(time - delay_time, signal, xout = time, rule = 2)$y
}

# -----------------------------------------------------------------------------
# Generate Chromatograms
# -----------------------------------------------------------------------------

generate_chromatogram <- function(sample_row, time_points) {
	# Extract sample properties
	sample_id <- sample_row$sample_id
	mw <- sample_row$known_mw
	dispersity <- sample_row$known_dispersity
	dn_dc <- sample_row$dn_dc
	ext_coef <- sample_row$extinction_coef

	# Calculate peak position
	peak_time <- mw_to_time(mw)

	# Generate base peak (concentration profile)
	conc_profile <- generate_peak(
		time_points,
		peak_time,
		dispersity,
		intensity = 1
	)

	# RI signal: proportional to concentration * dn/dc
	ri_base <- conc_profile * dn_dc

	# UV signal: proportional to concentration * extinction coefficient
	uv_base <- conc_profile * ext_coef

	# MALS signal: proportional to concentration * MW * (dn/dc)^2
	# Higher MW gives stronger MALS signal
	mals_base <- conc_profile * (mw / 100000) * dn_dc^2

	# Apply inter-detector delays
	# RI is reference (no delay)
	ri_signal <- ri_base
	# UV sees sample slightly before RI

	uv_signal <- apply_delay(uv_base, time_points, delay_times["uv"])
	# MALS sees sample after RI
	mals_signal <- apply_delay(mals_base, time_points, delay_times["mals"])

	# Add baseline drift
	ri_signal <- add_drift(ri_signal, time_points, baseline_drift * dn_dc)
	uv_signal <- add_drift(
		uv_signal,
		time_points,
		baseline_drift * max(ext_coef, 0.1)
	)
	mals_signal <- add_drift(mals_signal, time_points, baseline_drift * dn_dc^2)

	# Add noise
	ri_signal <- add_noise(ri_signal, noise_levels["ri"])
	uv_signal <- add_noise(uv_signal, noise_levels["uv"])
	mals_signal <- add_noise(mals_signal, noise_levels["mals"])

	# Ensure non-negative (detector can't go below zero)
	ri_signal <- pmax(ri_signal, 0)
	uv_signal <- pmax(uv_signal, 0)
	mals_signal <- pmax(mals_signal, 0)

	tibble(
		sample_id = sample_id,
		elution_time = time_points,
		ri_signal = ri_signal,
		uv_signal = uv_signal,
		mals_signal = mals_signal
	)
}

# Generate data for all samples
chromatograms <- samples |>
	split(seq_len(nrow(samples))) |>
	map(generate_chromatogram, time_points = time_points) |>
	bind_rows()

# Join with sample metadata
sec_triple_detect <- chromatograms |>
	left_join(
		select(
			samples,
			sample_id,
			sample_type,
			polymer_type,
			known_mw,
			known_dispersity,
			dn_dc,
			extinction_coef
		),
		by = "sample_id"
	) |>
	select(
		sample_id,
		sample_type,
		polymer_type,
		elution_time,
		ri_signal,
		uv_signal,
		mals_signal,
		known_mw,
		known_dispersity,
		dn_dc,
		extinction_coef
	)

# -----------------------------------------------------------------------------
# Save Dataset
# -----------------------------------------------------------------------------

usethis::use_data(sec_triple_detect, overwrite = TRUE)

# Print summary
cat("\n=== sec_triple_detect Dataset ===\n")
cat("Samples:", nrow(samples), "\n")
cat("Time points per sample:", n_points, "\n")
cat("Total rows:", nrow(sec_triple_detect), "\n")
cat("\nSample types:\n")
print(table(samples$sample_type))
cat("\nPolymer types:\n")
print(table(samples$polymer_type))
cat("\nColumn names:\n")
print(names(sec_triple_detect))
