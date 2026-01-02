# ==============================================================================
# Generate Additional SEC Example Datasets
#
# Creates specialized datasets for:
# - Copolymer composition analysis (sec_copolymer)
# - Protein SEC with aggregates (sec_protein)
# - Branched polymer analysis (sec_branched)
# - System suitability testing (sec_system_suitability)
# ==============================================================================

library(tibble)
library(dplyr)
library(purrr)

set.seed(42)

# =============================================================================
# Common Helper Functions
# =============================================================================

# Convert MW to elution time (log-linear calibration)
mw_to_time <- function(mw, slope = -0.8, intercept = 20) {
  intercept + slope * log10(mw)
}

# Generate a chromatogram peak (log-normal shape)
generate_peak <- function(time, peak_time, width = 0.3, intensity = 1) {
  if (peak_time <= 0 || is.na(peak_time)) {
    return(rep(0, length(time)))
  }

  mu <- log(peak_time)
  x_shifted <- time - (peak_time - exp(mu - width^2))
  x_shifted[x_shifted <= 0] <- 1e-10

  y <- dlnorm(x_shifted, meanlog = mu, sdlog = width)
  y <- y / max(y, na.rm = TRUE) * intensity
  y[is.na(y)] <- 0
  y
}

# Add Gaussian noise
add_noise <- function(signal, noise_level) {
  signal +
    rnorm(length(signal), 0, noise_level * max(abs(signal), na.rm = TRUE))
}

# Add baseline
add_baseline <- function(signal, offset = 0, drift = 0, time = NULL) {
  if (is.null(time)) {
    return(signal + offset)
  }
  baseline <- offset + drift * (time - min(time)) / (max(time) - min(time))
  signal + baseline
}

# =============================================================================
# 1. COPOLYMER DATASET (sec_copolymer)
# =============================================================================
# Demonstrates UV/RI ratio analysis for composition determination
# Styrene-acrylate copolymers with varying styrene content

cat("\n=== Generating sec_copolymer Dataset ===\n")

time_points <- seq(8, 22, by = 0.02)

# Copolymer samples with different compositions
copolymer_samples <- tribble(
  ~sample_id   , ~styrene_fraction , ~mw   , ~dispersity , ~description                ,
  "Copoly-20S" , 0.20              , 45000 , 1.8         , "20% styrene, 80% acrylate" ,
  "Copoly-40S" , 0.40              , 52000 , 1.9         , "40% styrene, 60% acrylate" ,
  "Copoly-60S" , 0.60              , 48000 , 1.7         , "60% styrene, 40% acrylate" ,
  "Copoly-80S" , 0.80              , 55000 , 2.0         , "80% styrene, 20% acrylate" ,
  "PS-Homo"    , 1.00              , 50000 , 1.5         , "Polystyrene homopolymer"   ,
  "PA-Homo"    , 0.00              , 47000 , 1.6         , "Polyacrylate homopolymer"
)

# dn/dc values
dndc_styrene <- 0.185 # Polystyrene
dndc_acrylate <- 0.068 # Polyacrylate

# UV extinction at 254 nm (styrene absorbs, acrylate doesn't)
ext_styrene <- 1.5
ext_acrylate <- 0.02

generate_copolymer_chrom <- function(sample, time_points) {
  peak_time <- mw_to_time(sample$mw)
  width <- 0.25 + 0.1 * (sample$dispersity - 1)

  # Concentration profile
  conc <- generate_peak(time_points, peak_time, width, intensity = 1)

  # Calculate composite properties
  f_sty <- sample$styrene_fraction
  dndc_mix <- f_sty * dndc_styrene + (1 - f_sty) * dndc_acrylate
  ext_mix <- f_sty * ext_styrene + (1 - f_sty) * ext_acrylate

  # RI signal (proportional to dn/dc)
  ri_signal <- conc * dndc_mix
  ri_signal <- add_baseline(
    ri_signal,
    offset = 0.001,
    drift = 0.0005,
    time = time_points
  )
  ri_signal <- add_noise(ri_signal, 0.003)
  ri_signal <- pmax(ri_signal, 0)

  # UV signal at 254 nm (styrene-selective)
  uv_signal <- conc * ext_mix
  uv_signal <- add_baseline(
    uv_signal,
    offset = 0.002,
    drift = 0.001,
    time = time_points
  )
  uv_signal <- add_noise(uv_signal, 0.004)
  uv_signal <- pmax(uv_signal, 0)

  tibble(
    sample_id = sample$sample_id,
    elution_time = time_points,
    ri_signal = ri_signal,
    uv_254_signal = uv_signal,
    styrene_fraction = sample$styrene_fraction,
    mw = sample$mw,
    dispersity = sample$dispersity,
    description = sample$description
  )
}

sec_copolymer <- copolymer_samples |>
  split(seq_len(nrow(copolymer_samples))) |>
  map(generate_copolymer_chrom, time_points = time_points) |>
  bind_rows()

cat("Rows:", nrow(sec_copolymer), "\n")
cat("Samples:", length(unique(sec_copolymer$sample_id)), "\n")

# =============================================================================
# 2. PROTEIN DATASET (sec_protein)
# =============================================================================
# Therapeutic protein SEC with monomer, dimer, aggregates, and fragments

cat("\n=== Generating sec_protein Dataset ===\n")

time_points_protein <- seq(5, 30, by = 0.02)

# Protein samples - mAb-like
protein_samples <- tribble(
  ~sample_id        , ~description          , ~monomer_pct , ~dimer_pct , ~hmw_pct , ~fragment_pct ,
  "mAb-Reference"   , "Reference standard"  , 98.5         , 1.0        , 0.3      , 0.2           ,
  "mAb-Stressed-1"  , "Heat stress 40C 1wk" , 94.0         , 3.5        , 1.5      , 1.0           ,
  "mAb-Stressed-2"  , "Heat stress 40C 2wk" , 88.0         , 6.0        , 4.0      , 2.0           ,
  "mAb-Aged"        , "12 month stability"  , 96.0         , 2.5        , 1.0      , 0.5           ,
  "mAb-Freeze-Thaw" , "5x freeze-thaw"      , 95.5         , 2.0        , 2.0      , 0.5
)

# Protein MW (mAb ~150 kDa)
mw_monomer <- 150000
mw_dimer <- 300000
mw_hmw <- 600000 # Higher aggregates
mw_fragment <- 50000 # Fab fragment

# Extinction coefficient at 280 nm (typical for mAb)
ext_280 <- 1.4

generate_protein_chrom <- function(sample, time_points) {
  # Peak positions
  peak_hmw <- mw_to_time(mw_hmw)
  peak_dimer <- mw_to_time(mw_dimer)
  peak_monomer <- mw_to_time(mw_monomer)
  peak_fragment <- mw_to_time(mw_fragment)

  # Generate each species peak (proteins are narrow)
  hmw_peak <- generate_peak(
    time_points,
    peak_hmw,
    width = 0.15,
    intensity = sample$hmw_pct / 100
  )
  dimer_peak <- generate_peak(
    time_points,
    peak_dimer,
    width = 0.12,
    intensity = sample$dimer_pct / 100
  )
  monomer_peak <- generate_peak(
    time_points,
    peak_monomer,
    width = 0.10,
    intensity = sample$monomer_pct / 100
  )
  fragment_peak <- generate_peak(
    time_points,
    peak_fragment,
    width = 0.18,
    intensity = sample$fragment_pct / 100
  )

  # Combined profile
  total_signal <- hmw_peak + dimer_peak + monomer_peak + fragment_peak

  # UV 280 nm signal
  uv_280 <- total_signal * ext_280
  uv_280 <- add_baseline(
    uv_280,
    offset = 0.005,
    drift = 0.002,
    time = time_points
  )
  uv_280 <- add_noise(uv_280, 0.002)
  uv_280 <- pmax(uv_280, 0)

  # UV 214 nm signal (peptide bond, more sensitive)
  uv_214 <- total_signal * ext_280 * 15 # Much higher response
  uv_214 <- add_baseline(
    uv_214,
    offset = 0.02,
    drift = 0.005,
    time = time_points
  )
  uv_214 <- add_noise(uv_214, 0.003)
  uv_214 <- pmax(uv_214, 0)

  tibble(
    sample_id = sample$sample_id,
    elution_time = time_points,
    uv_280_signal = uv_280,
    uv_214_signal = uv_214,
    description = sample$description,
    monomer_pct = sample$monomer_pct,
    dimer_pct = sample$dimer_pct,
    hmw_pct = sample$hmw_pct,
    fragment_pct = sample$fragment_pct
  )
}

sec_protein <- protein_samples |>
  split(seq_len(nrow(protein_samples))) |>
  map(generate_protein_chrom, time_points = time_points_protein) |>
  bind_rows()

cat("Rows:", nrow(sec_protein), "\n")
cat("Samples:", length(unique(sec_protein$sample_id)), "\n")

# =============================================================================
# 3. BRANCHED POLYMER DATASET (sec_branched)
# =============================================================================
# Linear vs branched polymers showing Mark-Houwink differences

cat("\n=== Generating sec_branched Dataset ===\n")

time_points_branched <- seq(8, 22, by = 0.02)

# Branched vs linear samples at similar MW
branched_samples <- tribble(
  ~sample_id    , ~topology  , ~mw    , ~branching_index , ~intrinsic_visc , ~rg  ,
  "Linear-50K"  , "linear"   ,  50000 , 1.00             ,  45.0           ,  8.5 ,
  "Linear-100K" , "linear"   , 100000 , 1.00             ,  72.0           , 12.0 ,
  "Linear-200K" , "linear"   , 200000 , 1.00             , 115.0           , 17.0 ,
  "Branch-50K"  , "branched" ,  50000 , 0.75             ,  34.0           ,  7.0 ,
  "Branch-100K" , "branched" , 100000 , 0.65             ,  47.0           ,  9.5 ,
  "Branch-200K" , "branched" , 200000 , 0.55             ,  63.0           , 12.5 ,
  "Star-50K"    , "star"     ,  50000 , 0.60             ,  27.0           ,  5.8 ,
  "Star-100K"   , "star"     , 100000 , 0.50             ,  36.0           ,  7.5
)

# dn/dc for polyethylene (common for branching studies)
dndc_pe <- 0.104

generate_branched_chrom <- function(sample, time_points) {
  # Branched polymers elute LATER than linear at same MW
  # (smaller hydrodynamic volume)
  base_time <- mw_to_time(sample$mw)

  # Shift based on branching (g' = branching index)
  # More branching = later elution
  time_shift <- -0.5 * log(sample$branching_index)
  peak_time <- base_time + time_shift

  width <- 0.30 # Moderate dispersity

  # Concentration profile
  conc <- generate_peak(time_points, peak_time, width, intensity = 1)

  # RI signal
  ri_signal <- conc * dndc_pe
  ri_signal <- add_baseline(
    ri_signal,
    offset = 0.0005,
    drift = 0.0002,
    time = time_points
  )
  ri_signal <- add_noise(ri_signal, 0.002)
  ri_signal <- pmax(ri_signal, 0)

  # Viscometer signal (proportional to intrinsic viscosity * concentration)
  visc_signal <- conc * sample$intrinsic_visc / 100
  visc_signal <- add_baseline(
    visc_signal,
    offset = 0.001,
    drift = 0.0005,
    time = time_points
  )
  visc_signal <- add_noise(visc_signal, 0.003)
  visc_signal <- pmax(visc_signal, 0)

  # MALS signal (same MW gives same light scattering)
  mals_signal <- conc * (sample$mw / 100000) * dndc_pe^2
  mals_signal <- add_baseline(
    mals_signal,
    offset = 0.0002,
    drift = 0.0001,
    time = time_points
  )
  mals_signal <- add_noise(mals_signal, 0.004)
  mals_signal <- pmax(mals_signal, 0)

  tibble(
    sample_id = sample$sample_id,
    elution_time = time_points,
    ri_signal = ri_signal,
    visc_signal = visc_signal,
    mals_signal = mals_signal,
    topology = sample$topology,
    mw = sample$mw,
    branching_index = sample$branching_index,
    intrinsic_visc = sample$intrinsic_visc,
    rg = sample$rg
  )
}

sec_branched <- branched_samples |>
  split(seq_len(nrow(branched_samples))) |>
  map(generate_branched_chrom, time_points = time_points_branched) |>
  bind_rows()

cat("Rows:", nrow(sec_branched), "\n")
cat("Samples:", length(unique(sec_branched$sample_id)), "\n")

# =============================================================================
# 4. SYSTEM SUITABILITY DATASET (sec_system_suitability)
# =============================================================================
# Column QC data with resolution test standards

cat("\n=== Generating sec_system_suitability Dataset ===\n")

time_points_sst <- seq(5, 25, by = 0.01)

# System suitability runs
sst_samples <- tribble(
  ~run_id      , ~column_age_days , ~plate_count , ~asymmetry , ~resolution ,
  "SST-Day0"   ,                0 ,        45000 , 1.05       , 2.8         ,
  "SST-Day30"  ,               30 ,        43000 , 1.08       , 2.7         ,
  "SST-Day60"  ,               60 ,        40000 , 1.12       , 2.5         ,
  "SST-Day90"  ,               90 ,        36000 , 1.18       , 2.3         ,
  "SST-Day120" ,              120 ,        32000 , 1.25       , 2.0
)

# Resolution standard: two narrow PS standards
mw_std1 <- 100000
mw_std2 <- 50000

generate_sst_chrom <- function(sample, time_points) {
  peak1_time <- mw_to_time(mw_std1)
  peak2_time <- mw_to_time(mw_std2)

  # Width increases as column ages (lower plate count)
  base_width <- 0.08
  width_factor <- sqrt(45000 / sample$plate_count)
  width <- base_width * width_factor

  # Asymmetry factor affects peak shape
  asym <- sample$asymmetry

  # Generate peaks
  peak1 <- generate_peak(time_points, peak1_time, width, intensity = 1)
  peak2 <- generate_peak(time_points, peak2_time, width, intensity = 1)

  # Apply asymmetry (tail on the right side)
  if (asym > 1.0) {
    tail_idx <- time_points > peak1_time
    peak1[tail_idx] <- peak1[tail_idx] *
      (1 + (asym - 1) * 0.3 * (time_points[tail_idx] - peak1_time))
    tail_idx2 <- time_points > peak2_time
    peak2[tail_idx2] <- peak2[tail_idx2] *
      (1 + (asym - 1) * 0.3 * (time_points[tail_idx2] - peak2_time))
  }

  # Normalize
  peak1 <- peak1 / max(peak1)
  peak2 <- peak2 / max(peak2)

  # Combined signal
  ri_signal <- (peak1 + peak2) * 0.185 # PS dn/dc
  ri_signal <- add_baseline(
    ri_signal,
    offset = 0.0002,
    drift = 0.0001,
    time = time_points
  )
  ri_signal <- add_noise(ri_signal, 0.001)
  ri_signal <- pmax(ri_signal, 0)

  tibble(
    run_id = sample$run_id,
    elution_time = time_points,
    ri_signal = ri_signal,
    column_age_days = sample$column_age_days,
    expected_plate_count = sample$plate_count,
    expected_asymmetry = sample$asymmetry,
    expected_resolution = sample$resolution
  )
}

sec_system_suitability <- sst_samples |>
  split(seq_len(nrow(sst_samples))) |>
  map(generate_sst_chrom, time_points = time_points_sst) |>
  bind_rows()

cat("Rows:", nrow(sec_system_suitability), "\n")
cat("Runs:", length(unique(sec_system_suitability$run_id)), "\n")

# =============================================================================
# Save All Datasets
# =============================================================================

cat("\n=== Saving Datasets ===\n")

usethis::use_data(sec_copolymer, overwrite = TRUE)
usethis::use_data(sec_protein, overwrite = TRUE)
usethis::use_data(sec_branched, overwrite = TRUE)
usethis::use_data(sec_system_suitability, overwrite = TRUE)

cat("\nDone! Created:\n")
cat("- sec_copolymer:", nrow(sec_copolymer), "rows\n")
cat("- sec_protein:", nrow(sec_protein), "rows\n")
cat("- sec_branched:", nrow(sec_branched), "rows\n")
cat("- sec_system_suitability:", nrow(sec_system_suitability), "rows\n")
