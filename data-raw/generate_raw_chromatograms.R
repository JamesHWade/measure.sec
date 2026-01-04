# ==============================================================================
# Generate Raw SEC Chromatogram Datasets
#
# Creates realistic "raw" chromatogram data that mimics instrument exports.
# Includes realistic noise, baseline drift, and typical lab conditions.
# ==============================================================================

library(tibble)
library(dplyr)
library(purrr)

set.seed(2026)

# -----------------------------------------------------------------------------
# Common Parameters
# -----------------------------------------------------------------------------

# Flow rate (mL/min)
flow_rate <- 1.0

# Time step (seconds converted to minutes)
time_step <- 0.1 / 60 # 0.1 second sampling = 10 Hz

# -----------------------------------------------------------------------------
# Helper Functions
# -----------------------------------------------------------------------------

# Convert MW to elution time (log-linear calibration)
# Higher MW elutes earlier in SEC
mw_to_time <- function(mw, slope = -0.75, intercept = 19) {
  intercept + slope * log10(mw)
}

# Generate a chromatogram peak (log-normal shape for SEC)
generate_peak <- function(time, peak_time, dispersity, intensity = 1) {
  sigma <- 0.25 + 0.12 * (dispersity - 1)
  if (peak_time <= 0) {
    return(rep(0, length(time)))
  }

  mu <- log(peak_time)
  x_shifted <- time - (peak_time - exp(mu - sigma^2))
  x_shifted[x_shifted <= 0] <- 1e-10

  y <- dlnorm(x_shifted, meanlog = mu, sdlog = sigma)
  y <- y / max(y) * intensity

  # Add asymmetry (tailing) - common in SEC
  tail_factor <- 1 + 0.15 * (time - peak_time) / (max(time) - min(time))
  tail_factor[tail_factor < 0.85] <- 0.85
  y * tail_factor
}

# Add realistic baseline drift (slow drift + injection artifact)
add_baseline <- function(signal, time, drift_amp, injection_time = NULL) {
  n <- length(time)
  t_range <- max(time) - min(time)

  # Slow sinusoidal drift (column temperature fluctuation)
  slow_drift <- drift_amp * sin(2 * pi * (time - min(time)) / t_range * 0.7)

  # Very slow linear drift (pump flow variation)
  linear_drift <- drift_amp * 0.5 * (time - min(time)) / t_range

  # Injection artifact (spike at injection time)
  injection_artifact <- rep(0, n)
  if (!is.null(injection_time)) {
    inj_idx <- which.min(abs(time - injection_time))
    # Sharp spike that decays
    decay_time <- 0.3 # minutes
    decay_idx <- which(
      time >= injection_time & time <= injection_time + decay_time
    )
    if (length(decay_idx) > 0) {
      t_decay <- time[decay_idx] - injection_time
      injection_artifact[decay_idx] <- drift_amp * 3 * exp(-t_decay / 0.05)
    }
  }

  signal + slow_drift + linear_drift + injection_artifact
}

# Add realistic noise (heteroscedastic - noise scales with signal)
add_noise <- function(signal, base_noise, hetero_factor = 0.02) {
  base <- rnorm(length(signal), 0, base_noise)
  hetero <- rnorm(length(signal), 0, hetero_factor * abs(signal))
  signal + base + hetero
}

# Apply detector delay (shift signal in time)
apply_delay <- function(signal, time, delay_time) {
  if (abs(delay_time) < 1e-6) {
    return(signal)
  }
  approx(time - delay_time, signal, xout = time, rule = 2)$y
}

# =============================================================================
# DATASET 1: Raw PS Standards for Conventional Calibration
# =============================================================================

generate_raw_ps_standards <- function() {
  cat("Generating sec_raw_standards...\n")

  # Time range (minutes) - typical SEC run
  time_points <- seq(4, 22, by = time_step)
  n_points <- length(time_points)

  # PS narrow standards (typical commercial kit)
  standards <- tribble(
    ~standard_name , ~mp    , ~dispersity ,
    "PS-580"       ,    580 , 1.06        ,
    "PS-1270"      ,   1270 , 1.04        ,
    "PS-2960"      ,   2960 , 1.03        ,
    "PS-5970"      ,   5970 , 1.02        ,
    "PS-9680"      ,   9680 , 1.02        ,
    "PS-19600"     ,  19600 , 1.02        ,
    "PS-33500"     ,  33500 , 1.02        ,
    "PS-67500"     ,  67500 , 1.02        ,
    "PS-135000"    , 135000 , 1.01        ,
    "PS-270000"    , 270000 , 1.01        ,
    "PS-495000"    , 495000 , 1.02        ,
    "PS-930000"    , 930000 , 1.02
  )

  # Generate chromatogram for each standard
  generate_standard_chrom <- function(std_row) {
    peak_time <- mw_to_time(std_row$mp)

    # Base peak (RI detector response)
    ri_signal <- generate_peak(
      time_points,
      peak_time,
      std_row$dispersity,
      intensity = 150 + runif(1, -20, 20) # mV-like values with injection variation
    )

    # Add realistic baseline and noise
    ri_signal <- add_baseline(
      ri_signal,
      time_points,
      drift_amp = 1.5 + runif(1, -0.5, 0.5),
      injection_time = min(time_points)
    )
    ri_signal <- add_noise(ri_signal, base_noise = 0.8)

    tibble(
      standard_name = std_row$standard_name,
      time_min = round(time_points, 4),
      ri_mv = round(ri_signal, 3)
    )
  }

  chromatograms <- standards |>
    split(seq_len(nrow(standards))) |>
    map(generate_standard_chrom) |>
    bind_rows()

  # Add metadata
  sec_raw_standards <- chromatograms |>
    left_join(
      standards |> mutate(log_mp = log10(mp)),
      by = "standard_name"
    ) |>
    select(standard_name, mp, log_mp, dispersity, time_min, ri_mv)

  sec_raw_standards
}

# =============================================================================
# DATASET 2: Raw Unknown Samples with Known MW for Validation
# =============================================================================

generate_raw_unknowns <- function() {
  cat("Generating sec_raw_unknowns...\n")

  time_points <- seq(4, 22, by = time_step)

  # Unknown samples with known true MW for validation
  unknowns <- tribble(
    ~sample_id        , ~true_mw , ~true_mn , ~true_mz , ~true_dispersity , ~description                    ,
    "Unknown-A"       ,    45000 ,    22000 ,    85000 , 2.05             , "Broad distribution PMMA-like"  ,
    "Unknown-B"       ,   125000 ,    95000 ,   165000 , 1.32             , "Medium dispersity PS-like"     ,
    "Unknown-C"       ,    82000 ,    75000 ,    92000 , 1.09             , "Narrow distribution reference" ,
    "Unknown-Bimodal" , NA       , NA       , NA       , NA               , "Bimodal mixture (50K + 200K)"  ,
    "Unknown-HMW"     ,  1500000 ,  1200000 ,  1900000 , 1.25             , "Very high MW with aggregates"  ,
    "Unknown-LMW"     ,     3500 ,     2800 ,     4500 , 1.25             , "Low MW oligomer region"
  )

  generate_unknown_chrom <- function(unk_row) {
    if (unk_row$sample_id == "Unknown-Bimodal") {
      # Bimodal: two peaks at 50K and 200K
      peak1 <- generate_peak(
        time_points,
        mw_to_time(50000),
        1.15,
        intensity = 80
      )
      peak2 <- generate_peak(
        time_points,
        mw_to_time(200000),
        1.2,
        intensity = 120
      )
      ri_signal <- peak1 + peak2
    } else if (unk_row$sample_id == "Unknown-HMW") {
      # High MW with aggregate shoulder
      main_peak <- generate_peak(
        time_points,
        mw_to_time(unk_row$true_mw),
        unk_row$true_dispersity,
        intensity = 100
      )
      # Small aggregate peak at higher MW (earlier elution)
      agg_peak <- generate_peak(
        time_points,
        mw_to_time(5000000),
        1.3,
        intensity = 8
      )
      ri_signal <- main_peak + agg_peak
    } else {
      ri_signal <- generate_peak(
        time_points,
        mw_to_time(unk_row$true_mw),
        unk_row$true_dispersity,
        intensity = 120 + runif(1, -30, 30)
      )
    }

    # More variable baselines for unknowns (real lab conditions)
    ri_signal <- add_baseline(
      ri_signal,
      time_points,
      drift_amp = 2.0 + runif(1, -1, 1),
      injection_time = min(time_points)
    )
    ri_signal <- add_noise(ri_signal, base_noise = 1.0, hetero_factor = 0.025)

    tibble(
      sample_id = unk_row$sample_id,
      time_min = round(time_points, 4),
      ri_mv = round(ri_signal, 3)
    )
  }

  chromatograms <- unknowns |>
    split(seq_len(nrow(unknowns))) |>
    map(generate_unknown_chrom) |>
    bind_rows()

  sec_raw_unknowns <- chromatograms |>
    left_join(unknowns, by = "sample_id") |>
    select(
      sample_id,
      description,
      true_mw,
      true_mn,
      true_mz,
      true_dispersity,
      time_min,
      ri_mv
    )

  sec_raw_unknowns
}

# =============================================================================
# DATASET 3: Raw Multi-Detector Data (RI, UV, MALS)
# =============================================================================

generate_raw_multidetector <- function() {
  cat("Generating sec_raw_multidetector...\n")

  time_points <- seq(4, 24, by = time_step)

  # Inter-detector delays (NOT corrected - this is raw data)
  # Typical order: UV -> RI -> MALS (UV first, MALS last)
  delay_uv <- -0.08 # UV is 0.08 mL before RI

  delay_mals <- 0.18 # MALS is 0.18 mL after RI

  # Samples for multi-detector analysis
  samples <- tribble(
    ~sample_id    , ~mw    , ~dispersity , ~dn_dc , ~ext_coef , ~description                        ,
    "PS-DelayStd" , 100000 , 1.02        , 0.185  , 1.2       , "Narrow PS for delay determination" ,
    "Sample-1"    ,  75000 , 1.85        , 0.185  , 1.1       , "PS sample with UV absorption"      ,
    "Sample-2"    , 150000 , 2.1         , 0.085  , 0.05      , "PMMA sample (weak UV)"             ,
    "Sample-3"    ,  45000 , 1.5         , 0.150  , 0.8       , "Copolymer sample"
  )

  generate_multidet_chrom <- function(smp_row) {
    peak_time <- mw_to_time(smp_row$mw)

    # RI signal (reference detector)
    conc_profile <- generate_peak(
      time_points,
      peak_time,
      smp_row$dispersity,
      intensity = 200 * smp_row$dn_dc
    )

    ri_signal <- conc_profile
    ri_signal <- add_baseline(
      ri_signal,
      time_points,
      drift_amp = 1.2,
      injection_time = min(time_points)
    )
    ri_signal <- add_noise(ri_signal, base_noise = 0.6)

    # UV signal (BEFORE RI in flow path - sees sample earlier)
    uv_base <- conc_profile * smp_row$ext_coef / smp_row$dn_dc
    uv_signal <- apply_delay(uv_base, time_points, delay_uv)
    uv_signal <- add_baseline(
      uv_signal,
      time_points,
      drift_amp = 0.8,
      injection_time = min(time_points) + delay_uv
    )
    uv_signal <- add_noise(uv_signal, base_noise = 0.4)

    # MALS signal (AFTER RI in flow path - sees sample later)
    # MALS scales with MW for absolute MW determination
    mals_base <- conc_profile * (smp_row$mw / 50000) * smp_row$dn_dc
    mals_signal <- apply_delay(mals_base, time_points, delay_mals)
    mals_signal <- add_baseline(
      mals_signal,
      time_points,
      drift_amp = 2.5,
      injection_time = min(time_points) + delay_mals
    )
    mals_signal <- add_noise(
      mals_signal,
      base_noise = 1.5,
      hetero_factor = 0.03
    )

    tibble(
      sample_id = smp_row$sample_id,
      time_min = round(time_points, 4),
      ri_mv = round(pmax(ri_signal, 0), 3),
      uv_au = round(pmax(uv_signal, 0) / 1000, 6), # AU units
      mals_mv = round(pmax(mals_signal, 0), 3)
    )
  }

  chromatograms <- samples |>
    split(seq_len(nrow(samples))) |>
    map(generate_multidet_chrom) |>
    bind_rows()

  sec_raw_multidetector <- chromatograms |>
    left_join(samples, by = "sample_id") |>
    mutate(
      delay_uv_ml = delay_uv,
      delay_mals_ml = delay_mals
    ) |>
    select(
      sample_id,
      description,
      mw,
      dispersity,
      dn_dc,
      ext_coef,
      time_min,
      ri_mv,
      uv_au,
      mals_mv,
      delay_uv_ml,
      delay_mals_ml
    )

  sec_raw_multidetector
}

# =============================================================================
# Generate and Save All Datasets
# =============================================================================

sec_raw_standards <- generate_raw_ps_standards()
sec_raw_unknowns <- generate_raw_unknowns()
sec_raw_multidetector <- generate_raw_multidetector()

usethis::use_data(sec_raw_standards, overwrite = TRUE)
usethis::use_data(sec_raw_unknowns, overwrite = TRUE)
usethis::use_data(sec_raw_multidetector, overwrite = TRUE)

# =============================================================================
# Summary
# =============================================================================

cat("\n=== Dataset Summary ===\n\n")

cat("sec_raw_standards:\n")
cat("  - 12 narrow PS standards (580 Da to 930K Da)\n")
cat("  - Columns:", paste(names(sec_raw_standards), collapse = ", "), "\n")
cat("  - Rows:", nrow(sec_raw_standards), "\n")
cat("  - For conventional calibration tutorials\n\n")

cat("sec_raw_unknowns:\n")
cat("  - 6 unknown samples with known true MW\n")
cat("  - Includes bimodal, HMW with aggregates, low MW\n")
cat("  - Columns:", paste(names(sec_raw_unknowns), collapse = ", "), "\n")
cat("  - Rows:", nrow(sec_raw_unknowns), "\n")
cat("  - For MW calculation validation\n\n")

cat("sec_raw_multidetector:\n")
cat("  - 4 samples with RI, UV, MALS signals\n")
cat("  - Includes delay standard for inter-detector offset\n")
cat("  - Delays NOT corrected (raw data)\n")
cat("  - Columns:", paste(names(sec_raw_multidetector), collapse = ", "), "\n")
cat("  - Rows:", nrow(sec_raw_multidetector), "\n")
cat("  - For triple detection tutorials\n")
