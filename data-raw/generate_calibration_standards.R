# ==============================================================================
# Generate SEC Calibration Standards Dataset
#
# Creates comprehensive calibration standard data based on commercial standard
# kits (EasiVial, ReadyCal) with realistic certificate values and retention data.
#
# Designed for separations scientists to:
# 1. Build conventional (narrow standard) calibration curves
# 2. Perform universal calibration with Mark-Houwink parameters
# 3. Compare calibrations across polymer types
# 4. Assess calibration quality with certificate uncertainties
# ==============================================================================

library(tibble)
library(dplyr)

set.seed(42)

# ==============================================================================
# Polystyrene Standards
# Based on typical commercial kits (e.g., Agilent EasiVial, PSS ReadyCal)
# THF mobile phase, 1.0 mL/min, dual PLgel Mixed-C columns
# ==============================================================================

# Certificate values based on typical narrow PS standards
# Mp values are what appear on most certificates
# Mn and Mw calculated from Mp and typical dispersities

ps_standards <- tribble(
  ~standard_name , ~mp     , ~mn     , ~mw     , ~dispersity , ~mp_uncertainty ,
  # Low MW range
  "PS-162"       ,     162 ,     160 ,     163 , 1.02        , 0.10            ,
  "PS-580"       ,     580 ,     570 ,     590 , 1.04        , 0.08            ,
  "PS-1050"      ,    1050 ,    1020 ,    1070 , 1.05        , 0.06            ,
  "PS-2970"      ,    2970 ,    2890 ,    3020 , 1.05        , 0.05            ,
  # Mid MW range
  "PS-5030"      ,    5030 ,    4900 ,    5120 , 1.04        , 0.04            ,
  "PS-9680"      ,    9680 ,    9470 ,    9820 , 1.04        , 0.04            ,
  "PS-19800"     ,   19800 ,   19400 ,   20100 , 1.04        , 0.03            ,
  "PS-33500"     ,   33500 ,   32800 ,   34000 , 1.04        , 0.03            ,
  # High MW range
  "PS-67500"     ,   67500 ,   66200 ,   68500 , 1.03        , 0.03            ,
  "PS-120000"    ,  120000 ,  118000 ,  122000 , 1.03        , 0.03            ,
  "PS-216000"    ,  216000 ,  212000 ,  219000 , 1.03        , 0.03            ,
  "PS-430000"    ,  430000 ,  422000 ,  437000 , 1.04        , 0.04            ,
  # Very high MW range
  "PS-630000"    ,  630000 ,  618000 ,  640000 , 1.04        , 0.05            ,
  "PS-1090000"   , 1090000 , 1070000 , 1100000 , 1.03        , 0.05            ,
  "PS-1870000"   , 1870000 , 1830000 , 1900000 , 1.04        , 0.06            ,
  "PS-3150000"   , 3150000 , 3080000 , 3200000 , 1.04        , 0.07
) |>
  mutate(
    polymer_type = "polystyrene",
    kit_name = "EasiVial PS-H Mix",
    # Mark-Houwink parameters for PS in THF at 35°C (well-established values)
    k_value = 0.000141, # mL/g
    a_value = 0.700,
    # dn/dc for PS in THF at 633 nm
    dn_dc = 0.185
  )

# ==============================================================================
# PMMA Standards
# For comparison/validation and demonstrating polymer-specific calibration
# ==============================================================================

pmma_standards <- tribble(
  ~standard_name , ~mp     , ~mn     , ~mw     , ~dispersity , ~mp_uncertainty ,
  "PMMA-602"     ,     602 ,     590 ,     612 , 1.04        , 0.08            ,
  "PMMA-1520"    ,    1520 ,    1490 ,    1550 , 1.04        , 0.06            ,
  "PMMA-4920"    ,    4920 ,    4800 ,    5010 , 1.04        , 0.04            ,
  "PMMA-10500"   ,   10500 ,   10300 ,   10700 , 1.04        , 0.04            ,
  "PMMA-30300"   ,   30300 ,   29700 ,   30800 , 1.04        , 0.03            ,
  "PMMA-67700"   ,   67700 ,   66400 ,   68800 , 1.04        , 0.03            ,
  "PMMA-137000"  ,  137000 ,  134000 ,  139000 , 1.04        , 0.04            ,
  "PMMA-342000"  ,  342000 ,  335000 ,  348000 , 1.04        , 0.04            ,
  "PMMA-675000"  ,  675000 ,  662000 ,  686000 , 1.04        , 0.05            ,
  "PMMA-1190000" , 1190000 , 1170000 , 1210000 , 1.03        , 0.06
) |>
  mutate(
    polymer_type = "pmma",
    kit_name = "ReadyCal PMMA Kit",
    # Mark-Houwink parameters for PMMA in THF at 35°C
    k_value = 0.000128, # mL/g
    a_value = 0.690,
    # dn/dc for PMMA in THF at 633 nm
    dn_dc = 0.084
  )

# ==============================================================================
# Validate Input Data
# ==============================================================================

validate_standards <- function(df, name) {
  errors <- character()

  # Check all MW values are positive
  if (any(df$mp <= 0) || any(df$mn <= 0) || any(df$mw <= 0)) {
    errors <- c(errors, paste0(name, ": Found non-positive MW values"))
  }

  # Check Mp is between Mn and Mw (physical requirement)
  if (any(df$mp < df$mn | df$mp > df$mw)) {
    bad_rows <- which(df$mp < df$mn | df$mp > df$mw)
    errors <- c(
      errors,
      paste0(
        name,
        ": Mp not between Mn and Mw for ",
        paste(df$standard_name[bad_rows], collapse = ", ")
      )
    )
  }

  # Check dispersity is approximately Mw/Mn
  calc_dispersity <- df$mw / df$mn
  if (any(abs(calc_dispersity - df$dispersity) > 0.02)) {
    bad_rows <- which(abs(calc_dispersity - df$dispersity) > 0.02)
    errors <- c(
      errors,
      paste0(
        name,
        ": Dispersity mismatch for ",
        paste(df$standard_name[bad_rows], collapse = ", ")
      )
    )
  }

  # Check uncertainty is in valid range
  if (any(df$mp_uncertainty <= 0 | df$mp_uncertainty > 1)) {
    errors <- c(
      errors,
      paste0(name, ": mp_uncertainty should be between 0 and 1")
    )
  }

  if (length(errors) > 0) {
    stop(paste(errors, collapse = "\n"), call. = FALSE)
  }

  cat(name, ": All", nrow(df), "standards passed validation\n")
}

validate_standards(ps_standards, "PS Standards")
validate_standards(pmma_standards, "PMMA Standards")

# ==============================================================================
# Generate Retention Times
# Using linear calibration model: log(Mp) = A + B*Vr
# This simulates a typical PLgel Mixed-C column set
# ==============================================================================

# Calibration curve parameters (simulating real column behavior)
# Retention volume range: ~12 to ~22 mL (dual 300x7.5mm columns)
# Exclusion limit: ~3,000,000 Da; Total permeation: ~100 Da

calculate_retention <- function(log_mp, flow_rate = 1.0) {
  # Input validation

  if (length(flow_rate) != 1 || flow_rate <= 0) {
    stop(
      "flow_rate must be a single positive value, got: ",
      flow_rate,
      call. = FALSE
    )
  }
  if (anyNA(log_mp) || !all(is.finite(log_mp))) {
    stop("log_mp contains NA or non-finite values", call. = FALSE)
  }

  # Polynomial coefficients for log(Mp) vs retention volume
  # Using linear approximation: log(Mp) = A + B*Vr
  # (Full model could include quadratic term, but linear is sufficient for
  # generating realistic synthetic data)
  A <- 11.5 # Intercept
  B <- -0.45 # Linear term (negative because higher MW elutes earlier)

  # Solve for retention volume given log(Mp)
  Vr_approx <- (log_mp - A) / B

  # Validate intermediate result
  if (any(Vr_approx <= 0)) {
    stop(
      "Calculated negative retention volumes. Input log_mp values may be ",
      "outside the valid calibration range.",
      call. = FALSE
    )
  }

  # Add small random variation (±0.5% CV) to simulate real injection variability
  Vr_with_noise <- Vr_approx * (1 + rnorm(length(Vr_approx), 0, 0.005))

  # Convert to retention time (minutes)
  retention_time <- Vr_with_noise / flow_rate

  # Also return retention volume
  list(
    retention_time = round(retention_time, 3),
    retention_volume = round(Vr_with_noise, 3)
  )
}

# Calculate retention data for PS standards
ps_retention <- calculate_retention(log10(ps_standards$mp))
ps_standards <- ps_standards |>
  mutate(
    retention_time = ps_retention$retention_time,
    retention_volume = ps_retention$retention_volume
  )

# Calculate retention data for PMMA standards
# PMMA has slightly smaller hydrodynamic volume than PS at same MW
# So it elutes ~0.2-0.5 mL later
pmma_retention <- calculate_retention(log10(pmma_standards$mp) - 0.05)
pmma_standards <- pmma_standards |>
  mutate(
    retention_time = pmma_retention$retention_time,
    retention_volume = pmma_retention$retention_volume
  )

# ==============================================================================
# Combine and Finalize Dataset
# ==============================================================================

sec_calibration_standards <- bind_rows(ps_standards, pmma_standards) |>
  mutate(
    # Add log values for convenience
    log_mp = log10(mp),
    log_mw = log10(mw),
    log_mn = log10(mn),
    # Calculate intrinsic viscosity from Mark-Houwink
    intrinsic_viscosity = k_value * mw^a_value,
    # Calculate hydrodynamic volume (M * [η])
    log_hydrodynamic_vol = log10(mw * intrinsic_viscosity),
    # Add notes
    notes = case_when(
      mp < 500 ~ "Oligomer range - verify column resolution",
      mp > 2000000 ~ "Near exclusion limit - verify peak shape",
      TRUE ~ NA_character_
    )
  ) |>
  select(
    # Identification
    standard_name,
    polymer_type,
    kit_name,
    # Certified MW values
    mp,
    mn,
    mw,
    dispersity,
    mp_uncertainty,
    # Pre-calculated log values
    log_mp,
    log_mw,
    log_mn,
    # Retention data
    retention_time,
    retention_volume,
    # Universal calibration parameters
    k_value,
    a_value,
    intrinsic_viscosity,
    log_hydrodynamic_vol,
    # Physical properties
    dn_dc,
    # Notes
    notes
  ) |>
  arrange(polymer_type, desc(mp))

# ==============================================================================
# Validate Combined Dataset
# ==============================================================================

validate_combined_data <- function(df) {
  errors <- character()

  # Check for NA/NaN/Inf in critical numeric columns
  numeric_cols <- c(
    "mp",
    "mn",
    "mw",
    "log_mp",
    "log_mw",
    "log_mn",
    "retention_time",
    "retention_volume",
    "intrinsic_viscosity",
    "log_hydrodynamic_vol"
  )
  for (col in numeric_cols) {
    if (any(is.na(df[[col]]) | !is.finite(df[[col]]))) {
      bad_count <- sum(is.na(df[[col]]) | !is.finite(df[[col]]))
      errors <- c(
        errors,
        paste0("Column '", col, "' has ", bad_count, " NA/NaN/Inf values")
      )
    }
  }

  # Verify intrinsic viscosity is positive and reasonable
  if (any(df$intrinsic_viscosity <= 0)) {
    errors <- c(errors, "Found non-positive intrinsic viscosity values")
  }

  # Verify monotonic relationship: higher MW = lower retention time (in SEC)
  for (polymer in unique(df$polymer_type)) {
    subset <- df[df$polymer_type == polymer, ]
    subset <- subset[order(subset$mp, decreasing = TRUE), ]
    if (!all(diff(subset$retention_time) >= -0.01)) {
      # Small tolerance for noise
      errors <- c(
        errors,
        paste0(
          polymer,
          ": Retention times not monotonically increasing with decreasing MW"
        )
      )
    }
  }

  # Verify row count
  if (nrow(df) != 26) {
    errors <- c(errors, paste0("Expected 26 standards, got ", nrow(df)))
  }

  if (length(errors) > 0) {
    stop(
      paste(c("Combined dataset validation failed:", errors), collapse = "\n"),
      call. = FALSE
    )
  }

  cat("Combined dataset: All", nrow(df), "standards passed validation\n")
}

validate_combined_data(sec_calibration_standards)

# ==============================================================================
# Create Quick-Reference Subsets
# ==============================================================================
# These are commonly used subsets for specific applications

# PS-only standards for conventional calibration (most common use case)
sec_ps_standards <- sec_calibration_standards |>
  filter(polymer_type == "polystyrene") |>
  select(
    standard_name,
    mp,
    log_mp,
    retention_time,
    retention_volume,
    mn,
    mw,
    dispersity,
    mp_uncertainty,
    k_value,
    a_value,
    dn_dc
  )

# PMMA-only standards
sec_pmma_standards <- sec_calibration_standards |>
  filter(polymer_type == "pmma") |>
  select(
    standard_name,
    mp,
    log_mp,
    retention_time,
    retention_volume,
    mn,
    mw,
    dispersity,
    mp_uncertainty,
    k_value,
    a_value,
    dn_dc
  )

# ==============================================================================
# Save Datasets
# ==============================================================================

tryCatch(
  {
    usethis::use_data(sec_calibration_standards, overwrite = TRUE)
    usethis::use_data(sec_ps_standards, overwrite = TRUE)
    usethis::use_data(sec_pmma_standards, overwrite = TRUE)

    # Verify files were created
    expected_files <- c(
      "data/sec_calibration_standards.rda",
      "data/sec_ps_standards.rda",
      "data/sec_pmma_standards.rda"
    )
    missing <- expected_files[!file.exists(expected_files)]
    if (length(missing) > 0) {
      stop("Expected files not created: ", paste(missing, collapse = ", "))
    }
    cat("\nAll 3 datasets saved successfully to data/\n")
  },
  error = function(e) {
    stop("Failed to save datasets: ", conditionMessage(e), call. = FALSE)
  }
)

# ==============================================================================
# Summary Output
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("SEC Calibration Standards Datasets Generated\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")

cat("\n=== sec_calibration_standards (full dataset) ===\n")
cat("Total standards:", nrow(sec_calibration_standards), "\n")
cat(
  "Polymer types:",
  paste(unique(sec_calibration_standards$polymer_type), collapse = ", "),
  "\n"
)
cat(
  "MW range:",
  format(min(sec_calibration_standards$mp), big.mark = ","),
  "to",
  format(max(sec_calibration_standards$mp), big.mark = ","),
  "Da\n"
)
cat(
  "Retention time range:",
  round(min(sec_calibration_standards$retention_time), 1),
  "to",
  round(max(sec_calibration_standards$retention_time), 1),
  "min\n"
)

cat("\n=== sec_ps_standards (polystyrene only) ===\n")
cat("Standards:", nrow(sec_ps_standards), "\n")
cat(
  "MW range:",
  format(min(sec_ps_standards$mp), big.mark = ","),
  "to",
  format(max(sec_ps_standards$mp), big.mark = ","),
  "Da\n"
)
cat("For use with: step_sec_conventional_cal()\n")

cat("\n=== sec_pmma_standards (PMMA only) ===\n")
cat("Standards:", nrow(sec_pmma_standards), "\n")
cat(
  "MW range:",
  format(min(sec_pmma_standards$mp), big.mark = ","),
  "to",
  format(max(sec_pmma_standards$mp), big.mark = ","),
  "Da\n"
)

cat("\n=== Column Names ===\n")
print(names(sec_calibration_standards))

cat("\n=== Preview: PS Standards ===\n")
print(
  sec_ps_standards |>
    select(standard_name, mp, log_mp, retention_time) |>
    head(8)
)

cat("\n=== Usage Example ===\n")
cat(
  '
library(measure.sec)
data(sec_ps_standards)

# For step_sec_conventional_cal, use retention_time and log_mp:
standards_for_cal <- sec_ps_standards |>
  select(retention = retention_time, log_mw = log_mp)

rec <- recipe(~., data = my_data) |>
  step_sec_conventional_cal(
    standards = standards_for_cal,
    fit_type = "cubic"
  )
'
)
