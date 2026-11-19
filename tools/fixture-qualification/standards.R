#!/usr/bin/env Rscript

# Qualify the raw SEC standards fixture without regenerating package data.
#
# Usage:
#   Rscript tools/fixture-qualification/standards.R <measure.sec> <output>
#
# The generator file contains calls to usethis::use_data() at its end.  This
# script evaluates only the helper definitions and the standards generator,
# then writes all evidence to the caller-supplied output directory.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript standards.R <measure.sec repository> <output directory>")
}

repo <- normalizePath(args[[1]], mustWork = TRUE)
output <- normalizePath(args[[2]], mustWork = FALSE)
dir.create(output, recursive = TRUE, showWarnings = FALSE)

generator_file <- file.path(repo, "data-raw", "generate_raw_chromatograms.R")
data_file <- file.path(repo, "data", "sec_raw_standards.rda")
if (!file.exists(generator_file) || !file.exists(data_file)) {
  stop(
    "Expected measure.sec generator and data files were not found.",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(tibble)
})

# Evaluate only definitions needed for Dataset 1. In particular, do not source
# the generator: its final expressions overwrite the shipped .rda files.
generator <- new.env(parent = globalenv())
allowed <- c(
  "flow_rate",
  "time_step",
  "mw_to_time",
  "generate_peak",
  "add_baseline",
  "add_noise",
  "generate_raw_ps_standards"
)
expressions <- parse(file = generator_file, keep.source = FALSE)
for (expression in expressions) {
  if (!is.call(expression) || !identical(expression[[1L]], as.name("<-"))) {
    next
  }
  name <- expression[[2L]]
  if (is.symbol(name) && as.character(name) %in% allowed) {
    eval(expression, envir = generator)
  }
}
missing_definitions <- setdiff(allowed, ls(generator))
if (length(missing_definitions) > 0L) {
  stop(
    "Generator definitions were not found: ",
    paste(missing_definitions, collapse = ", "),
    call. = FALSE
  )
}

# Capture the signal at each generator stage while preserving every numerical
# operation. The wrappers live in the generator's environment, so the
# generator's internal function resolves them just as it would the originals.
original_peak <- generator$generate_peak
original_baseline <- generator$add_baseline
original_noise <- generator$add_noise
stage <- new.env(parent = emptyenv())
stage$records <- list()
stage$current <- NULL
generator$generate_peak <- function(...) {
  value <- original_peak(...)
  stage$current <- list(shape = value)
  value
}
generator$add_baseline <- function(...) {
  value <- original_baseline(...)
  stage$current$drift <- value
  value
}
generator$add_noise <- function(...) {
  value <- original_noise(...)
  stage$current$recorded <- value
  stage$records[[length(stage$records) + 1L]] <- stage$current
  value
}

shipped_environment <- new.env(parent = emptyenv())
load(data_file, envir = shipped_environment)
shipped <- shipped_environment$sec_raw_standards

set.seed(2026)
generated <- generator$generate_raw_ps_standards()
if (!identical(generated, shipped)) {
  stop(
    "The standards generator no longer reproduces the shipped fixture exactly.",
    call. = FALSE
  )
}

if (length(stage$records) != length(unique(shipped$standard_name))) {
  stop(
    "The generator stage capture did not observe every standard.",
    call. = FALSE
  )
}

baseline_from_edges <- function(
  signal,
  left_fraction = 0.05,
  right_fraction = 0.05
) {
  index <- seq_along(signal)
  left <- head(index, floor(length(index) * left_fraction))
  right <- tail(index, floor(length(index) * right_fraction))
  left_mean <- mean(signal[left])
  right_mean <- mean(signal[right])
  left_center <- mean(left)
  right_center <- mean(right)
  left_mean +
    (index - left_center) *
      (right_mean - left_mean) /
      (right_center - left_center)
}

peak_time <- function(time, signal) time[[which.max(signal)]]
fit_summary <- function(name, time) {
  ordered <- order(summary_data$mp_Da)
  x <- time[ordered]
  y <- log10(summary_data$mp_Da[ordered])
  fit <- stats::lm(y ~ stats::poly(x, 3L, raw = TRUE))
  residual <- stats::residuals(fit)
  grid <- seq(min(x), max(x), length.out = 200L)
  fitted_grid <- stats::predict(fit, newdata = data.frame(x = grid))
  data.frame(
    stage = name,
    mass_min_Da = min(summary_data$mp_Da),
    mass_max_Da = max(summary_data$mp_Da),
    retention_min_min = min(x),
    retention_max_min = max(x),
    adjacent_reversals = sum(diff(x) > 0),
    increasing_fit_intervals_of_199 = sum(diff(fitted_grid) > 0),
    max_percent_over_assignment = max(100 * expm1(-residual * log(10))),
    r_squared = 1 - sum(residual^2) / sum((y - mean(y))^2),
    rmse_log10_Da = sqrt(mean(residual^2)),
    row.names = NULL
  )
}

standard_names <- unique(shipped$standard_name)
records <- lapply(seq_along(standard_names), function(index) {
  name <- standard_names[[index]]
  data <- shipped[shipped$standard_name == name, , drop = FALSE]
  capture <- stage$records[[index]]
  time <- data$time_min
  corrected <- data$ri_mv - baseline_from_edges(data$ri_mv)
  edge_count <- floor(length(time) * 0.05)
  right_edge <- tail(seq_along(time), edge_count)
  noise <- capture$recorded - capture$drift
  data.frame(
    standard_name = name,
    mp_Da = data$mp[[1L]],
    dispersity = data$dispersity[[1L]],
    n_points = length(time),
    time_min = min(time),
    time_max = max(time),
    nominal_peak_min = generator$mw_to_time(data$mp[[1L]]),
    shape_peak_min = peak_time(time, capture$shape),
    drift_peak_min = peak_time(time, capture$drift),
    recorded_peak_min = peak_time(time, data$ri_mv),
    corrected_peak_min = peak_time(time, corrected),
    shape_right_edge_fraction = mean(capture$shape[right_edge]) /
      max(capture$shape),
    drift_right_edge_fraction = mean(capture$drift[right_edge]) /
      max(capture$drift),
    noise_sd_mV = stats::sd(noise),
    noise_rmse_mV = sqrt(mean(noise^2)),
    finite_recorded = sum(is.finite(data$ri_mv)),
    row.names = NULL
  )
})
summary_data <- do.call(rbind, records)

stage_times <- list(
  nominal = summary_data$nominal_peak_min,
  shape = summary_data$shape_peak_min,
  drift = summary_data$drift_peak_min,
  recorded = summary_data$recorded_peak_min,
  corrected = summary_data$corrected_peak_min
)
fits <- do.call(
  rbind,
  Map(fit_summary, names(stage_times), stage_times)
)

# The qualification is intentionally strict about the fixture's declared
# behavior. Noiseless generator stages must preserve the SEC direction; the
# noisy shipped signal is measured and reported, not silently repaired.
stopifnot(
  nrow(summary_data) == 12L,
  all(summary_data$finite_recorded == summary_data$n_points),
  fits$adjacent_reversals[fits$stage == "nominal"] == 0,
  fits$adjacent_reversals[fits$stage == "shape"] == 0,
  fits$adjacent_reversals[fits$stage == "drift"] == 0,
  all(summary_data$noise_sd_mV > 0),
  all(summary_data$shape_right_edge_fraction > 0),
  all(summary_data$shape_right_edge_fraction < 1)
)

write.csv(
  summary_data,
  file.path(output, "standards-stage-summary.csv"),
  row.names = FALSE
)
write.csv(
  fits,
  file.path(output, "standards-fit-summary.csv"),
  row.names = FALSE
)

pdf(file.path(output, "standards-qualification.pdf"), width = 8.1, height = 5.3)
par(mar = c(4, 4, 3, 1))
plot(
  summary_data$mp_Da,
  summary_data$nominal_peak_min,
  log = "x",
  type = "l",
  lwd = 2,
  xlab = "Assigned peak molecular weight (Da)",
  ylab = "Peak time (min)",
  ylim = range(unlist(stage_times)),
  main = "Generator stages and observed standard maxima"
)
points(summary_data$mp_Da, summary_data$shape_peak_min, pch = 1)
points(summary_data$mp_Da, summary_data$drift_peak_min, pch = 2)
points(summary_data$mp_Da, summary_data$recorded_peak_min, pch = 16)
points(summary_data$mp_Da, summary_data$corrected_peak_min, pch = 4)
legend(
  "topright",
  c("Nominal", "Noiseless shape", "After drift", "Recorded", "Edge-corrected"),
  lty = c(1, NA, NA, NA, NA),
  pch = c(NA, 1, 2, 16, 4),
  bty = "n"
)
dev.off()

cat("Reproduced shipped fixture exactly: yes\n")
cat("Standards:", nrow(summary_data), "\n")
cat(
  "Mass range (Da):",
  min(summary_data$mp_Da),
  "to",
  max(summary_data$mp_Da),
  "\n"
)
cat("Observed peak reversals:\n")
print(
  fits[, c("stage", "adjacent_reversals", "r_squared", "rmse_log10_Da")],
  row.names = FALSE
)
