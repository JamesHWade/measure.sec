#!/usr/bin/env Rscript

# Qualify the packaged SEC unknowns without changing the fixture.
#
# Usage:
#   Rscript tools/fixture-qualification/unknowns.R <measure.sec> <output-dir>
#
# The script deliberately keeps the fitted calibration and the generator's
# nominal retention mapping separate. The former is the observable recipe
# result. The latter is a sensitivity calculation on the same supported
# interval and is not an experimental calibration.

options(warn = 1) # Print the calibration diagnostics alongside each recipe.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop(
    "Usage: Rscript tools/fixture-qualification/unknowns.R ",
    "<measure.sec-repository> <output-directory>",
    call. = FALSE
  )
}

repo <- normalizePath(args[[1]], mustWork = TRUE)
output <- normalizePath(args[[2]], mustWork = FALSE)
if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
}

suppressPackageStartupMessages({
  library(measure.sec)
  library(measure)
  library(recipes)
})

installed_fixture <- function(name) {
  src <- new.env()
  installed <- new.env()
  load(file.path(repo, "data", paste0(name, ".rda")), envir = src)
  utils::data(list = name, package = "measure.sec", envir = installed)
  stopifnot(identical(src[[name]], installed[[name]]), !is.null(src[[name]]))
  installed[[name]]
}
unknowns <- installed_fixture("sec_raw_unknowns")
standards_raw <- installed_fixture("sec_raw_standards")

# The tutorial estimates each standard's apex after subtracting the default
# straight-line edge baseline.  This small independent implementation keeps
# the calibration input auditable instead of reading a fitted recipe step.
independent_baseline <- function(data) {
  index <- seq_len(nrow(data))
  edge_size <- max(1L, floor(0.05 * length(index)))
  left <- head(index, edge_size)
  right <- tail(index, edge_size)
  signal <- data$ri_mv
  baseline <- mean(signal[left]) +
    (index - mean(left)) *
      (mean(signal[right]) - mean(signal[left])) /
      (mean(right) - mean(left))
  list(location = data$time_min, signal = signal - baseline)
}

standard_apex <- function(data) {
  corrected <- independent_baseline(data)
  apex <- which.max(corrected$signal)
  data.frame(
    standard_name = data$standard_name[[1]],
    mp = data$mp[[1]],
    log_mw = data$log_mp[[1]],
    retention = corrected$location[[apex]],
    stringsAsFactors = FALSE
  )
}

standard_groups <- split(standards_raw, standards_raw$standard_name)
standard_apices <- do.call(rbind, lapply(standard_groups, standard_apex))
standard_apices <- standard_apices[order(standard_apices$retention), ]
row.names(standard_apices) <- NULL

calibration <- stats::lm(
  log_mw ~ stats::poly(retention, 3),
  data = standard_apices
)
calibration_range <- range(standard_apices$retention)
# All six unknowns use the standards' sampling grid. Work only on actual
# readings; fail if a window endpoint ever requires interpolation.
slice_weights <- function(x, signal) {
  dx <- diff(x)
  stopifnot(length(dx) >= 1L, all(dx > 0), all(is.finite(signal)))
  pmax(signal, 0) * (c(0, dx) + c(dx, 0)) / 2
}

moments <- function(weight, mass) {
  stopifnot(all(is.finite(mass) & mass > 0), sum(weight) > 0)
  total <- sum(weight)
  c(
    mn = total / sum(weight / mass),
    mw = sum(weight * mass) / total,
    mz = sum(weight * mass^2) / sum(weight * mass)
  )
}

recipe_for_unknown <- function(data) {
  recipes::recipe(ri_mv + time_min ~ sample_id, data = data) |>
    recipes::update_role(sample_id, new_role = "id") |>
    measure::step_measure_input_long(
      ri_mv,
      location = dplyr::vars(time_min),
      col_name = "ri"
    ) |>
    step_sec_baseline(measures = "ri") |>
    step_sec_conventional_cal(
      standards = standard_apices[c("retention", "log_mw")],
      fit_type = "cubic",
      extrapolation = "none"
    ) |>
    step_sec_mw_averages(integration_range = calibration_range)
}

qualify_unknown <- function(data) {
  prepped <- recipes::prep(recipe_for_unknown(data))
  baked <- recipes::bake(prepped, new_data = NULL)
  signal <- baked$ri[[1]]
  location <- signal$location
  corrected <- signal$value
  independent <- independent_baseline(data)
  if (
    !isTRUE(all.equal(location, independent$location, tolerance = 1e-12)) ||
      !isTRUE(all.equal(corrected, independent$signal, tolerance = 1e-10))
  ) {
    stop(
      "The public baseline step differs from the independent edge-line control.",
      call. = FALSE
    )
  }
  x <- independent$location
  y <- independent$signal
  keep <- x >= calibration_range[1] & x <= calibration_range[2]
  stopifnot(identical(range(x[keep]), calibration_range))
  weights <- slice_weights(x[keep], y[keep])
  fitted_mass <- 10^stats::predict(calibration, data.frame(retention = x[keep]))
  supported <- moments(weights, fitted_mass)
  nominal_supported <- moments(weights, 10^((19 - x[keep]) / 0.75))
  area_ratio <- sum(weights) / sum(slice_weights(x, y))
  recipe_values <- unlist(baked[c("mw_mn", "mw_mw", "mw_mz")])
  if (
    !isTRUE(all.equal(
      unname(recipe_values),
      unname(supported),
      tolerance = 1e-8
    ))
  ) {
    stop(
      "Independent integration disagrees with the installed recipe for ",
      data$sample_id[[1L]]
    )
  }
  data.frame(
    sample_id = data$sample_id[[1]],
    assigned_mn = data$true_mn[[1]],
    assigned_mw = data$true_mw[[1]],
    assigned_mz = data$true_mz[[1]],
    assigned_dispersity = data$true_dispersity[[1]],
    recipe_mn = recipe_values[["mw_mn"]],
    recipe_mw = recipe_values[["mw_mw"]],
    recipe_mz = recipe_values[["mw_mz"]],
    direct_supported_mn = supported[["mn"]],
    direct_supported_mw = supported[["mw"]],
    direct_supported_mz = supported[["mz"]],
    nominal_supported_mn = nominal_supported[["mn"]],
    nominal_supported_mw = nominal_supported[["mw"]],
    nominal_supported_mz = nominal_supported[["mz"]],
    supported_positive_signal_area_fraction = area_ratio,
    independent_baseline_max_abs_difference = max(
      abs(corrected - independent$signal)
    ),
    stringsAsFactors = FALSE
  )
}

unknown_groups <- split(unknowns, unknowns$sample_id)
results <- do.call(rbind, lapply(unknown_groups, qualify_unknown))
row.names(results) <- NULL
write.csv(
  results,
  file.path(output, "unknown-moment-qualification.csv"),
  row.names = FALSE
)

# A second control checks the arithmetic with a deliberately known distribution
# on the same direct calculator. Three parts of the signal have MW 5,000 and
# one part has MW 20,000, so the expected moments are analytic.
control_signal <- c(0, 3, 0, 0, 1, 0)
control_mass <- c(5000, 5000, 5000, 20000, 20000, 20000)
control <- moments(
  slice_weights(seq_along(control_signal), control_signal),
  control_mass
)
expected <- c(
  mn = 80000 / 13,
  mw = 8750,
  mz = 118750000 / 8750
)
if (
  !isTRUE(all.equal(control[c("mn", "mw", "mz")], expected, tolerance = 1e-10))
) {
  stop("Known-weight moment control failed.", call. = FALSE)
}
control_output <- data.frame(
  moment = names(expected),
  expected = unname(expected),
  direct_observed = unname(control[names(expected)]),
  stringsAsFactors = FALSE
)
control_data <- data.frame(id = "known-weight-control")
control_data$signal <- measure::new_measure_list(list(
  measure::new_measure_tbl(seq_along(control_signal), control_signal)
))
control_data$mass <- measure::new_measure_list(list(
  measure::new_measure_tbl(seq_along(control_mass), control_mass)
))
control_recipe <- recipes::recipe(~., data = control_data) |>
  step_sec_mw_averages(
    measures = "signal",
    mw_col = "mass",
    mw_scale = "absolute",
    output_cols = c("mn", "mw", "mz")
  )
control_baked <- recipes::bake(
  recipes::prep(control_recipe),
  new_data = NULL
)
public_observed <- unname(as.numeric(control_baked[
  1,
  c(
    "mw_mn",
    "mw_mw",
    "mw_mz"
  )
]))
if (!isTRUE(all.equal(public_observed, unname(expected), tolerance = 1e-10))) {
  stop("Public known-weight moment control failed.", call. = FALSE)
}
control_output$public_recipe_observed <- public_observed
write.csv(
  control_output,
  file.path(output, "known-weight-control.csv"),
  row.names = FALSE
)

# Keep the output useful without requiring a plotting package. The upper panel
# shows the reported Mw against the assigned metadata; the lower panel shows
# how much positive signal lies inside the fitted calibration range.
png(
  file.path(output, "unknown-moment-qualification.png"),
  width = 1200,
  height = 900,
  res = 140
)
par(mfrow = c(2, 1), mar = c(3.5, 4.5, 1, 1))
finite <- is.finite(results$assigned_mw)
parity_limits <- range(
  c(
    results$assigned_mw[finite],
    results$recipe_mw[finite]
  ),
  finite = TRUE
)
plot(
  results$assigned_mw[finite],
  results$recipe_mw[finite],
  log = "xy",
  xlim = parity_limits,
  ylim = parity_limits,
  pch = 19,
  col = "#0072B2",
  xlab = "Assigned Mw (Da)",
  ylab = "Recipe Mw (Da)"
)
abline(0, 1, lty = 2, col = "grey50")
par(mar = c(5.7, 4.5, 1, 1))
barplot(
  results$supported_positive_signal_area_fraction,
  names.arg = results$sample_id,
  las = 2,
  ylim = c(0, 1),
  ylab = "Fraction of positive area in range",
  col = "#A44200"
)
dev.off()

metadata <- data.frame(
  installed_package = find.package("measure.sec"),
  calibration_mapping = "cubic fit to independently baseline-corrected raw-standard apices",
  nominal_generator_mapping = "log10(MW) = (19 - retention_min) / 0.75; used only on the supported interval",
  supported_interval = paste(calibration_range, collapse = " to "),
  area_definition = paste(
    "sum(pmax(independent baseline-corrected signal, 0) * trapezoid widths)",
    "inside the fitted standard-apex interval divided by the same quantity",
    "over the complete observed trace"
  ),
  stringsAsFactors = FALSE
)
write.csv(
  metadata,
  file.path(output, "qualification-metadata.csv"),
  row.names = FALSE
)

cat("Qualified", nrow(results), "unknowns.\n")
cat("Calibration interval (min):", calibration_range, "\n")
print(
  results[, c(
    "sample_id",
    "recipe_mw",
    "assigned_mw",
    "nominal_supported_mw",
    "supported_positive_signal_area_fraction"
  )],
  row.names = FALSE
)
cat("Known-weight moment control passed.\n")
