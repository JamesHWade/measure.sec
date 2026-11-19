# Usage: Rscript tools/fixture-qualification/delays.R REPO OUTPUT
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
repo <- normalizePath(args[1])
out <- args[2]
dir.create(out, showWarnings = FALSE, recursive = TRUE)

# Stop at the requested object. The generator's final use_data() calls write
# package files and must not be executed as part of a qualification run.
replay <- function(script, last) {
  env <- new.env(parent = globalenv())
  for (expr in parse(file.path(repo, "data-raw", script))) {
    head <- expr[[1]]
    if (
      !is.symbol(head) ||
        !as.character(head) %in% c("<-", "library", "set.seed")
    ) {
      stop("Unexpected top-level expression before ", last, " in ", script)
    }
    eval(expr, env)
    if (exists(last, env, inherits = FALSE)) return(env)
  }
  stop("Did not generate ", last)
}
old <- replay("generate_datasets.R", "sec_triple_detect")
raw <- replay("generate_raw_chromatograms.R", "sec_raw_multidetector")
stopifnot(
  old$delay_convention == "legacy_earlier_v1",
  raw$delay_convention == "target_later_v1"
)
compare <- function(env, name) {
  src <- new.env()
  load(file.path(repo, "data", paste0(name, ".rda")), envir = src)
  # Unrounded legacy traces can differ at machine precision between libm builds.
  actual <- env[[name]]
  expected <- src[[name]]
  stopifnot(!is.null(expected), identical(names(expected), names(actual)))
  stopifnot(isTRUE(all.equal(actual, expected, tolerance = 1e-13)))
  data.frame(
    fixture = name,
    rows = nrow(expected),
    exactly_identical = identical(actual, expected)
  )
}
fixtures <- rbind(
  compare(old, "sec_triple_detect"),
  compare(raw, "sec_raw_standards"),
  compare(raw, "sec_raw_unknowns"),
  compare(raw, "sec_raw_multidetector")
)
write.csv(
  fixtures,
  file.path(out, "delay-fixture-replay.csv"),
  row.names = FALSE
)
print(fixtures)

# At the shipped flow of 1 mL/min, a missing volume-to-time conversion would
# pass the replay. Capture the generator's actual time arguments at another flow.
flow <- 0.5
raw$flow_rate <- flow
observations <- new.env()
observations$delays <- numeric()
construct <- raw$apply_delay
raw$apply_delay <- function(signal, time, delay_time) {
  observations$delays <- c(observations$delays, delay_time)
  construct(signal, time, delay_time)
}
alternate <- raw$generate_raw_multidetector()
declared <- unique(alternate[, c("sample_id", "delay_uv_ml", "delay_mals_ml")])
expected_delays <- as.vector(t(as.matrix(declared[, c(
  "delay_uv_ml",
  "delay_mals_ml"
)]))) /
  flow
stopifnot(
  identical(observations$delays, expected_delays),
  all(declared$delay_uv_ml == -0.08),
  all(declared$delay_mals_ml == 0.18)
)

# A binary-exact grid makes the known shifts integer sample displacements.
# Indexing (not another call to approx) defines the reference and its missing
# boundaries. Noisy cases translate the same realization, not independent noise.
time <- seq(0, 6, by = 1 / 64)
impulse <- as.numeric(time == 3)
smooth <- exp(-0.5 * ((time - 3) / 0.45)^2)
set.seed(88)
signals <- list(
  impulse = impulse,
  broad = smooth,
  noisy = smooth + rnorm(length(time), sd = 0.02)
)
cases <- expand.grid(
  signal = names(signals),
  offset = c(-16L, 16L),
  axis = c("min", "mL"),
  stringsAsFactors = FALSE
)
trace <- function(x, y) {
  measure::new_measure_list(list(measure::new_measure_tbl(x, y)))
}
run <- function(i) {
  case <- cases[i, ]
  k <- case$offset
  delta <- k / 64
  y <- signals[[case$signal]]
  idx <- seq_along(y)
  target <- construct(y, time, delta)
  stopifnot(identical(target, y[pmin(length(y), pmax(1L, idx - k))]))
  # The archived example displaces the same signal in the opposite direction.
  stopifnot(identical(
    old$apply_delay_legacy_v1(y, time, delta),
    y[pmin(length(y), pmax(1L, idx + k))]
  ))
  x <- if (case$axis == "min") time else time * flow
  data <- tibble::tibble(reference = trace(x, y), target = trace(x, target))
  recipe <- measure.sec::step_sec_detector_delay(
    recipes::recipe(data),
    reference = "reference",
    targets = "target",
    flow_rate = flow,
    elution_units = case$axis,
    delay_times = if (case$axis == "mL") c(target = delta),
    delay_volumes = if (case$axis == "min") c(target = delta * flow),
    method = "interpolate"
  )
  aligned <- recipes::bake(recipes::prep(recipe), new_data = data)$target[[1]]
  supported <- idx + k >= 1L & idx + k <= length(y)
  stopifnot(
    identical(aligned$location, x),
    identical(is.finite(aligned$value), supported),
    isTRUE(all.equal(aligned$value[supported], y[supported], tolerance = 1e-13))
  )
  data.frame(
    case,
    delay_minutes = delta,
    delay_mL = delta * flow,
    finite = sum(supported),
    max_abs_error = max(abs(aligned$value[supported] - y[supported]))
  )
}
results <- do.call(rbind, lapply(seq_len(nrow(cases)), run))
write.csv(results, file.path(out, "delay-controls.csv"), row.names = FALSE)
print(results)
