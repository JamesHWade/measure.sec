#!/usr/bin/env Rscript

qualification_sha256 <- c(
  fixtures = "dbd3c3e15015f81264a83743b5989da1fe0e43712fda71e3a823ab46e98b6680",
  declaration = "26f2810af74c8a0e80eda4b0baaf29192f3d12e6c6adccf2d6eb68335a7e1316"
)

json_scalar <- function(x) {
  jsonlite::unbox(x)
}

json_array <- function(x, encode = json_scalar) {
  lapply(as.list(x), encode)
}

write_json <- function(x, path, pretty = TRUE) {
  encoded <- jsonlite::toJSON(
    x,
    auto_unbox = FALSE,
    null = "null",
    na = "null",
    digits = 17,
    pretty = pretty,
    force = TRUE
  )
  writeLines(as.character(encoded), con = path, useBytes = TRUE)
  invisible(path)
}

write_new_json <- function(x, path) {
  link_target <- Sys.readlink(path)
  has_link <- length(link_target) == 1L &&
    !is.na(link_target) &&
    nzchar(link_target)
  if (file.exists(path) || dir.exists(path) || has_link) {
    stop(
      "Refusing to overwrite an existing qualification file: ",
      path,
      call. = FALSE
    )
  }
  write_json(x, path)
}

decode_rscript_path <- function(path) {
  gsub("~+~", " ", path, fixed = TRUE)
}

normalize_input_path <- function(path) {
  normalizePath(
    path.expand(decode_rscript_path(path)),
    winslash = "/",
    mustWork = FALSE
  )
}

sha256_file <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

read_json <- function(path) {
  jsonlite::fromJSON(path, simplifyVector = FALSE)
}

json_number <- function(x) {
  if (is.null(x)) {
    return(NA_real_)
  }
  if (length(x) != 1L || !is.numeric(x)) {
    stop(
      "A frozen numeric fixture value is not a numeric scalar.",
      call. = FALSE
    )
  }
  as.numeric(x)
}

json_numeric_vector <- function(x) {
  if (length(x) == 0L) {
    return(double())
  }
  vapply(x, json_number, numeric(1), USE.NAMES = FALSE)
}

json_value <- function(x) {
  if (length(x) != 1L || !is.numeric(x)) {
    stop("A measurement value was not a numeric scalar.", call. = FALSE)
  }
  if (is.nan(x)) {
    return(list(state = json_scalar("nan"), number = NULL))
  }
  if (is.na(x)) {
    return(list(state = json_scalar("na"), number = NULL))
  }
  if (is.infinite(x)) {
    state <- if (x > 0) "+inf" else "-inf"
    return(list(state = json_scalar(state), number = NULL))
  }
  list(state = json_scalar("number"), number = json_scalar(as.numeric(x)))
}

curve_json <- function(curve) {
  if (!inherits(curve, "measure_tbl")) {
    stop("A baked SEC curve did not retain class measure_tbl.", call. = FALSE)
  }
  location <- curve[["location"]]
  value <- curve[["value"]]
  if (
    !is.numeric(location) ||
      !is.null(dim(location)) ||
      !all(is.finite(location))
  ) {
    stop(
      "A baked curve location must be a finite numeric vector.",
      call. = FALSE
    )
  }
  if (
    !is.numeric(value) ||
      !is.null(dim(value)) ||
      length(location) != length(value)
  ) {
    stop(
      "A baked curve value must be a matching numeric vector.",
      call. = FALSE
    )
  }
  list(
    location = json_array(as.numeric(location)),
    value = lapply(as.list(value), json_value)
  )
}

measure_column <- function(
  ids,
  trace_map,
  signal_field,
  axis_per_minute,
  reverse = FALSE,
  concentration_unit = NULL
) {
  curves <- lapply(ids, function(id) {
    trace <- trace_map[[id]]
    location <- json_numeric_vector(trace$time_min) * axis_per_minute
    value <- json_numeric_vector(trace[[signal_field]])
    if (isTRUE(reverse)) {
      location <- rev(location)
      value <- rev(value)
    }
    measure::new_measure_tbl(location = location, value = value)
  })
  column <- measure::new_measure_list(curves)
  if (!is.null(concentration_unit)) {
    attr(column, "sec_concentration_unit") <- concentration_unit
  }
  column
}

mass_data <- function(ids, fixture_map, axis_per_minute, reverse = FALSE) {
  ids <- unlist(ids, use.names = FALSE)
  data <- data.frame(fixture_id = ids, stringsAsFactors = FALSE)
  data$ri <- measure_column(
    ids,
    fixture_map,
    signal_field = "signal",
    axis_per_minute = axis_per_minute,
    reverse = reverse
  )
  data
}

delay_data <- function(
  fixture,
  axis_per_minute,
  reverse = FALSE,
  reference_location = NULL,
  target_location = NULL,
  target_signal = NULL,
  reference_signal = NULL
) {
  location <- json_numeric_vector(fixture$time_min) * axis_per_minute
  if (is.null(target_signal)) {
    target_signal <- json_numeric_vector(fixture$signal)
  }
  if (is.null(reference_signal)) {
    reference_signal <- json_numeric_vector(fixture$reference)
  }
  if (is.null(reference_location)) {
    reference_location <- location
  }
  if (is.null(target_location)) {
    target_location <- location
  }
  if (isTRUE(reverse)) {
    target_location <- rev(target_location)
    reference_location <- rev(reference_location)
    target_signal <- rev(target_signal)
    reference_signal <- rev(reference_signal)
  }
  reference <- measure::new_measure_tbl(
    location = reference_location,
    value = reference_signal
  )
  target <- measure::new_measure_tbl(
    location = target_location,
    value = target_signal
  )
  data <- data.frame(row_id = 1L)
  data$reference <- measure::new_measure_list(list(reference))
  data$target <- measure::new_measure_list(list(target))
  attr(data$target, "sec_concentration_unit") <- "mg/mL"
  data
}

recipe_for <- function(data) {
  recipes::recipe(~., data = data)
}

stable_step_id <- function(prefix, id) {
  paste0(prefix, "_", gsub("[^A-Za-z0-9_]+", "_", id))
}

set_step_field <- function(recipe, field, value, keep_null = FALSE) {
  index <- length(recipe$steps)
  step <- recipe$steps[[index]]
  if (keep_null && is.null(value)) {
    step[field] <- list(NULL)
  } else if (is.null(value)) {
    step[field] <- NULL
  } else {
    step[field] <- list(value)
  }
  recipe$steps[[index]] <- step
  unserialize(serialize(recipe, connection = NULL))
}

remove_step_fields <- function(recipe, fields) {
  index <- length(recipe$steps)
  step <- recipe$steps[[index]]
  for (field in fields) {
    step[field] <- NULL
  }
  recipe$steps[[index]] <- step
  unserialize(serialize(recipe, connection = NULL))
}

run_with_warnings <- function(expr, context, allow_clipping = FALSE) {
  warning_records <- new.env(parent = emptyenv())
  warning_records$items <- list()
  value <- withCallingHandlers(
    force(expr),
    warning = function(warning) {
      if (
        allow_clipping &&
          inherits(warning, "measure_sec_negative_signal_clipped")
      ) {
        warning_records$items[[length(warning_records$items) + 1L]] <- warning
        invokeRestart("muffleWarning")
      }
      message <- paste0(
        "Unexpected warning during ",
        context,
        ": ",
        conditionMessage(warning)
      )
      stop(structure(
        list(message = message, call = NULL, warning = warning),
        class = c("unexpected_qualification_warning", "error", "condition")
      ))
    }
  )
  list(value = value, warnings = warning_records$items)
}

clipping_diagnostics <- function(warnings, context) {
  if (length(warnings) != 1L) {
    stop(
      context,
      " must emit exactly one clipping condition; observed ",
      length(warnings),
      call. = FALSE
    )
  }
  condition <- warnings[[1L]]
  diagnostics <- condition$diagnostics
  required <- c(
    "row",
    "negative_points",
    "points",
    "positive_area",
    "negative_area",
    "excluded_fraction"
  )
  if (
    !inherits(diagnostics, "data.frame") ||
      !identical(names(diagnostics), required) ||
      nrow(diagnostics) != 4L
  ) {
    stop(
      context,
      " clipping condition did not carry the four current rows.",
      call. = FALSE
    )
  }
  list(
    column = as.character(condition$column),
    elution_units = as.character(condition$elution_units),
    diagnostics = lapply(seq_len(nrow(diagnostics)), function(row) {
      record <- diagnostics[row, , drop = FALSE]
      list(
        row = json_scalar(as.integer(record$row[[1L]])),
        negative_points = json_scalar(as.integer(record$negative_points[[1L]])),
        points = json_scalar(as.integer(record$points[[1L]])),
        positive_area = json_scalar(as.numeric(record$positive_area[[1L]])),
        negative_area = json_scalar(as.numeric(record$negative_area[[1L]])),
        excluded_fraction = json_scalar(as.numeric(record$excluded_fraction[[
          1L
        ]]))
      )
    }),
    source = diagnostics
  )
}

package_provenance <- function(package) {
  installed_path <- find.package(package, quiet = TRUE)
  if (length(installed_path) != 1L || !nzchar(installed_path)) {
    stop("Installed package not found: ", package, call. = FALSE)
  }
  installed_path <- normalizePath(
    installed_path,
    winslash = "/",
    mustWork = TRUE
  )
  metadata_path <- file.path(installed_path, "Meta", "package.rds")
  if (!file.exists(metadata_path)) {
    stop("Installed package lacks Meta/package.rds: ", package, call. = FALSE)
  }
  namespace <- loadNamespace(package)
  namespace_path <- normalizePath(
    getNamespaceInfo(namespace, "path"),
    winslash = "/",
    mustWork = TRUE
  )
  if (!identical(namespace_path, installed_path)) {
    stop(
      "Loaded namespace path differs from find.package() for ",
      package,
      ": ",
      namespace_path,
      " != ",
      installed_path,
      call. = FALSE
    )
  }
  list(
    path = json_scalar(installed_path),
    namespace_path = json_scalar(namespace_path),
    version = json_scalar(as.character(utils::packageVersion(package)))
  )
}

tool_checkout_head <- function(package_root) {
  if (
    !file.exists(file.path(package_root, ".git")) &&
      !dir.exists(file.path(package_root, ".git"))
  ) {
    return(NULL)
  }
  description_path <- file.path(package_root, "DESCRIPTION")
  if (!file.exists(description_path)) {
    return(NULL)
  }
  package_name <- tryCatch(
    read.dcf(description_path, fields = "Package")[[1L, 1L]],
    error = function(error) ""
  )
  if (!identical(package_name, "measure.sec")) {
    return(NULL)
  }
  top <- suppressWarnings(system2(
    "git",
    c("-C", shQuote(package_root), "rev-parse", "--show-toplevel"),
    stdout = TRUE,
    stderr = FALSE
  ))
  if (!is.null(attr(top, "status")) && attr(top, "status") != 0L) {
    return(NULL)
  }
  if (length(top) != 1L || !nzchar(top[[1L]])) {
    return(NULL)
  }
  top_path <- normalizePath(top[[1L]], winslash = "/", mustWork = TRUE)
  if (!identical(top_path, package_root)) {
    return(NULL)
  }
  head <- suppressWarnings(system2(
    "git",
    c("-C", shQuote(package_root), "rev-parse", "HEAD"),
    stdout = TRUE,
    stderr = FALSE
  ))
  if (!is.null(attr(head, "status")) && attr(head, "status") != 0L) {
    return(NULL)
  }
  if (length(head) != 1L || !nzchar(head[[1L]])) {
    return(NULL)
  }
  as.character(head[[1L]])
}

run_mass_case <- function(case, fixture_map) {
  training <- mass_data(
    case$training_order,
    fixture_map,
    axis_per_minute = case$axis_per_minute,
    reverse = case$reverse
  )
  bake_data <- mass_data(
    case$bake_order,
    fixture_map,
    axis_per_minute = case$axis_per_minute,
    reverse = case$reverse
  )
  selected_masses <- vapply(
    case$bake_order,
    function(id) {
      as.numeric(fixture_map[[id]]$injection_mass_mg)
    },
    numeric(1),
    USE.NAMES = FALSE
  )
  unique_masses <- unique(selected_masses)
  by_mass <- vector("list", length(unique_masses))
  names(by_mass) <- format(unique_masses, scientific = TRUE, digits = 17)
  reference_diagnostics <- NULL

  for (mass_index in seq_along(unique_masses)) {
    mass <- unique_masses[[mass_index]]
    recipe <- recipe_for(training)
    recipe <- measure.sec::step_sec_concentration(
      recipe,
      measures = "ri",
      injection_mass = mass,
      flow_rate = case$flow_ml_min,
      concentration_units = case$concentration_units,
      elution_units = case$elution_units,
      id = stable_step_id("qual_mass", paste0(case$id, "_", mass_index))
    )
    if (isTRUE(case$legacy)) {
      recipe <- remove_step_fields(recipe, "elution_units")
    }
    prepared <- run_with_warnings(
      recipes::prep(recipe, training = training),
      paste0(case$id, " prep at mass ", format(mass, digits = 17)),
      allow_clipping = TRUE
    )
    baked <- run_with_warnings(
      recipes::bake(prepared$value, new_data = bake_data),
      paste0(case$id, " current bake at mass ", format(mass, digits = 17)),
      allow_clipping = TRUE
    )
    condition <- clipping_diagnostics(
      baked$warnings,
      paste0(case$id, " current bake at mass ", format(mass, digits = 17))
    )
    if (
      !identical(condition$column, "ri") ||
        !identical(condition$elution_units, case$elution_units)
    ) {
      stop(
        "Clipping warning fields did not retain the requested column and units.",
        call. = FALSE
      )
    }
    if (is.null(reference_diagnostics)) {
      reference_diagnostics <- condition$source
    } else if (!identical(reference_diagnostics, condition$source)) {
      stop(
        "Raw clipping diagnostics changed with injected mass.",
        call. = FALSE
      )
    }
    tag <- attr(baked$value$ri, "sec_concentration_unit", exact = TRUE)
    if (!identical(tag, case$concentration_units)) {
      stop(
        "Mass-normalized curve lost the requested concentration unit tag.",
        call. = FALSE
      )
    }
    by_mass[[mass_index]] <- list(data = baked$value, condition = condition)
  }

  members <- lapply(seq_along(case$bake_order), function(row) {
    fixture_id <- case$bake_order[[row]]
    mass <- selected_masses[[row]]
    mass_index <- match(mass, unique_masses)
    baked <- by_mass[[mass_index]]$data
    list(
      fixture_id = json_scalar(fixture_id),
      row = json_scalar(as.integer(row)),
      injection_mass_mg = json_scalar(mass),
      curve = curve_json(baked$ri[[row]])
    )
  })

  list(
    id = json_scalar(case$id),
    concentration_unit = json_scalar(case$concentration_units),
    members = members,
    condition = list(
      column = json_scalar("ri"),
      elution_units = json_scalar(case$elution_units),
      diagnostics = lapply(seq_len(nrow(reference_diagnostics)), function(row) {
        record <- reference_diagnostics[row, , drop = FALSE]
        list(
          row = json_scalar(as.integer(record$row[[1L]])),
          negative_points = json_scalar(as.integer(record$negative_points[[
            1L
          ]])),
          points = json_scalar(as.integer(record$points[[1L]])),
          positive_area = json_scalar(as.numeric(record$positive_area[[1L]])),
          negative_area = json_scalar(as.numeric(record$negative_area[[1L]])),
          excluded_fraction = json_scalar(as.numeric(record$excluded_fraction[[
            1L
          ]]))
        )
      })
    )
  )
}

delay_recipe <- function(
  data,
  case,
  delay_value,
  delay_kind,
  legacy = FALSE,
  id = case$id
) {
  recipe <- recipe_for(data)
  args <- list(
    recipe = recipe,
    reference = "reference",
    targets = "target",
    flow_rate = case$flow_ml_min,
    method = case$method,
    elution_units = case$elution_units,
    id = stable_step_id("qual_delay", id)
  )
  if (identical(delay_kind, "time")) {
    if (legacy) {
      delay_volume <- delay_value * case$flow_ml_min
      args$delay_volumes <- stats::setNames(delay_volume, "target")
    } else {
      args$delay_times <- stats::setNames(delay_value, "target")
    }
  } else {
    args$delay_volumes <- stats::setNames(delay_value, "target")
  }
  recipe <- do.call(measure.sec::step_sec_detector_delay, args)
  if (legacy) {
    recipe <- remove_step_fields(recipe, c("delay_times", "elution_units"))
  }
  recipe
}

run_delay_case <- function(case, fixture_map) {
  data <- delay_data(
    fixture_map[[case$fixture]],
    axis_per_minute = case$axis_per_minute,
    reverse = case$reverse
  )
  recipe <- delay_recipe(
    data,
    case,
    delay_value = case$delay_value,
    delay_kind = case$delay_kind,
    legacy = case$legacy
  )
  prepared <- run_with_warnings(
    recipes::prep(recipe, training = data),
    paste0(case$id, " prep")
  )
  baked <- run_with_warnings(
    recipes::bake(prepared$value, new_data = data),
    paste0(case$id, " current bake")
  )$value
  tag <- attr(baked$target, "sec_concentration_unit", exact = TRUE)
  if (!identical(tag, "mg/mL")) {
    stop(
      "Detector delay did not retain the target concentration tag.",
      call. = FALSE
    )
  }
  list(
    id = json_scalar(case$id),
    target = curve_json(baked$target[[1L]]),
    reference = curve_json(baked$reference[[1L]]),
    concentration_unit = json_scalar(tag)
  )
}

run_noop_case <- function(case, noop_trace) {
  location <- json_numeric_vector(noop_trace$time_min)
  value <- json_numeric_vector(noop_trace$signal)
  curve <- measure::new_measure_tbl(location = location, value = value)
  data <- data.frame(row_id = 1L)
  data$ri <- measure::new_measure_list(list(curve))
  attr(data$ri, "sec_concentration_unit") <- noop_trace$concentration_unit
  recipe <- recipe_for(data)
  recipe <- measure.sec::step_sec_concentration(
    recipe,
    measures = "ri",
    flow_rate = case$flow_ml_min,
    normalize_to_mass = FALSE,
    elution_units = case$elution_units,
    id = stable_step_id("qual_noop", case$id)
  )
  prepared <- run_with_warnings(
    recipes::prep(recipe, training = data),
    paste0(case$id, " prep")
  )
  baked <- run_with_warnings(
    recipes::bake(prepared$value, new_data = data),
    paste0(case$id, " bake")
  )$value
  tag <- attr(baked$ri, "sec_concentration_unit", exact = TRUE)
  if (!identical(tag, noop_trace$concentration_unit)) {
    stop(
      "No-op normalization did not preserve the existing concentration tag.",
      call. = FALSE
    )
  }
  list(
    id = json_scalar(case$id),
    curve = curve_json(baked$ri[[1L]]),
    concentration_unit = json_scalar(tag)
  )
}

run_lals_case <- function(case, lals_fixture) {
  location <- json_numeric_vector(lals_fixture$time_min) * case$axis_per_minute
  ri <- measure::new_measure_tbl(
    location = location,
    value = json_numeric_vector(lals_fixture$ri_signal)
  )
  lals <- measure::new_measure_tbl(
    location = location,
    value = json_numeric_vector(lals_fixture$scattering_signal)
  )
  data <- data.frame(row_id = 1L)
  data$ri <- measure::new_measure_list(list(ri))
  data$lals <- measure::new_measure_list(list(lals))

  recipe <- recipe_for(data)
  recipe <- measure.sec::step_sec_concentration(
    recipe,
    measures = "ri",
    injection_mass = lals_fixture$injection_mass_mg,
    flow_rate = lals_fixture$flow_ml_min,
    concentration_units = case$concentration_units,
    elution_units = case$elution_units,
    id = stable_step_id("qual_lals_mass", case$id)
  )
  recipe <- measure.sec::step_sec_detector_delay(
    recipe,
    reference = "lals",
    targets = "ri",
    delay_times = stats::setNames(lals_fixture$delay_time_min, "ri"),
    flow_rate = lals_fixture$flow_ml_min,
    elution_units = case$elution_units,
    id = stable_step_id("qual_lals_delay", case$id)
  )
  recipe <- measure.sec::step_sec_lals(
    recipe,
    measures = "lals",
    concentration_col = "ri",
    optical_constant = lals_fixture$optical_constant,
    calibration_constant = lals_fixture$calibration_constant,
    min_signal = 0,
    id = stable_step_id("qual_lals", case$id)
  )
  prepared <- run_with_warnings(
    recipes::prep(recipe, training = data),
    paste0(case$id, " prep")
  )
  baked <- run_with_warnings(
    recipes::bake(prepared$value, new_data = data),
    paste0(case$id, " bake")
  )$value
  tag <- attr(baked$ri, "sec_concentration_unit", exact = TRUE)
  if (!identical(tag, case$concentration_units)) {
    stop(
      "Downstream LALS input lost the current concentration unit tag.",
      call. = FALSE
    )
  }
  list(
    id = json_scalar(case$id),
    concentration = curve_json(baked$ri[[1L]]),
    scattering = curve_json(baked$lals[[1L]]),
    mw = curve_json(baked$mw_lals[[1L]]),
    concentration_unit = json_scalar(tag)
  )
}

simple_mass_recipe <- function(
  data,
  units = "min",
  flow = 0.5,
  injection_mass = 2,
  id = "invalid_mass"
) {
  recipe <- recipe_for(data)
  measure.sec::step_sec_concentration(
    recipe,
    measures = "ri",
    injection_mass = injection_mass,
    flow_rate = flow,
    elution_units = units,
    id = stable_step_id("qual_invalid_mass", id)
  )
}

simple_delay_case <- function(
  location = c(0, 1, 2),
  signal = c(0, 1, 0),
  reference_location = location,
  reference_signal = c(10, 11, 12),
  method = "interpolate",
  units = "min",
  flow = 0.5,
  delay_times = c(target = 0.5),
  delay_volumes = NULL,
  targets = "target",
  reference = "reference",
  id = "invalid_delay"
) {
  target <- measure::new_measure_tbl(
    location = location,
    value = signal
  )
  reference_curve <- measure::new_measure_tbl(
    location = reference_location,
    value = reference_signal
  )
  data <- data.frame(row_id = 1L)
  data$reference <- measure::new_measure_list(list(reference_curve))
  data$target <- measure::new_measure_list(list(target))
  attr(data$target, "sec_concentration_unit") <- "mg/mL"
  recipe <- recipe_for(data)
  recipe <- measure.sec::step_sec_detector_delay(
    recipe,
    reference = reference,
    targets = targets,
    delay_times = delay_times,
    delay_volumes = delay_volumes,
    flow_rate = flow,
    method = method,
    elution_units = units,
    id = stable_step_id("qual_invalid_delay", id)
  )
  list(recipe = recipe, data = data)
}

invalid_data <- function(location = c(0, 1, 2), value = c(0, 1, 0)) {
  data <- data.frame(row_id = 1L)
  data$ri <- measure::new_measure_list(list(
    measure::new_measure_tbl(location = location, value = value)
  ))
  data
}

expect_prep_bake_error <- function(recipe, training = NULL, new_data = NULL) {
  if (is.list(recipe) && all(c("recipe", "data") %in% names(recipe))) {
    training <- recipe$data
    new_data <- training
    recipe <- recipe$recipe
  }
  if (is.null(new_data)) {
    new_data <- training
  }
  prepared <- run_with_warnings(
    recipes::prep(recipe, training = training),
    "invalid-case prep"
  )$value
  run_with_warnings(
    recipes::bake(prepared, new_data = new_data),
    "invalid-case bake"
  )
  invisible(NULL)
}

run_invalid_case <- function(id) {
  switch(
    id,
    "concentration-invalid-unit" = {
      data <- invalid_data()
      recipe <- recipe_for(data)
      recipe <- measure.sec::step_sec_concentration(
        recipe,
        measures = "ri",
        injection_mass = 2,
        flow_rate = 0.5,
        elution_units = "minute",
        id = "reject_invalid_concentration_unit"
      )
      expect_prep_bake_error(recipe, data)
    },
    "concentration-null-time-flow" = {
      data <- invalid_data()
      recipe <- simple_mass_recipe(data, flow = NULL, id = id)
      expect_prep_bake_error(recipe, data)
    },
    "concentration-nonpositive-flow" = {
      data <- invalid_data()
      recipe <- simple_mass_recipe(data, flow = 0, id = id)
      expect_prep_bake_error(recipe, data)
    },
    "concentration-duplicate-axis" = {
      data <- invalid_data(location = c(0, 1, 1), value = c(0, 1, 0))
      recipe <- simple_mass_recipe(data, id = id)
      expect_prep_bake_error(recipe, data)
    },
    "concentration-nonfinite-axis" = {
      data <- invalid_data(location = c(0, Inf, 2), value = c(0, 1, 0))
      recipe <- simple_mass_recipe(data, id = id)
      expect_prep_bake_error(recipe, data)
    },
    "concentration-matrix-signal" = {
      data <- invalid_data(
        location = c(0, 1, 2),
        value = matrix(c(0, 1, 0), ncol = 1L)
      )
      recipe <- simple_mass_recipe(data, id = id)
      expect_prep_bake_error(recipe, data)
    },
    "concentration-no-positive-area" = {
      data <- invalid_data(location = c(0, 1, 2), value = c(0, 0, 0))
      recipe <- simple_mass_recipe(data, id = id)
      expect_prep_bake_error(recipe, data)
    },
    "concentration-unit-explicit-null-mutated" = {
      data <- invalid_data()
      recipe <- simple_mass_recipe(data, id = id)
      recipe <- set_step_field(recipe, "elution_units", NULL, keep_null = TRUE)
      expect_prep_bake_error(recipe, data)
    },
    "delay-invalid-unit" = {
      data <- delay_data(
        list(
          time_min = list(0, 1, 2),
          signal = list(0, 1, 0),
          reference = list(10, 11, 12)
        ),
        axis_per_minute = 1
      )
      recipe <- delay_recipe(
        data,
        list(
          id = id,
          method = "interpolate",
          flow_ml_min = 0.5,
          elution_units = "minute"
        ),
        0.5,
        "time"
      )
      expect_prep_bake_error(recipe, data)
    },
    "delay-null-cross-dimension-flow" = {
      attempt <- simple_delay_case(
        location = c(0, 1, 2),
        signal = c(0, 1, 0),
        flow = NULL,
        delay_times = NULL,
        delay_volumes = c(target = 0.5),
        id = id
      )
      expect_prep_bake_error(attempt$recipe, attempt$data)
    },
    "delay-duplicate-names" = {
      attempt <- simple_delay_case(
        location = c(0, 1, 2),
        signal = c(0, 1, 0),
        delay_times = c(target = 0.5, target = 0.6),
        id = id
      )
      expect_prep_bake_error(attempt)
    },
    "delay-blank-name" = {
      attempt <- simple_delay_case(
        location = c(0, 1, 2),
        signal = c(0, 1, 0),
        delay_times = stats::setNames(0.5, ""),
        id = id
      )
      expect_prep_bake_error(attempt)
    },
    "delay-reference-is-target" = {
      attempt <- simple_delay_case(
        location = c(0, 1, 2),
        signal = c(0, 1, 0),
        targets = c("reference"),
        delay_times = c(reference = 0.5),
        id = id
      )
      expect_prep_bake_error(attempt)
    },
    "delay-missing-current-reference" = {
      attempt <- simple_delay_case(id = id)
      prepared <- run_with_warnings(
        recipes::prep(attempt$recipe, training = attempt$data),
        "invalid-case prep"
      )$value
      current <- data.frame(row_id = 1L)
      current$target <- attempt$data$target
      run_with_warnings(
        recipes::bake(prepared, new_data = current),
        "invalid-case missing-reference bake"
      )
    },
    "delay-missing-current-target" = {
      attempt <- simple_delay_case(id = id)
      prepared <- run_with_warnings(
        recipes::prep(attempt$recipe, training = attempt$data),
        "invalid-case prep"
      )$value
      current <- data.frame(row_id = 1L)
      current$reference <- attempt$data$reference
      run_with_warnings(
        recipes::bake(prepared, new_data = current),
        "invalid-case missing-target bake"
      )
    },
    "delay-mismatched-current-grid" = {
      attempt <- simple_delay_case(id = id)
      prepared <- run_with_warnings(
        recipes::prep(attempt$recipe, training = attempt$data),
        "invalid-case prep"
      )$value
      current <- delay_data(
        list(
          time_min = list(0, 1, 2),
          signal = list(0, 1, 0),
          reference = list(10, 11, 12)
        ),
        axis_per_minute = 1,
        reference_location = c(0, 1, 2),
        target_location = c(0, 1.01, 2)
      )
      run_with_warnings(
        recipes::bake(prepared, new_data = current),
        "invalid-case mismatched-grid bake"
      )
    },
    "delay-nonzero-irregular-index" = {
      attempt <- simple_delay_case(
        location = c(0, 1, 3, 4),
        signal = c(0, 1, 1, 0),
        reference_location = c(0, 1, 3, 4),
        reference_signal = c(10, 11, 13, 14),
        method = "shift",
        delay_times = c(target = 1),
        id = id
      )
      expect_prep_bake_error(attempt)
    },
    "delay-short-nonzero" = {
      attempt <- simple_delay_case(
        location = 0,
        signal = 1,
        reference_location = 0,
        reference_signal = 10,
        id = id
      )
      expect_prep_bake_error(attempt)
    },
    "delay-nonfinite-axis" = {
      attempt <- simple_delay_case(
        location = c(0, Inf, 2),
        signal = c(0, 1, 0),
        reference_location = c(0, Inf, 2),
        reference_signal = c(10, 11, 12),
        id = id
      )
      expect_prep_bake_error(attempt)
    },
    "delay-nan-signal" = {
      attempt <- simple_delay_case(
        location = c(0, 1, 2),
        signal = c(0, NaN, 0),
        id = id
      )
      expect_prep_bake_error(attempt)
    },
    "delay-infinite-signal" = {
      attempt <- simple_delay_case(
        location = c(0, 1, 2),
        signal = c(0, Inf, 0),
        id = id
      )
      expect_prep_bake_error(attempt)
    },
    "delay-matrix-signal" = {
      attempt <- simple_delay_case(
        location = c(0, 1, 2),
        signal = matrix(c(0, 1, 0), ncol = 1L),
        id = id
      )
      expect_prep_bake_error(attempt)
    },
    "delay-conversion-overflow" = {
      attempt <- simple_delay_case(
        location = c(0, 1, 2),
        signal = c(0, 1, 0),
        units = "mL",
        flow = 1e308,
        delay_times = c(target = 1e308),
        id = id
      )
      expect_prep_bake_error(attempt)
    },
    "delay-conversion-underflow" = {
      attempt <- simple_delay_case(
        location = c(0, 1, 2),
        signal = c(0, 1, 0),
        units = "mL",
        flow = 1e-300,
        delay_times = c(target = 1e-300),
        id = id
      )
      expect_prep_bake_error(attempt)
    },
    "delay-translation-cancellation" = {
      location <- c(1e16, 1e16 + 4, 1e16 + 8)
      attempt <- simple_delay_case(
        location = location,
        signal = c(0, 1, 0),
        reference_location = location,
        reference_signal = c(10, 11, 12),
        delay_times = c(target = 0.25),
        id = id
      )
      expect_prep_bake_error(attempt)
    },
    "delay-unit-explicit-null-mutated" = {
      attempt <- simple_delay_case(id = id)
      attempt$recipe <- set_step_field(
        attempt$recipe,
        "elution_units",
        NULL,
        keep_null = TRUE
      )
      expect_prep_bake_error(attempt)
    },
    stop("Unknown frozen invalid case id: ", id, call. = FALSE)
  )
}

reject_case <- function(id) {
  error <- tryCatch(
    {
      run_invalid_case(id)
      NULL
    },
    unexpected_qualification_warning = function(condition) {
      stop(condition)
    },
    error = function(condition) condition
  )
  if (is.null(error)) {
    stop("Expected rejection case did not error: ", id, call. = FALSE)
  }
  message <- conditionMessage(error)
  if (length(message) != 1L || is.na(message) || !nzchar(message)) {
    stop("Expected rejection case had an empty diagnostic: ", id, call. = FALSE)
  }
  list(
    id = json_scalar(id),
    error = json_scalar(TRUE),
    message = json_scalar(message)
  )
}

run_qualification <- function() {
  trailing <- commandArgs(trailingOnly = TRUE)
  if (length(trailing) != 1L || !nzchar(trailing[[1L]])) {
    stop(
      "Usage: Rscript tools/sec-elution/qualify.R NEW_OUTPUT_ROOT",
      call. = FALSE
    )
  }
  script_arg <- Filter(
    function(argument) startsWith(argument, "--file="),
    commandArgs(trailingOnly = FALSE)
  )
  if (length(script_arg) != 1L) {
    stop("Could not determine the qualification script path.", call. = FALSE)
  }
  script_path <- normalizePath(
    decode_rscript_path(sub("^--file=", "", script_arg[[1L]])),
    winslash = "/",
    mustWork = TRUE
  )
  tool_directory <- dirname(script_path)
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop(
      "Preflight failed: jsonlite is required to record run receipts; ",
      "the output root was not reserved.",
      call. = FALSE
    )
  }
  package_root <- normalizePath(
    dirname(dirname(tool_directory)),
    winslash = "/",
    mustWork = TRUE
  )
  output <- normalize_input_path(trailing[[1L]])
  link_target <- Sys.readlink(output)
  has_link <- length(link_target) == 1L &&
    !is.na(link_target) &&
    nzchar(link_target)
  if (file.exists(output) || dir.exists(output) || has_link) {
    stop(
      "Output root already exists; choose a fresh, non-existing directory: ",
      output,
      call. = FALSE
    )
  }
  if (!dir.exists(dirname(output))) {
    stop(
      "The parent directory for the output root must already exist: ",
      dirname(output),
      call. = FALSE
    )
  }
  if (!dir.create(output, recursive = FALSE, showWarnings = FALSE)) {
    stop("Could not reserve fresh output root: ", output, call. = FALSE)
  }

  run <- new.env(parent = emptyenv())
  run$output <- output
  run$stage <- "initialize"
  run$current_case <- NULL
  run$error <- NULL
  run$counts <- list(
    mass = 0L,
    delay = 0L,
    noop = 0L,
    lals = 0L,
    invalid_cases = 0L
  )
  run$completed <- FALSE
  status_path <- file.path(output, "run-status.json")
  progress_path <- file.path(output, "qualification-progress.json")

  write_status <- function() {
    write_json(
      list(
        complete = json_scalar(run$completed),
        stage = json_scalar(run$stage),
        error = if (is.null(run$error)) NULL else json_scalar(run$error),
        updated_utc = json_scalar(format(Sys.time(), tz = "UTC", usetz = TRUE))
      ),
      status_path
    )
  }
  write_progress <- function() {
    write_json(
      list(
        complete = json_scalar(run$completed),
        stage = json_scalar(run$stage),
        current_case = if (is.null(run$current_case)) {
          NULL
        } else {
          json_scalar(run$current_case)
        },
        completed_counts = lapply(run$counts, json_scalar),
        error = if (is.null(run$error)) NULL else json_scalar(run$error),
        updated_utc = json_scalar(format(Sys.time(), tz = "UTC", usetz = TRUE))
      ),
      progress_path
    )
  }
  advance <- function(stage, current_case = NULL) {
    run$stage <- stage
    run$current_case <- current_case
    write_progress()
    write_status()
  }

  write_status()
  write_progress()
  tryCatch(
    {
      advance("load-inputs")
      for (dependency in c("digest", "recipes")) {
        if (!requireNamespace(dependency, quietly = TRUE)) {
          stop(
            "Qualification prerequisite is not installed: ",
            dependency,
            call. = FALSE
          )
        }
      }

      input_paths <- c(
        fixtures = file.path(tool_directory, "fixtures.json"),
        declaration = file.path(tool_directory, "declaration.json")
      )
      if (!all(file.exists(input_paths))) {
        stop("Frozen SEC qualification input file is missing.", call. = FALSE)
      }
      input_hashes <- vapply(input_paths, sha256_file, character(1))
      if (!identical(unname(input_hashes), unname(qualification_sha256))) {
        stop("Frozen SEC qualification input SHA-256 mismatch.", call. = FALSE)
      }
      fixtures <- read_json(input_paths[["fixtures"]])
      declaration <- read_json(input_paths[["declaration"]])
      if (
        !identical(fixtures$profile, declaration$profile) ||
          !identical(fixtures$profile, "sec-elution-operations/1")
      ) {
        stop("Frozen SEC qualification profile mismatch.", call. = FALSE)
      }

      advance("verify-installed-packages")
      installed <- list(
        measure = package_provenance("measure"),
        measure.sec = package_provenance("measure.sec")
      )
      dependency_versions <- list(
        digest = as.character(utils::packageVersion("digest")),
        jsonlite = as.character(utils::packageVersion("jsonlite")),
        recipes = as.character(utils::packageVersion("recipes")),
        measure = as.character(utils::packageVersion("measure")),
        measure.sec = as.character(utils::packageVersion("measure.sec"))
      )
      checkout_head <- tool_checkout_head(package_root)
      fixture_map <- setNames(
        fixtures$mass_traces,
        vapply(
          fixtures$mass_traces,
          function(item) item$id,
          character(1)
        )
      )
      delay_fixture_map <- setNames(
        fixtures$delay_traces,
        vapply(
          fixtures$delay_traces,
          function(item) item$id,
          character(1)
        )
      )
      results <- list(
        profile = json_scalar(declaration$profile),
        input_files_sha256 = list(
          `fixtures.json` = json_scalar(input_hashes[["fixtures"]]),
          `declaration.json` = json_scalar(input_hashes[["declaration"]])
        ),
        installed_packages = installed,
        R = json_scalar(R.version.string),
        tool_checkout_head = if (is.null(checkout_head)) {
          NULL
        } else {
          json_scalar(checkout_head)
        },
        mass = list(),
        delay = list(),
        noop = list(),
        lals = list(),
        invalid_cases = list()
      )

      advance("mass")
      results$mass <- lapply(declaration$mass, function(case) {
        advance("mass", case$id)
        result <- run_mass_case(case, fixture_map)
        run$counts$mass <- run$counts$mass + 1L
        write_progress()
        result
      })

      advance("delay")
      results$delay <- lapply(declaration$delay, function(case) {
        advance("delay", case$id)
        result <- run_delay_case(case, delay_fixture_map)
        run$counts$delay <- run$counts$delay + 1L
        write_progress()
        result
      })

      advance("noop")
      results$noop <- lapply(declaration$noop, function(case) {
        advance("noop", case$id)
        result <- run_noop_case(case, fixtures$noop_trace)
        run$counts$noop <- run$counts$noop + 1L
        write_progress()
        result
      })

      advance("lals")
      results$lals <- lapply(declaration$lals, function(case) {
        advance("lals", case$id)
        result <- run_lals_case(case, fixtures$lals)
        run$counts$lals <- run$counts$lals + 1L
        write_progress()
        result
      })

      advance("invalid-cases")
      results$invalid_cases <- lapply(declaration$invalid_cases, function(id) {
        advance("invalid-cases", id)
        result <- reject_case(id)
        run$counts$invalid_cases <- run$counts$invalid_cases + 1L
        write_progress()
        result
      })

      write_new_json(results, file.path(output, "results.json"))
      runtime <- list(
        created_utc = json_scalar(format(Sys.time(), tz = "UTC", usetz = TRUE)),
        session_info = lapply(
          as.list(utils::capture.output(utils::sessionInfo())),
          json_scalar
        ),
        dependency_versions = lapply(dependency_versions, json_scalar)
      )
      write_new_json(runtime, file.path(output, "runtime.json"))
      run$stage <- "complete"
      run$current_case <- NULL
      run$completed <- TRUE
      run$error <- NULL
      write_progress()
      write_status()
      invisible(output)
    },
    error = function(condition) {
      run$error <- conditionMessage(condition)
      if (is.null(run$error) || !nzchar(run$error)) {
        run$error <- paste0(
          "Qualification failed with ",
          paste(class(condition), collapse = "/")
        )
      }
      try(write_progress(), silent = TRUE)
      try(write_status(), silent = TRUE)
      stop(condition)
    }
  )
}

withCallingHandlers(
  run_qualification(),
  warning = function(warning) {
    message <- paste0(
      "Unexpected warning outside an allowed case handler: ",
      conditionMessage(warning)
    )
    stop(structure(
      list(message = message, call = NULL, warning = warning),
      class = c("unexpected_qualification_warning", "error", "condition")
    ))
  }
)
