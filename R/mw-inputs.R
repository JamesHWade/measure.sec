# Helpers for pairing a concentration signal with a per-slice MW curve.
# The tag describes how a curve was generated; the locations are still checked
# for every sample at bake time.
.sec_tag_mw_curve <- function(curves, source, scale) {
  attr(curves, "sec_mw_curve") <- list(source = source, scale = scale)
  curves
}

.sec_resolve_mw_input <- function(
  training,
  measures,
  mw_col,
  mw_scale,
  calibration
) {
  available <- find_measure_cols(training)
  tags <- lapply(
    training[available],
    attr,
    which = "sec_mw_curve",
    exact = TRUE
  )

  if (is.null(measures)) {
    measures <- available[vapply(tags, is.null, logical(1))]
  }
  if (length(measures) != 1 || !measures %in% available) {
    cli::cli_abort(
      "Select exactly one concentration measure column with {.arg measures}."
    )
  }
  if (!is.null(mw_col) && !is.null(calibration)) {
    cli::cli_abort(
      "Specify either {.arg mw_col} or {.arg calibration}, not both."
    )
  }

  # A calibration step records both its input signal and output scale. Only an
  # unambiguous matching curve is selected automatically.
  if (is.null(mw_col) && is.null(calibration)) {
    matches <- names(tags)[vapply(
      tags,
      function(tag) {
        is.list(tag) && identical(tag$source, measures)
      },
      logical(1)
    )]
    if (length(matches) > 1) {
      cli::cli_abort(
        "Multiple molecular-weight curves match {.val {measures}}; specify {.arg mw_col}."
      )
    }
    if (length(matches) == 1) {
      mw_col <- matches[[1]]
    } else {
      uncalibrated <- vapply(
        training[available],
        function(curve) {
          identical(
            attr(curve, "sec_scattering_basis", exact = TRUE),
            "relative_uncalibrated"
          ) &&
            identical(
              attr(curve, "sec_scattering_source", exact = TRUE),
              measures
            )
        },
        logical(1)
      )
      if (any(uncalibrated)) {
        cli::cli_abort(
          "The matching scattering curve is uncalibrated and cannot provide molecular weights."
        )
      }
    }
  }

  if (!is.null(mw_col)) {
    if (
      !is.character(mw_col) || length(mw_col) != 1 || !mw_col %in% available
    ) {
      cli::cli_abort(
        "{.arg mw_col} must name one molecular-weight measure column."
      )
    }
    if (identical(mw_col, measures)) {
      cli::cli_abort(
        "{.arg mw_col} and the concentration measure must be different columns."
      )
    }
    if (
      identical(
        attr(training[[mw_col]], "sec_scattering_basis", exact = TRUE),
        "relative_uncalibrated"
      )
    ) {
      cli::cli_abort(
        "An uncalibrated scattering curve cannot provide molecular weights."
      )
    }
    tag <- attr(training[[mw_col]], "sec_mw_curve", exact = TRUE)
    if (is.null(mw_scale) && is.list(tag)) {
      mw_scale <- tag$scale
    }
    if (
      !is.character(mw_scale) ||
        length(mw_scale) != 1 ||
        is.na(mw_scale) ||
        !mw_scale %in% c("log10", "absolute")
    ) {
      cli::cli_abort(
        "Specify {.arg mw_scale} as {.val log10} or {.val absolute} for {.arg mw_col}."
      )
    }
    if (is.list(tag) && !identical(tag$scale, mw_scale)) {
      cli::cli_abort(
        "{.arg mw_scale} conflicts with the recorded scale of {.val {mw_col}}."
      )
    }
  } else if (!is.null(mw_scale)) {
    cli::cli_abort(
      "{.arg mw_scale} requires a molecular-weight curve in {.arg mw_col}."
    )
  }

  list(measures = measures, mw_col = mw_col, mw_scale = mw_scale)
}

.sec_validate_integration_range <- function(integration_range) {
  if (is.null(integration_range)) {
    return(NULL)
  }

  if (
    !is.numeric(integration_range) ||
      length(integration_range) != 2 ||
      !all(is.finite(integration_range)) ||
      integration_range[1] >= integration_range[2]
  ) {
    cli::cli_abort(
      "{.arg integration_range} must contain two increasing finite values."
    )
  }

  as.numeric(integration_range)
}

.sec_empty_excluded_regions <- function() {
  tibble::tibble(
    start = numeric(),
    end = numeric(),
    purpose = character()
  )
}

.sec_validate_excluded_regions <- function(excluded_regions) {
  if (is.null(excluded_regions)) {
    return(.sec_empty_excluded_regions())
  }

  if (!is.data.frame(excluded_regions)) {
    cli::cli_abort("{.arg .excluded_regions} must be a data frame or NULL.")
  }
  required <- c("start", "end", "purpose")
  if (!all(required %in% names(excluded_regions))) {
    cli::cli_abort(
      "{.arg .excluded_regions} must contain numeric {.field start} and {.field end} columns and a {.field purpose} column."
    )
  }
  if (
    !is.numeric(excluded_regions$start) ||
      !is.numeric(excluded_regions$end) ||
      !is.character(excluded_regions$purpose)
  ) {
    cli::cli_abort(
      "{.arg .excluded_regions} has invalid column types for {.field start}, {.field end}, or {.field purpose}."
    )
  }
  if (
    !all(is.finite(excluded_regions$start)) ||
      !all(is.finite(excluded_regions$end)) ||
      any(excluded_regions$start >= excluded_regions$end) ||
      anyNA(excluded_regions$purpose) ||
      !all(excluded_regions$purpose %in% c("baseline", "integration", "both"))
  ) {
    cli::cli_abort(
      "{.arg .excluded_regions} contains malformed region metadata."
    )
  }

  excluded_regions <- excluded_regions[
    excluded_regions$purpose %in% c("integration", "both"),
    required,
    drop = FALSE
  ]
  tibble::as_tibble(excluded_regions)
}

.sec_row_integration_metadata <- function(
  new_data,
  i,
  integration_range = NULL
) {
  if (!is.data.frame(new_data)) {
    cli::cli_abort("Integration metadata requires tabular {.arg new_data}.")
  }
  if (
    length(i) != 1 ||
      is.na(i) ||
      !is.numeric(i) ||
      i < 1 ||
      i > nrow(new_data) ||
      i != as.integer(i)
  ) {
    cli::cli_abort("The row index {.arg i} is outside {.arg new_data}.")
  }
  i <- as.integer(i)

  resolved_range <- .sec_validate_integration_range(integration_range)
  if (is.null(resolved_range) && ".integration_window" %in% names(new_data)) {
    windows <- new_data[[".integration_window"]]
    if (!is.list(windows) || length(windows) != nrow(new_data)) {
      cli::cli_abort(
        "Column {.field .integration_window} must be a list column with one entry per row."
      )
    }
    window <- windows[[i]]
    if (!is.null(window)) {
      if (
        !is.data.frame(window) ||
          nrow(window) != 1 ||
          !all(c("start", "end") %in% names(window)) ||
          !is.numeric(window$start) ||
          !is.numeric(window$end)
      ) {
        cli::cli_abort(
          "Row {i} of {.field .integration_window} must be NULL or contain one numeric {.field start}/{.field end} row."
        )
      }
      resolved_range <- .sec_validate_integration_range(
        c(window$start[[1]], window$end[[1]])
      )
    }
  }

  excluded_regions <- .sec_empty_excluded_regions()
  if (".excluded_regions" %in% names(new_data)) {
    exclusions <- new_data[[".excluded_regions"]]
    if (!is.list(exclusions) || length(exclusions) != nrow(new_data)) {
      cli::cli_abort(
        "Column {.field .excluded_regions} must be a list column with one entry per row."
      )
    }
    excluded_regions <- .sec_validate_excluded_regions(exclusions[[i]])
  }

  list(
    integration_range = resolved_range,
    excluded_regions = excluded_regions
  )
}

.sec_validate_mw_axis <- function(location, label) {
  if (!is.numeric(location) || !all(is.finite(location))) {
    cli::cli_abort("{label} locations must be finite numeric values.")
  }
  if (
    length(location) > 1 &&
      !(all(diff(location) > 0) || all(diff(location) < 0))
  ) {
    cli::cli_abort("{label} locations must be strictly monotonic and unique.")
  }
  invisible(location)
}

.sec_interpolate_linear <- function(location, value, at) {
  if (length(at) == 0) {
    return(numeric())
  }
  if (length(location) == 0) {
    return(rep(NA_real_, length(at)))
  }

  order <- order(location)
  x <- location[order]
  y <- value[order]
  exact <- match(at, x)
  result <- rep(NA_real_, length(at))
  on_axis <- !is.na(exact)
  result[on_axis] <- y[exact[on_axis]]

  between <- which(
    !on_axis & at >= x[1] & at <= x[length(x)]
  )
  if (length(between) == 0 || length(x) < 2) {
    return(result)
  }

  left <- findInterval(at[between], x)
  left <- pmax(1L, pmin(left, length(x) - 1L))
  right <- left + 1L
  valid <- is.finite(y[left]) & is.finite(y[right])
  if (any(valid)) {
    indices <- between[valid]
    x0 <- x[left[valid]]
    x1 <- x[right[valid]]
    y0 <- y[left[valid]]
    y1 <- y[right[valid]]
    result[indices] <- y0 + (at[indices] - x0) * (y1 - y0) / (x1 - x0)
  }

  result
}

.sec_merge_intervals <- function(regions, low, high) {
  if (nrow(regions) == 0 || low == high) {
    return(matrix(numeric(), ncol = 2))
  }
  overlaps <- regions$end > low & regions$start < high
  regions <- regions[overlaps, , drop = FALSE]
  if (nrow(regions) == 0) {
    return(matrix(numeric(), ncol = 2))
  }

  regions$start <- pmax(regions$start, low)
  regions$end <- pmin(regions$end, high)
  regions <- regions[order(regions$start, regions$end), , drop = FALSE]
  merged <- list()
  current <- c(regions$start[1], regions$end[1])
  if (nrow(regions) > 1) {
    for (j in 2:nrow(regions)) {
      if (regions$start[j] <= current[2]) {
        current[2] <- max(current[2], regions$end[j])
      } else {
        merged[[length(merged) + 1]] <- current
        current <- c(regions$start[j], regions$end[j])
      }
    }
  }
  merged[[length(merged) + 1]] <- current
  do.call(rbind, merged)
}

.sec_support_intervals <- function(low, high, excluded_regions) {
  if (low == high) {
    excluded <- nrow(excluded_regions) > 0 &&
      any(
        excluded_regions$start <= low & excluded_regions$end >= high
      )
    if (excluded) {
      return(matrix(numeric(), ncol = 2))
    }
    return(matrix(c(low, high), nrow = 1))
  }

  excluded <- .sec_merge_intervals(excluded_regions, low, high)
  intervals <- list()
  cursor <- low
  if (nrow(excluded) > 0) {
    for (j in seq_len(nrow(excluded))) {
      if (excluded[j, 1] > cursor) {
        intervals[[length(intervals) + 1]] <- c(cursor, excluded[j, 1])
      }
      cursor <- max(cursor, excluded[j, 2])
    }
  }
  if (cursor < high) {
    intervals[[length(intervals) + 1]] <- c(cursor, high)
  }
  if (length(intervals) == 0) {
    return(matrix(numeric(), ncol = 2))
  }
  do.call(rbind, intervals)
}

.sec_segment_widths <- function(location) {
  if (length(location) < 2) {
    return(rep(0, length(location)))
  }
  step <- abs(diff(location))
  c(
    step[1] / 2,
    if (length(step) > 1) (step[-1] + step[-length(step)]) / 2,
    step[length(step)] / 2
  )
}

.sec_mw_slices <- function(
  signal,
  mw_curve = NULL,
  mw_scale = NULL,
  calibration = NULL,
  integration_range = NULL,
  excluded_regions = NULL
) {
  if (!is.list(signal) || is.null(signal$location) || is.null(signal$value)) {
    cli::cli_abort("The concentration curve must contain locations and values.")
  }
  location <- signal$location
  value <- signal$value
  .sec_validate_mw_axis(location, "Concentration curve")
  if (!is.numeric(value) || length(location) != length(value)) {
    cli::cli_abort(
      "The concentration curve must have finite locations and matching values."
    )
  }

  if (!is.null(mw_curve)) {
    if (
      !is.list(mw_curve) ||
        is.null(mw_curve$location) ||
        is.null(mw_curve$value)
    ) {
      cli::cli_abort(
        "The molecular-weight curve must contain locations and values."
      )
    }
    .sec_validate_mw_axis(mw_curve$location, "Molecular-weight curve")
    if (
      !isTRUE(all.equal(
        location,
        mw_curve$location,
        tolerance = 1e-9,
        check.attributes = FALSE
      ))
    ) {
      cli::cli_abort(
        "Molecular-weight and concentration curve locations must align."
      )
    }
    mw_values <- mw_curve$value
    if (!is.numeric(mw_values) || length(mw_values) != length(location)) {
      cli::cli_abort(
        "Molecular-weight and concentration curves must have the same numeric values length."
      )
    }
    if (
      !is.null(mw_scale) &&
        (!is.character(mw_scale) ||
          length(mw_scale) != 1 ||
          is.na(mw_scale) ||
          !mw_scale %in% c("log10", "absolute"))
    ) {
      cli::cli_abort(
        "{.arg mw_scale} must be {.val log10} or {.val absolute}."
      )
    }
    mw_scale <- mw_scale %||% "log10"
  } else {
    mw_values <- NULL
  }

  integration_range <- .sec_validate_integration_range(integration_range)
  excluded_regions <- .sec_validate_excluded_regions(excluded_regions)
  if (!all(is.finite(location))) {
    cli::cli_abort(
      "The concentration curve must have finite locations and matching values."
    )
  }
  if (length(location) == 0) {
    return(list(
      location = numeric(),
      signal = numeric(),
      mw = numeric(),
      weight = numeric(),
      segment = integer()
    ))
  }

  hull <- range(location)
  low <- max(hull[1], integration_range[1] %||% hull[1])
  high <- min(hull[2], integration_range[2] %||% hull[2])
  if (low > high) {
    return(list(
      location = numeric(),
      signal = numeric(),
      mw = numeric(),
      weight = numeric(),
      segment = integer()
    ))
  }

  intervals <- .sec_support_intervals(low, high, excluded_regions)
  if (nrow(intervals) == 0) {
    return(list(
      location = numeric(),
      signal = numeric(),
      mw = numeric(),
      weight = numeric(),
      segment = integer()
    ))
  }

  descending <- length(location) > 1 && location[1] > location[length(location)]
  if (descending) {
    intervals <- intervals[rev(seq_len(nrow(intervals))), , drop = FALSE]
  }
  locations <- vector("list", nrow(intervals))
  segment_ids <- vector("list", nrow(intervals))
  for (j in seq_len(nrow(intervals))) {
    bounds <- intervals[j, ]
    points <- sort(unique(c(
      bounds[1],
      location[location > bounds[1] & location < bounds[2]],
      bounds[2]
    )))
    if (descending) {
      points <- rev(points)
    }
    locations[[j]] <- points
    segment_ids[[j]] <- rep.int(j, length(points))
  }
  selected_location <- unlist(locations, use.names = FALSE)
  segment <- as.integer(unlist(segment_ids, use.names = FALSE))
  selected_value <- .sec_interpolate_linear(location, value, selected_location)
  if (!all(is.finite(selected_value))) {
    cli::cli_abort(
      "The selected concentration signal contains missing or non-finite values."
    )
  }

  if (!is.null(mw_values)) {
    mw_as_absolute <- if (identical(mw_scale, "absolute")) {
      mw_values
    } else {
      10^mw_values
    }
    invalid_source_mw <- !is.finite(mw_as_absolute) | mw_as_absolute <= 0
    removed_from_support <- location < low | location > high
    if (nrow(excluded_regions) > 0) {
      removed_from_support <- removed_from_support |
        vapply(
          location,
          function(x) {
            any(excluded_regions$start <= x & excluded_regions$end >= x)
          },
          logical(1)
        )
    }
    ignore_invalid <- invalid_source_mw & removed_from_support
    mw_interpolation_location <- location[!ignore_invalid]
    mw_interpolation_values <- mw_values[!ignore_invalid]
    selected_mw_values <- .sec_interpolate_linear(
      mw_interpolation_location,
      mw_interpolation_values,
      selected_location
    )
    mw <- if (identical(mw_scale, "absolute")) {
      selected_mw_values
    } else {
      10^selected_mw_values
    }
  } else if (is.null(calibration)) {
    # Backward-compatible mode: the signal axis is log10(MW).
    mw <- 10^selected_location
  } else if (
    is.numeric(calibration) &&
      length(calibration) == 2 &&
      all(is.finite(calibration))
  ) {
    mw <- 10^(calibration[1] * selected_location + calibration[2])
  } else if (identical(calibration, "auto")) {
    if (length(selected_location) >= 2) {
      cli::cli_warn(
        "The {.val auto} calibration is synthetic and must not be treated as an experimental calibration."
      )
    }
    synthetic_log_mw <- seq(7, 2, length.out = length(location))
    selected_log_mw <- .sec_interpolate_linear(
      location,
      synthetic_log_mw,
      selected_location
    )
    mw <- 10^selected_log_mw
  } else {
    cli::cli_abort(
      "{.arg calibration} must be NULL, a finite slope/intercept pair, or {.val auto}."
    )
  }

  signal_value <- pmax(selected_value, 0)
  widths <- unlist(lapply(locations, .sec_segment_widths), use.names = FALSE)
  weight <- signal_value * widths
  bad_mw <- !is.finite(mw) | mw <= 0
  if (any(bad_mw & weight > 0)) {
    cli::cli_abort(
      "Positive signal includes slices without a valid molecular-weight calibration; adjust the integration range or calibration."
    )
  }
  if (any(selected_value < 0)) {
    cli::cli_warn(
      "Negative concentration values were clipped to zero for molecular-weight calculations."
    )
  }

  list(
    location = selected_location,
    signal = signal_value,
    mw = mw,
    weight = weight,
    segment = segment
  )
}
