# ==============================================================================
# step_sec_dls
#
# Dynamic light scattering processing
# ==============================================================================

#' Dynamic Light Scattering Processing for SEC
#'
#' `step_sec_dls()` creates a *specification* of a recipe step that processes
#' dynamic light scattering (DLS) correlation data to estimate diffusion
#' coefficient and hydrodynamic radius (Rh).
#'
#' @param recipe A recipe object.
#' @param measures Character vector of DLS measure columns. If `NULL`, the
#'   step searches for measure columns containing "dls" or "corr".
#' @param temperature Temperature in degrees Celsius.
#' @param viscosity Solvent viscosity in mPa*s. If `NULL`, uses a water
#'   approximation based on temperature.
#' @param laser_wavelength Laser wavelength in nm. Default is 633.
#' @param angle Scattering angle in degrees. Default is 90.
#' @param solvent_ri Solvent refractive index. Default is 1.333 (water).
#' @param rg_col Optional Rg measure column from MALS for Rg/Rh ratio.
#' @param output_cols Output column names for Rh and diffusion coefficient.
#' @param role Role for generated columns.
#' @param trained Logical indicating if the step has been trained.
#' @param skip Logical. Should the step be skipped when baking?
#' @param id Unique step identifier.
#'
#' @return An updated recipe with the new step added.
#'
#' @details
#' The DLS correlation function g2(tau) is approximated as:
#' \deqn{g2(\tau) = 1 + \beta e^{-2 \Gamma \tau}}
#' where Gamma is the decay rate. The diffusion coefficient is calculated as
#' D = Gamma / q^2 with q determined from the scattering angle and wavelength.
#' Rh is then obtained from the Stokes-Einstein equation.
#'
#' @family sec-detectors
#' @export
#'
#' @examples
#' \dontrun{
#' library(recipes)
#' library(measure)
#'
#' rec <- recipe(~., data = sec_triple_detect) |>
#'   step_measure_input_long(dls_corr, location = vars(lag_time), col_name = "dls") |>
#'   step_sec_dls(measures = "dls", temperature = 25, viscosity = 0.89) |>
#'   prep()
#' }
step_sec_dls <- function(
  recipe,
  measures = NULL,
  temperature = 25,
  viscosity = NULL,
  laser_wavelength = 633,
  angle = 90,
  solvent_ri = 1.333,
  rg_col = "rg",
  output_cols = c("rh", "diffusion_coef"),
  role = NA,
  trained = FALSE,
  skip = FALSE,
  id = recipes::rand_id("sec_dls")
) {
  if (!is.numeric(temperature) || temperature <= -273.15) {
    cli::cli_abort("{.arg temperature} must be above absolute zero.")
  }

  if (!is.null(viscosity) && (!is.numeric(viscosity) || viscosity <= 0)) {
    cli::cli_abort("{.arg viscosity} must be a positive numeric value.")
  }

  if (!is.numeric(angle) || angle <= 0 || angle >= 180) {
    cli::cli_abort("{.arg angle} must be between 0 and 180 degrees.")
  }

  if (!is.character(output_cols)) {
    cli::cli_abort("{.arg output_cols} must be a character vector.")
  }

  if (length(output_cols) < 2) {
    cli::cli_abort("{.arg output_cols} must include Rh and diffusion columns.")
  }

  recipes::add_step(
    recipe,
    step_sec_dls_new(
      measures = measures,
      temperature = temperature,
      viscosity = viscosity,
      laser_wavelength = laser_wavelength,
      angle = angle,
      solvent_ri = solvent_ri,
      rg_col = rg_col,
      output_cols = output_cols,
      role = role,
      trained = trained,
      skip = skip,
      id = id
    )
  )
}

step_sec_dls_new <- function(
  measures,
  temperature,
  viscosity,
  laser_wavelength,
  angle,
  solvent_ri,
  rg_col,
  output_cols,
  role,
  trained,
  skip,
  id
) {
  recipes::step(
    subclass = "sec_dls",
    measures = measures,
    temperature = temperature,
    viscosity = viscosity,
    laser_wavelength = laser_wavelength,
    angle = angle,
    solvent_ri = solvent_ri,
    rg_col = rg_col,
    output_cols = output_cols,
    role = role,
    trained = trained,
    skip = skip,
    id = id
  )
}

#' @export
prep.step_sec_dls <- function(x, training, info = NULL, ...) {
  check_for_measure(training)

  measure_cols <- find_measure_cols(training)

  if (is.null(x$measures)) {
    dls_cols <- measure_cols[grepl(
      "dls|corr",
      measure_cols,
      ignore.case = TRUE
    )]
    if (length(dls_cols) == 0) {
      cli::cli_abort(
        "No DLS correlation column found. Specify {.arg measures} explicitly."
      )
    }
    measures <- dls_cols[1]
  } else {
    measures <- x$measures
  }

  if (length(measures) != 1) {
    cli::cli_abort("{.arg measures} must specify a single DLS column.")
  }

  if (!measures %in% measure_cols) {
    cli::cli_abort(
      "DLS column {.val {measures}} not found in measure columns."
    )
  }

  step_sec_dls_new(
    measures = measures,
    temperature = x$temperature,
    viscosity = x$viscosity,
    laser_wavelength = x$laser_wavelength,
    angle = x$angle,
    solvent_ri = x$solvent_ri,
    rg_col = x$rg_col,
    output_cols = x$output_cols,
    role = x$role,
    trained = TRUE,
    skip = x$skip,
    id = x$id
  )
}

#' @export
bake.step_sec_dls <- function(object, new_data, ...) {
  dls_col <- object$measures
  temperature <- object$temperature
  viscosity <- object$viscosity
  laser_wavelength <- object$laser_wavelength
  angle <- object$angle
  solvent_ri <- object$solvent_ri
  rg_col <- object$rg_col
  output_cols <- object$output_cols

  if (is.null(viscosity)) {
    viscosity <- .water_viscosity_pa_s(temperature)
  } else {
    viscosity <- viscosity * 1e-3
  }

  q_val <- .dls_scattering_vector(laser_wavelength, angle, solvent_ri)

  dls_results <- purrr::map(new_data[[dls_col]], function(m) {
    .dls_fit(m$location, m$value, q_val, temperature, viscosity)
  })

  rh_col <- output_cols[[1]]
  diff_col <- output_cols[[2]]

  new_data[[rh_col]] <- new_measure_list(
    purrr::map2(new_data[[dls_col]], dls_results, function(m, res) {
      new_measure_tbl(
        location = m$location,
        value = rep(res$rh_nm, length(m$location))
      )
    })
  )

  new_data[[diff_col]] <- new_measure_list(
    purrr::map2(new_data[[dls_col]], dls_results, function(m, res) {
      new_measure_tbl(
        location = m$location,
        value = rep(res$diffusion, length(m$location))
      )
    })
  )

  if (!is.null(rg_col) && rg_col %in% names(new_data)) {
    new_data[["rg_rh"]] <- new_measure_list(
      purrr::map2(new_data[[rg_col]], new_data[[rh_col]], function(rg_m, rh_m) {
        rh_at_rg <- stats::approx(
          x = rh_m$location,
          y = rh_m$value,
          xout = rg_m$location,
          rule = 2
        )$y
        ratio_values <- rg_m$value / rh_at_rg
        new_measure_tbl(location = rg_m$location, value = ratio_values)
      })
    )
  }

  tibble::as_tibble(new_data)
}

#' @export
print.step_sec_dls <- function(
  x,
  width = max(20, options()$width - 30),
  ...
) {
  title <- "SEC DLS processing"
  if (x$trained) {
    cat(title, " on ", x$measures, sep = "")
  } else {
    cat(title)
  }
  cat("\n")
  invisible(x)
}

#' @rdname tidy.step_sec
#' @export
#' @keywords internal
tidy.step_sec_dls <- function(x, ...) {
  tibble::tibble(
    dls_col = x$measures %||% NA_character_,
    temperature = x$temperature,
    viscosity = x$viscosity %||% NA_real_,
    laser_wavelength = x$laser_wavelength,
    angle = x$angle,
    output_cols = list(x$output_cols),
    id = x$id
  )
}

#' @rdname required_pkgs.step_sec
#' @export
#' @keywords internal
required_pkgs.step_sec_dls <- function(x, ...) {
  c("measure.sec", "measure")
}

.dls_scattering_vector <- function(laser_wavelength, angle, solvent_ri) {
  lambda_m <- laser_wavelength * 1e-9
  angle_rad <- angle * pi / 180
  (4 * pi * solvent_ri / lambda_m) * sin(angle_rad / 2)
}

.dls_fit <- function(tau, g2, q_val, temperature, viscosity_pa_s) {
  valid <- !is.na(tau) & !is.na(g2)
  tau <- tau[valid]
  g2 <- g2[valid]

  if (length(tau) < 5) {
    return(list(diffusion = NA_real_, rh_nm = NA_real_))
  }

  beta <- max(g2, na.rm = TRUE) - 1
  if (!is.finite(beta) || beta <= 0) {
    return(list(diffusion = NA_real_, rh_nm = NA_real_))
  }

  g1 <- sqrt(pmax((g2 - 1) / beta, 1e-12))

  order_idx <- order(tau)
  tau <- tau[order_idx]
  g1 <- g1[order_idx]

  n_fit <- max(5, floor(length(tau) * 0.3))
  tau_fit <- tau[seq_len(n_fit)]
  g1_fit <- g1[seq_len(n_fit)]

  fit <- stats::lm(log(g1_fit) ~ tau_fit)
  slope <- stats::coef(fit)[[2]]

  gamma <- -slope
  if (!is.finite(gamma) || gamma <= 0) {
    return(list(diffusion = NA_real_, rh_nm = NA_real_))
  }

  diffusion <- gamma / (q_val^2)

  k_b <- 1.380649e-23
  temp_k <- temperature + 273.15
  rh_m <- k_b * temp_k / (6 * pi * viscosity_pa_s * diffusion)

  list(diffusion = diffusion, rh_nm = rh_m * 1e9)
}

.water_viscosity_pa_s <- function(temperature_c) {
  # Vogel-Fulcher-Tammann equation approximation for water
  2.414e-5 * 10^(247.8 / (temperature_c + 133.15))
}
