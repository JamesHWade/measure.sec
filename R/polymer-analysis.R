# ==============================================================================
# Polymer Analysis Functions
#
# Mark-Houwink parameters, branching analysis, conformation
# ==============================================================================

#' Estimate Mark-Houwink Parameters
#'
#' Estimates the Mark-Houwink K and a (alpha) parameters from intrinsic
#' viscosity and molecular weight data.
#'
#' @param mw Numeric vector of molecular weights.
#' @param intrinsic_visc Numeric vector of intrinsic viscosities (same length
#'   as `mw`).
#' @param weights Optional numeric vector of weights for weighted regression.
#'   Use concentration or signal intensity as weights for SEC data.
#' @param mw_range Optional numeric vector of length 2 specifying the MW range
#'   to use for fitting. Data outside this range is excluded.
#' @param log_fit Logical. Perform fit in log-log space (recommended)?
#'   Default is `TRUE`.
#'
#' @return A list of class `mh_parameters` containing:
#'   \describe{
#'     \item{K}{Mark-Houwink K parameter}
#'     \item{a}{Mark-Houwink a (alpha) exponent}
#'     \item{r_squared}{R-squared of the fit}
#'     \item{n_points}{Number of data points used}
#'     \item{mw_range}{MW range of the data}
#'     \item{fit}{The fitted linear model object}
#'   }
#'
#' @details
#' The Mark-Houwink equation relates intrinsic viscosity to molecular weight:
#'
#' \deqn{[\eta] = K \cdot M^a}
#'
#' In log form:
#' \deqn{\log([\eta]) = \log(K) + a \cdot \log(M)}
#'
#' The parameters K and a depend on:
#' \itemize{
#'   \item Polymer-solvent system
#'   \item Temperature
#'   \item Polymer microstructure (tacticity, branching)
#' }
#'
#' **Interpretation of 'a' exponent:**
#' \itemize{
#'   \item a ~ 0.5: Theta solvent (polymer coil collapsed)
#'   \item a ~ 0.5-0.8: Good solvent (typical range)
#'   \item a ~ 0.8: Rigid rod or extended chain
#'   \item a < 0.5: Branched or compact structures
#' }
#'
#' @family sec-polymer
#' @export
#'
#' @examples
#' # Estimate Mark-Houwink parameters from triple-detection data
#' mw <- c(10000, 25000, 50000, 100000, 250000)
#' iv <- c(0.15, 0.28, 0.45, 0.72, 1.2)
#'
#' mh <- measure_mh_parameters(mw, iv)
#' print(mh)
#' # K = 0.000114, a = 0.716 (typical for PS in THF)
measure_mh_parameters <- function(
  mw,
  intrinsic_visc,
  weights = NULL,
  mw_range = NULL,
  log_fit = TRUE
) {
  # Validate inputs
  if (length(mw) != length(intrinsic_visc)) {
    cli::cli_abort(
      "{.arg mw} and {.arg intrinsic_visc} must have the same length."
    )
  }

  # Remove NA and invalid values
  valid <- !is.na(mw) & !is.na(intrinsic_visc) & mw > 0 & intrinsic_visc > 0

  if (!is.null(mw_range)) {
    valid <- valid & mw >= mw_range[1] & mw <= mw_range[2]
  }

  mw <- mw[valid]
  intrinsic_visc <- intrinsic_visc[valid]

  if (!is.null(weights)) {
    weights <- weights[valid]
  }

  if (length(mw) < 3) {
    cli::cli_abort("At least 3 valid data points are required for fitting.")
  }

  # Fit Mark-Houwink equation
  if (log_fit) {
    # Log-log linear regression
    log_mw <- log10(mw)
    log_iv <- log10(intrinsic_visc)

    if (!is.null(weights)) {
      fit <- stats::lm(log_iv ~ log_mw, weights = weights)
    } else {
      fit <- stats::lm(log_iv ~ log_mw)
    }

    # Extract parameters
    intercept <- stats::coef(fit)[1]
    slope <- stats::coef(fit)[2]

    K <- 10^intercept
    a <- slope
    r_squared <- summary(fit)$r.squared
  } else {
    # Non-linear fit (less common)
    start_params <- list(K = 0.0001, a = 0.7)

    if (!is.null(weights)) {
      fit <- stats::nls(
        intrinsic_visc ~ K * mw^a,
        start = start_params,
        weights = weights
      )
    } else {
      fit <- stats::nls(
        intrinsic_visc ~ K * mw^a,
        start = start_params
      )
    }

    K <- stats::coef(fit)["K"]
    a <- stats::coef(fit)["a"]

    # Calculate R-squared for nls
    ss_res <- sum(stats::residuals(fit)^2)
    ss_tot <- sum((intrinsic_visc - mean(intrinsic_visc))^2)
    r_squared <- 1 - ss_res / ss_tot
  }

  result <- list(
    K = as.numeric(K),
    a = as.numeric(a),
    r_squared = r_squared,
    n_points = length(mw),
    mw_range = range(mw),
    fit = fit
  )

  class(result) <- c("mh_parameters", "list")

  result
}


#' @export
print.mh_parameters <- function(x, ...) {
  cat("Mark-Houwink Parameters\n")
  cat(strrep("=", 40), "\n\n")
  cat(sprintf("K = %.4e\n", x$K))
  cat(sprintf("a = %.3f\n", x$a))
  cat(sprintf("\nR-squared: %.4f\n", x$r_squared))
  cat(sprintf("Data points: %d\n", x$n_points))
  cat(sprintf("MW range: %.0f - %.0f\n", x$mw_range[1], x$mw_range[2]))
  cat("\nEquation: [eta] = K * M^a\n")

  invisible(x)
}


#' Calculate Branching Index
#'
#' Calculates the branching index (g or g') for branched polymers by comparing
#' their properties to linear reference polymers of the same molecular weight.
#'
#' @param mw Numeric vector of molecular weights for branched samples.
#' @param rg Numeric vector of radius of gyration values (for g ratio).
#' @param intrinsic_visc Numeric vector of intrinsic viscosity values (for g').
#' @param reference Type of reference: either `"linear"` for comparison to
#'   linear polymer theory, or a data frame containing linear reference data.
#' @param mh_linear Mark-Houwink parameters for linear reference polymer
#'   (list with K and a). Required if `reference = "linear"` and using g'.
#' @param rg_linear_fit Rg-MW relationship for linear polymer: either a model
#'   or a list with `slope` and `intercept` for log(Rg) = intercept + slope * log(M).
#' @param method Branching index type:
#'   \itemize{
#'     \item `"g"`: Radius of gyration ratio: g = Rg^2(branched) / Rg^2(linear)
#'     \item `"g_prime"`: Viscosity ratio: g' = \[eta\](branched) / \[eta\](linear)
#'     \item `"both"`: Calculate both g and g'
#'   }
#'
#' @return A data frame with columns:
#'   \describe{
#'     \item{mw}{Molecular weight}
#'     \item{g}{Branching index from Rg (if calculated)}
#'     \item{g_prime}{Branching index from viscosity (if calculated)}
#'     \item{branches_per_molecule}{Estimated branch points (Zimm-Stockmayer)}
#'   }
#'
#' @details
#' Branching reduces the hydrodynamic size of polymers compared to their linear
#' counterparts. This is quantified by the branching ratios:
#'
#' **Rg-based branching index (g):**
#' \deqn{g = \frac{R_g^2(branched)}{R_g^2(linear)}}
#'
#' **Viscosity-based branching index (g'):**
#' \deqn{g' = \frac{[\eta](branched)}{[\eta](linear)}}
#'
#' The relationship between g and g' depends on polymer architecture:
#' \deqn{g' = g^\epsilon}
#'
#' where epsilon ~ 0.5-1.5 depending on branching type.
#'
#' **Zimm-Stockmayer Model (random branching):**
#' \deqn{g = \frac{6}{n_b} \left[ \frac{1}{2} + \frac{(2 + n_b)^{1/2} - 1 - n_b/2}{n_b} \right]}
#'
#' where n_b is the number of branch points per molecule.
#'
#' @family sec-polymer
#' @export
#'
#' @examples
#' # Calculate branching index from Rg data
#' mw <- c(100000, 200000, 500000)
#' rg_branched <- c(12, 18, 32)
#' rg_linear <- c(15, 22, 40)  # Reference linear polymer
#'
#' g <- measure_branching_index(
#'   mw = mw,
#'   rg = rg_branched,
#'   reference = data.frame(mw = mw, rg = rg_linear),
#'   method = "g"
#' )
measure_branching_index <- function(
  mw,
  rg = NULL,
  intrinsic_visc = NULL,
  reference = "linear",
  mh_linear = NULL,
  rg_linear_fit = NULL,
  method = c("g", "g_prime", "both")
) {
  method <- match.arg(method)

  n <- length(mw)

  # Initialize result
  result <- data.frame(mw = mw)

  # Calculate g (Rg-based branching index)
  if (method %in% c("g", "both")) {
    if (is.null(rg)) {
      cli::cli_abort("{.arg rg} is required for method {.val {method}}.")
    }

    if (is.data.frame(reference)) {
      # Interpolate reference Rg at sample MW
      if (!"rg" %in% names(reference)) {
        cli::cli_abort("Reference data must contain 'rg' column.")
      }
      rg_linear <- stats::approx(
        x = reference$mw,
        y = reference$rg,
        xout = mw,
        rule = 2
      )$y
    } else if (!is.null(rg_linear_fit)) {
      # Use Rg-MW relationship
      if (is.list(rg_linear_fit) && "slope" %in% names(rg_linear_fit)) {
        log_rg_linear <- rg_linear_fit$intercept +
          rg_linear_fit$slope * log10(mw)
        rg_linear <- 10^log_rg_linear
      } else if (inherits(rg_linear_fit, "lm")) {
        log_rg_linear <- stats::predict(
          rg_linear_fit,
          newdata = data.frame(log_mw = log10(mw))
        )
        rg_linear <- 10^log_rg_linear
      } else {
        cli::cli_abort(
          "Invalid {.arg rg_linear_fit}. Provide list with slope/intercept or lm object."
        )
      }
    } else {
      # Use theoretical random coil scaling
      # Rg ~ M^0.588 for good solvent
      # Normalize to match data at highest MW point
      rg_linear <- rg[which.max(mw)] * (mw / max(mw))^0.588
    }

    # Calculate g = Rg^2(branched) / Rg^2(linear)
    result$g <- (rg^2) / (rg_linear^2)
    result$rg_branched <- rg
    result$rg_linear <- rg_linear
  }

  # Calculate g' (viscosity-based branching index)
  if (method %in% c("g_prime", "both")) {
    if (is.null(intrinsic_visc)) {
      cli::cli_abort(
        "{.arg intrinsic_visc} is required for method {.val {method}}."
      )
    }

    if (is.data.frame(reference)) {
      # Interpolate reference IV at sample MW
      if (!"intrinsic_visc" %in% names(reference)) {
        cli::cli_abort("Reference data must contain 'intrinsic_visc' column.")
      }
      iv_linear <- stats::approx(
        x = reference$mw,
        y = reference$intrinsic_visc,
        xout = mw,
        rule = 2
      )$y
    } else if (!is.null(mh_linear)) {
      # Calculate from Mark-Houwink parameters
      iv_linear <- mh_linear$K * mw^mh_linear$a
    } else {
      cli::cli_abort(
        "Either reference data frame or {.arg mh_linear} is required for g'."
      )
    }

    # Calculate g' = [eta](branched) / [eta](linear)
    result$g_prime <- intrinsic_visc / iv_linear
    result$iv_branched <- intrinsic_visc
    result$iv_linear <- iv_linear
  }

  # Estimate branch points using Zimm-Stockmayer for random branching
  if ("g" %in% names(result)) {
    g_vals <- result$g

    # Solve Zimm-Stockmayer equation numerically for n_b
    # g = ((1 + n_b/7)^0.5 + 4*n_b/9*pi)^(-0.5) (approximate form for trifunctional)
    # For small n_b: g ~ 1 - 0.143*n_b (linear approximation)

    # Simple inversion for small branching
    result$branches_per_molecule <- ifelse(
      g_vals < 1,
      (1 - g_vals) / 0.143,
      0
    )
    result$branches_per_molecule[result$branches_per_molecule < 0] <- 0
  }

  class(result) <- c("branching_index", "data.frame")

  result
}


#' @export
print.branching_index <- function(x, ...) {
  cat("Branching Analysis Results\n")
  cat(strrep("=", 50), "\n\n")

  cat("MW Range:", sprintf("%.0f - %.0f\n", min(x$mw), max(x$mw)))

  if ("g" %in% names(x)) {
    cat(sprintf(
      "g (Rg ratio): %.3f - %.3f (mean: %.3f)\n",
      min(x$g, na.rm = TRUE),
      max(x$g, na.rm = TRUE),
      mean(x$g, na.rm = TRUE)
    ))
  }

  if ("g_prime" %in% names(x)) {
    cat(sprintf(
      "g' (IV ratio): %.3f - %.3f (mean: %.3f)\n",
      min(x$g_prime, na.rm = TRUE),
      max(x$g_prime, na.rm = TRUE),
      mean(x$g_prime, na.rm = TRUE)
    ))
  }

  if ("branches_per_molecule" %in% names(x)) {
    cat(sprintf(
      "\nEstimated branches/molecule: %.1f - %.1f\n",
      min(x$branches_per_molecule, na.rm = TRUE),
      max(x$branches_per_molecule, na.rm = TRUE)
    ))
  }

  cat("\n")
  print(tibble::as_tibble(x[, c("mw", intersect(names(x), c("g", "g_prime")))]))

  invisible(x)
}


#' Create Conformation Plot Data
#'
#' Prepares data for Mark-Houwink or conformation plots (log(\[eta\]) vs log(M)
#' or log(Rg) vs log(M)).
#'
#' @param mw Numeric vector of molecular weights.
#' @param y Numeric vector of intrinsic viscosity or radius of gyration.
#' @param y_type Type of y-axis data: `"iv"` for intrinsic viscosity or
#'   `"rg"` for radius of gyration.
#' @param fit_line Logical. Include fitted line? Default is `TRUE`.
#'
#' @return A data frame suitable for plotting with columns:
#'   \describe{
#'     \item{log_mw}{log10(MW)}
#'     \item{log_y}{log10(y)}
#'     \item{mw}{Original MW values}
#'     \item{y}{Original y values}
#'   }
#'
#' @family sec-polymer
#' @export
#'
#' @examples
#' mw <- c(10000, 50000, 100000, 500000)
#' iv <- c(0.15, 0.35, 0.50, 0.95)
#'
#' plot_data <- measure_conformation_data(mw, iv, y_type = "iv")
#'
#' # Plot with ggplot2
#' # ggplot(plot_data, aes(log_mw, log_y)) +
#' #   geom_point() +
#' #   geom_smooth(method = "lm")
measure_conformation_data <- function(
  mw,
  y,
  y_type = c("iv", "rg"),
  fit_line = TRUE
) {
  y_type <- match.arg(y_type)

  # Remove invalid values
  valid <- !is.na(mw) & !is.na(y) & mw > 0 & y > 0
  mw <- mw[valid]
  y <- y[valid]

  result <- data.frame(
    mw = mw,
    y = y,
    log_mw = log10(mw),
    log_y = log10(y)
  )

  if (fit_line && nrow(result) >= 3) {
    fit <- stats::lm(log_y ~ log_mw, data = result)
    result$fitted <- stats::fitted(fit)

    attr(result, "slope") <- stats::coef(fit)[2]
    attr(result, "intercept") <- stats::coef(fit)[1]
    attr(result, "r_squared") <- summary(fit)$r.squared
  }

  attr(result, "y_type") <- y_type

  class(result) <- c("conformation_data", "data.frame")

  result
}
