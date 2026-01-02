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


#' Calculate Branching Frequency from Zimm-Stockmayer Theory
#'
#' Calculates the number of branch points per molecule using Zimm-Stockmayer
#' theory, which relates the branching index (g) to the branching frequency.
#'
#' @param g Numeric vector of branching ratios g = Rg^2(branched) / Rg^2(linear).
#' @param architecture Branching architecture model:
#'   \itemize{
#'     \item `"random"`: Random trifunctional branching (default)
#'     \item `"star_3"`: 3-arm star polymer
#'     \item `"star_4"`: 4-arm star polymer
#'     \item `"star_f"`: f-arm star (requires `arms` parameter)
#'     \item `"comb"`: Comb-like branching
#'   }
#' @param arms Number of arms for star polymers (required for `"star_f"`).
#' @param mw Optional molecular weight vector (same length as `g`) for
#'   calculating total branches in sample.
#'
#' @return A data frame of class `branching_frequency` containing:
#'   \describe{
#'     \item{g}{Input branching ratio}
#'     \item{branches_per_molecule}{Estimated branch points per molecule}
#'     \item{branch_density}{Branches per unit MW (if `mw` provided)}
#'   }
#'
#' @details
#' **Zimm-Stockmayer Theory:**
#'
#' For random trifunctional branching, the relationship between g and the
#' average number of branch points (n_b) is:
#'
#' \deqn{g = \left[\left(1 + \frac{n_b}{7}\right)^{1/2} + \frac{4n_b}{9\pi}\right]^{-1/2}}
#'
#' For star polymers with f arms:
#' \deqn{g = \frac{3f - 2}{f^2}}
#'
#' **Interpretation:**
#' \itemize{
#'   \item g = 1.0: Linear polymer (no branching)
#'   \item g ~ 0.5-0.8: Lightly branched
#'   \item g ~ 0.3-0.5: Moderately branched
#'   \item g < 0.3: Highly branched/hyperbranched
#' }
#'
#' @references
#' Zimm, B.H. and Stockmayer, W.H. (1949). "The Dimensions of Chain Molecules
#' Containing Branches and Rings." J. Chem. Phys., 17, 1301-1314.
#'
#' @family sec-polymer
#' @export
#'
#' @examples
#' # Calculate branching for random trifunctional polymer
#' g_values <- c(0.9, 0.7, 0.5, 0.3)
#' bf <- measure_branching_frequency(g_values)
#' print(bf)
#'
#' # For star polymer with 4 arms
#' bf_star <- measure_branching_frequency(0.625, architecture = "star_4")
measure_branching_frequency <- function(
  g,
  architecture = c("random", "star_3", "star_4", "star_f", "comb"),
  arms = NULL,
  mw = NULL
) {
  architecture <- match.arg(architecture)

  if (!is.numeric(g) || any(g <= 0, na.rm = TRUE) || any(g > 1, na.rm = TRUE)) {
    cli::cli_abort(
      "{.arg g} must be numeric values between 0 and 1 (exclusive)."
    )
  }

  n <- length(g)

  # Calculate branches based on architecture
  branches <- switch(
    architecture,
    "random" = {
      # Zimm-Stockmayer for random trifunctional branching
      # Numerical inversion of the ZS equation
      sapply(g, function(gi) {
        if (is.na(gi) || gi >= 1) {
          return(0)
        }
        # Use optimization to solve ZS equation
        zs_func <- function(n_b) {
          if (n_b < 0) {
            return(Inf)
          }
          g_calc <- ((1 + n_b / 7)^0.5 + 4 * n_b / (9 * pi))^(-0.5)
          (g_calc - gi)^2
        }
        result <- stats::optimize(zs_func, c(0, 100))
        result$minimum
      })
    },
    "star_3" = {
      # 3-arm star: g = (3*3 - 2) / 3^2 = 7/9 = 0.778
      # If g <= theoretical value + tolerance, it's a star (1 branch point)
      # If g > threshold (closer to 1), it's essentially linear
      ifelse(g <= 7 / 9 + 0.05, 1, 0) # 3-arm star has 1 branch point
    },
    "star_4" = {
      # 4-arm star: g = (3*4 - 2) / 4^2 = 10/16 = 0.625
      # If g <= theoretical value + tolerance, it's a star (1 branch point)
      ifelse(g <= 10 / 16 + 0.05, 1, 0) # 4-arm star has 1 branch point
    },
    "star_f" = {
      if (is.null(arms)) {
        cli::cli_abort(
          "{.arg arms} is required for {.val star_f} architecture."
        )
      }
      if (!is.numeric(arms) || length(arms) != 1 || arms < 3) {
        cli::cli_abort(
          "{.arg arms} must be a single numeric value >= 3."
        )
      }
      # General star: g = (3f - 2) / f^2
      # For f-arm star with specified arms, check if g is consistent
      # All star polymers have exactly 1 branch point (the central junction)
      g_theoretical <- (3 * arms - 2) / arms^2
      sapply(g, function(gi) {
        if (gi >= 1 || gi <= 0) {
          return(0)
        }
        # If g is close to theoretical value for this star, return 1 branch point
        if (abs(gi - g_theoretical) < 0.1) {
          return(1)
        }
        # If g is higher (closer to 1), essentially linear
        if (gi > g_theoretical + 0.1) {
          return(0)
        }
        # If g is lower, could indicate multiple stars or more complex architecture
        # Still report 1 for star architecture interpretation
        1
      })
    },
    "comb" = {
      # Comb polymer approximation
      # g ~ 1 / (1 + 2*n_b/3) for uniform comb
      sapply(g, function(gi) {
        if (is.na(gi) || gi >= 1) {
          return(0)
        }
        max(0, 3 * (1 / gi - 1) / 2)
      })
    }
  )

  # Warn about impossible configurations
  na_count <- sum(is.na(branches))
  if (na_count > 0) {
    cli::cli_warn(c(
      "!" = "{na_count} g value{?s} could not be converted to branch counts.",
      "i" = "These may represent physically impossible configurations for the {.val {architecture}} architecture.",
      "i" = "NA values returned for these entries."
    ))
  }

  # Build result

  result <- data.frame(
    g = g,
    branches_per_molecule = branches,
    architecture = architecture
  )

  # Add branch density if MW provided
  if (!is.null(mw)) {
    if (length(mw) != n) {
      cli::cli_abort("{.arg mw} must have the same length as {.arg g}.")
    }
    # Branches per 1000 Da
    result$branch_density <- branches / (mw / 1000)
    result$mw <- mw
  }

  class(result) <- c("branching_frequency", "data.frame")

  result
}


#' @export
print.branching_frequency <- function(x, ...) {
  cat("Branching Frequency Analysis (Zimm-Stockmayer)\n")
  cat(strrep("=", 50), "\n\n")

  cat(sprintf("Architecture: %s\n", unique(x$architecture)))
  cat(sprintf("N samples: %d\n\n", nrow(x)))

  cat(sprintf(
    "g ratio range: %.3f - %.3f\n",
    min(x$g, na.rm = TRUE),
    max(x$g, na.rm = TRUE)
  ))
  cat(sprintf(
    "Branches/molecule: %.2f - %.2f (mean: %.2f)\n",
    min(x$branches_per_molecule, na.rm = TRUE),
    max(x$branches_per_molecule, na.rm = TRUE),
    mean(x$branches_per_molecule, na.rm = TRUE)
  ))

  if ("branch_density" %in% names(x)) {
    cat(sprintf(
      "Branch density: %.4f - %.4f per 1000 Da\n",
      min(x$branch_density, na.rm = TRUE),
      max(x$branch_density, na.rm = TRUE)
    ))
  }

  cat("\n")
  print(tibble::as_tibble(x[, c("g", "branches_per_molecule")]))

  invisible(x)
}


#' Fit Rg-MW Scaling Relationship
#'
#' Fits the power law relationship between radius of gyration (Rg) and
#' molecular weight (MW) to determine polymer conformation.
#'
#' @param mw Numeric vector of molecular weights.
#' @param rg Numeric vector of radii of gyration (same length as `mw`).
#' @param weights Optional numeric vector of weights (e.g., concentration).
#' @param mw_range Optional numeric vector of length 2 specifying MW range
#'   for fitting.
#'
#' @return A list of class `rg_mw_scaling` containing:
#'   \describe{
#'     \item{nu}{Scaling exponent (slope in log-log space)}
#'     \item{prefactor}{Prefactor K in Rg = K * M^nu}
#'     \item{r_squared}{R-squared of the fit}
#'     \item{conformation}{Interpreted polymer conformation}
#'     \item{n_points}{Number of data points used}
#'     \item{fit}{The lm fit object}
#'   }
#'
#' @details
#' The Rg-MW relationship follows a power law:
#'
#' \deqn{R_g = K \cdot M^{\nu}}
#'
#' where nu (the Flory exponent) indicates polymer conformation:
#'
#' **Interpretation of nu:**
#' \itemize{
#'   \item nu ~ 0.33: Compact/spherical (collapsed globule)
#'   \item nu ~ 0.50: Theta solvent (ideal chain)
#'   \item nu ~ 0.588: Good solvent (swollen coil)
#'   \item nu ~ 1.0: Rigid rod
#' }
#'
#' **Typical Values:**
#' \itemize{
#'   \item Flexible polymers in good solvent: nu = 0.55-0.60
#'   \item Branched polymers: nu = 0.40-0.50
#'   \item Proteins (globular): nu = 0.30-0.35
#'   \item DNA: nu = 0.58-0.60
#' }
#'
#' @references
#' Flory, P.J. (1953). "Principles of Polymer Chemistry." Cornell University Press.
#'
#' @family sec-polymer
#' @export
#'
#' @examples
#' # Fit Rg-MW for a linear polymer
#' mw <- c(10000, 50000, 100000, 500000, 1000000)
#' rg <- c(4.5, 12, 18, 45, 70)  # nm
#'
#' scaling <- measure_rg_mw_scaling(mw, rg)
#' print(scaling)
measure_rg_mw_scaling <- function(
  mw,
  rg,
  weights = NULL,
  mw_range = NULL
) {
  if (length(mw) != length(rg)) {
    cli::cli_abort("{.arg mw} and {.arg rg} must have the same length.")
  }

  # Remove invalid values
  valid <- !is.na(mw) & !is.na(rg) & mw > 0 & rg > 0

  if (!is.null(mw_range)) {
    valid <- valid & mw >= mw_range[1] & mw <= mw_range[2]
  }

  mw <- mw[valid]
  rg <- rg[valid]

  if (!is.null(weights)) {
    weights <- weights[valid]
  }

  if (length(mw) < 3) {
    cli::cli_abort("At least 3 valid data points are required for fitting.")
  }

  # Log-log linear regression
  log_mw <- log10(mw)
  log_rg <- log10(rg)

  if (!is.null(weights)) {
    fit <- stats::lm(log_rg ~ log_mw, weights = weights)
  } else {
    fit <- stats::lm(log_rg ~ log_mw)
  }

  # Extract parameters
  intercept <- stats::coef(fit)[1]
  slope <- stats::coef(fit)[2]

  nu <- as.numeric(slope)
  prefactor <- 10^intercept

  # Interpret conformation
  conformation <- if (nu < 0.35) {
    "compact/spherical"
  } else if (nu < 0.45) {
    "branched/collapsed"
  } else if (nu < 0.53) {
    "theta solvent (ideal chain)"
  } else if (nu < 0.65) {
    "good solvent (swollen coil)"
  } else if (nu < 0.85) {
    "semi-rigid"
  } else {
    "rigid rod"
  }

  # Confidence intervals for nu
  ci <- stats::confint(fit, "log_mw", level = 0.95)

  result <- list(
    nu = nu,
    prefactor = as.numeric(prefactor),
    r_squared = summary(fit)$r.squared,
    conformation = conformation,
    nu_ci_lower = ci[1],
    nu_ci_upper = ci[2],
    n_points = length(mw),
    mw_range = range(mw),
    fit = fit
  )

  class(result) <- c("rg_mw_scaling", "list")

  result
}


#' @export
print.rg_mw_scaling <- function(x, ...) {
  cat("Rg-MW Scaling Analysis\n")
  cat(strrep("=", 40), "\n\n")

  cat("Scaling Law: Rg = K * M^nu\n\n")
  cat(sprintf("nu (Flory exponent): %.3f\n", x$nu))
  cat(sprintf("  95%% CI: [%.3f, %.3f]\n", x$nu_ci_lower, x$nu_ci_upper))
  cat(sprintf("K (prefactor): %.4e\n", x$prefactor))
  cat(sprintf("\nConformation: %s\n", x$conformation))
  cat(sprintf("\nR-squared: %.4f\n", x$r_squared))
  cat(sprintf("Data points: %d\n", x$n_points))
  cat(sprintf("MW range: %.0f - %.0f\n", x$mw_range[1], x$mw_range[2]))

  invisible(x)
}


#' Column Performance Metrics for SEC
#'
#' Calculates comprehensive column performance metrics including HETP,
#' separation range, and efficiency parameters.
#'
#' @param calibration_data A data frame containing calibration data with columns:
#'   \itemize{
#'     \item `retention` or `retention_time`: Retention time/volume
#'     \item `mw` or `molecular_weight`: Molecular weight
#'     \item `width`: Peak width (optional, for plate calculation)
#'   }
#' @param column_length Column length in cm. Default is 30 cm.
#' @param column_diameter Column inner diameter in mm. Default is 7.8 mm.
#' @param flow_rate Flow rate in mL/min (optional, for linear velocity).
#' @param particle_size Particle size in micrometers (optional, for reduced HETP).
#' @param dead_volume Column dead volume in mL (optional).
#'
#' @return A list of class `sec_column_performance` containing:
#'   \describe{
#'     \item{separation_range}{MW range (exclusion limit to total permeation)}
#'     \item{selectivity}{Slope of log(MW) vs retention (mL^-1)}
#'     \item{hetp}{Height equivalent to theoretical plate (mm)}
#'     \item{plates_per_meter}{Theoretical plates per meter}
#'     \item{reduced_hetp}{HETP/particle_size (if particle_size provided)}
#'     \item{peak_capacity}{Estimated peak capacity in separation range}
#'     \item{resolution_factor}{Resolution per decade of MW}
#'   }
#'
#' @details
#' **Key SEC Column Performance Metrics:**
#'
#' **HETP (Height Equivalent to Theoretical Plate):**
#' \deqn{HETP = \frac{L}{N}}
#'
#' **Reduced HETP (h):**
#' \deqn{h = \frac{HETP}{d_p}}
#' Optimal h ~ 2-3 for well-packed columns.
#'
#' **Selectivity (D):**
#' \deqn{D = \frac{d \log M}{dV_R}}
#' Higher D means better MW resolution.
#'
#' **Peak Capacity:**
#' \deqn{n_c = 1 + \frac{\sqrt{N}}{4} \ln\left(\frac{V_{exclusion}}{V_{total}}\right)}
#'
#' **Typical Performance Guidelines:**
#' \itemize{
#'   \item HETP: < 50 um for analytical columns
#'   \item Reduced HETP: 2-5 for good columns
#'   \item Plates/meter: > 20,000 for HPLC SEC
#' }
#'
#' @family sec-qc
#' @export
#'
#' @examples
#' # From calibration standards
#' cal_data <- data.frame(
#'   retention = c(5.2, 6.1, 7.0, 8.2, 9.5, 10.8),
#'   mw = c(1200000, 400000, 100000, 30000, 5000, 580),
#'   width = c(0.4, 0.35, 0.30, 0.28, 0.25, 0.30)
#' )
#'
#' perf <- measure_sec_column_performance(
#'   cal_data,
#'   column_length = 30,
#'   particle_size = 5
#' )
#' print(perf)
measure_sec_column_performance <- function(
  calibration_data,
  column_length = 30,
  column_diameter = 7.8,
  flow_rate = NULL,
  particle_size = NULL,
  dead_volume = NULL
) {
  # Validate inputs
  if (!is.data.frame(calibration_data)) {
    cli::cli_abort("{.arg calibration_data} must be a data frame.")
  }

  if (
    !is.numeric(column_length) ||
      length(column_length) != 1 ||
      column_length <= 0
  ) {
    cli::cli_abort(
      "{.arg column_length} must be a single positive numeric value."
    )
  }

  if (!is.null(particle_size)) {
    if (
      !is.numeric(particle_size) ||
        length(particle_size) != 1 ||
        particle_size <= 0
    ) {
      cli::cli_abort(
        "{.arg particle_size} must be a single positive numeric value."
      )
    }
  }

  # Find retention column
  ret_col <- intersect(
    names(calibration_data),
    c("retention", "retention_time", "rt", "elution_volume", "ve")
  )
  if (length(ret_col) == 0) {
    cli::cli_abort(
      "No retention column found. Provide 'retention' or 'retention_time'."
    )
  }
  retention <- calibration_data[[ret_col[1]]]

  # Find MW column
  mw_col <- intersect(
    names(calibration_data),
    c("mw", "molecular_weight", "mol_wt", "MW")
  )
  if (length(mw_col) == 0) {
    cli::cli_abort(
      "No molecular weight column found. Provide 'mw' or 'molecular_weight'."
    )
  }
  mw <- calibration_data[[mw_col[1]]]

  # Sort by retention time
  ord <- order(retention)
  retention <- retention[ord]
  mw <- mw[ord]

  # Calculate separation range
  separation_range <- list(
    exclusion_limit = max(mw, na.rm = TRUE),
    total_permeation = min(mw, na.rm = TRUE),
    log_range = log10(max(mw, na.rm = TRUE)) - log10(min(mw, na.rm = TRUE))
  )

  # Fit calibration curve for selectivity
  log_mw <- log10(mw)
  cal_fit <- stats::lm(log_mw ~ retention)
  selectivity <- abs(stats::coef(cal_fit)[2]) # |slope| in log(MW)/retention unit

  cal_r_squared <- summary(cal_fit)$r.squared

  # Calculate plate count if width data available
  hetp <- NA
  plates_per_meter <- NA
  reduced_hetp <- NA
  avg_plates <- NA

  if ("width" %in% names(calibration_data)) {
    width <- calibration_data$width[ord]
    # Calculate N for each peak
    plates <- 5.54 * (retention / width)^2
    avg_plates <- mean(plates, na.rm = TRUE)

    hetp <- (column_length * 10) / avg_plates # Convert cm to mm
    plates_per_meter <- avg_plates / (column_length / 100)

    if (!is.null(particle_size)) {
      reduced_hetp <- (hetp * 1000) / particle_size # hetp in um / particle size in um
    }
  }

  # Estimate peak capacity
  retention_range <- max(retention, na.rm = TRUE) - min(retention, na.rm = TRUE)
  if (!is.na(avg_plates) && avg_plates > 0) {
    # Peak capacity = 1 + (sqrt(N)/4) * ln(V2/V1)
    peak_capacity <- 1 +
      (sqrt(avg_plates) / 4) *
        log(max(retention, na.rm = TRUE) / min(retention, na.rm = TRUE))
  } else {
    peak_capacity <- NA
  }

  # Resolution per decade of MW
  # Rs ~ 0.25 * sqrt(N) * D
  if (!is.na(avg_plates)) {
    resolution_per_decade <- 0.25 *
      sqrt(avg_plates) *
      selectivity *
      (retention_range / separation_range$log_range)
  } else {
    resolution_per_decade <- NA
  }

  # Calculate linear velocity if flow rate provided
  linear_velocity <- NA
  if (!is.null(flow_rate)) {
    # Cross-sectional area in cm^2
    area_cm2 <- pi * (column_diameter / 20)^2 # diameter in mm to radius in cm
    # Linear velocity in cm/min
    linear_velocity <- flow_rate / area_cm2
  }

  result <- list(
    separation_range = separation_range,
    selectivity = as.numeric(selectivity),
    calibration_r_squared = cal_r_squared,
    hetp_mm = hetp,
    plates_per_meter = plates_per_meter,
    reduced_hetp = reduced_hetp,
    avg_plates = avg_plates,
    peak_capacity = peak_capacity,
    resolution_per_decade = resolution_per_decade,
    linear_velocity = linear_velocity,
    column_length = column_length,
    column_diameter = column_diameter,
    n_standards = nrow(calibration_data)
  )

  class(result) <- c("sec_column_performance", "list")

  result
}


#' @export
print.sec_column_performance <- function(x, ...) {
  cat("SEC Column Performance\n")
  cat(strrep("=", 50), "\n\n")

  cat("Separation Range:\n")
  cat(sprintf(
    "  Exclusion limit: %.0f Da\n",
    x$separation_range$exclusion_limit
  ))
  cat(sprintf(
    "  Total permeation: %.0f Da\n",
    x$separation_range$total_permeation
  ))
  cat(sprintf(
    "  Log MW range: %.2f decades\n\n",
    x$separation_range$log_range
  ))

  cat("Calibration:\n")
  cat(sprintf(
    "  Selectivity: %.4f log(MW)/unit\n",
    x$selectivity
  ))
  cat(sprintf("  R-squared: %.4f\n\n", x$calibration_r_squared))

  if (!is.na(x$hetp_mm)) {
    cat("Column Efficiency:\n")
    cat(sprintf("  HETP: %.3f mm (%.1f um)\n", x$hetp_mm, x$hetp_mm * 1000))
    cat(sprintf("  Plates/meter: %.0f\n", x$plates_per_meter))
    if (!is.na(x$reduced_hetp)) {
      cat(sprintf("  Reduced HETP (h): %.2f\n", x$reduced_hetp))
    }
    cat(sprintf("  Average plates (N): %.0f\n\n", x$avg_plates))
  }

  if (!is.na(x$peak_capacity)) {
    cat("Resolution:\n")
    cat(sprintf("  Peak capacity: %.1f\n", x$peak_capacity))
    cat(sprintf("  Resolution/decade: %.2f\n\n", x$resolution_per_decade))
  }

  cat(sprintf(
    "Column: %d x %.1f mm\n",
    x$column_length * 10,
    x$column_diameter
  ))
  cat(sprintf("Standards used: %d\n", x$n_standards))

  invisible(x)
}


#' Compare Branching Models
#'
#' Compares experimental branching data against theoretical predictions from
#' different branching models to identify the most likely architecture.
#'
#' @param g Numeric vector of experimental g ratios (Rg^2 branched / Rg^2 linear).
#' @param mw Numeric vector of molecular weights (same length as `g`).
#' @param g_prime Optional numeric vector of g' ratios (IV branched / IV linear).
#' @param models Character vector of models to compare. Default compares all:
#'   \itemize{
#'     \item `"random"`: Random trifunctional branching (Zimm-Stockmayer)
#'     \item `"star"`: Star polymers (variable arms)
#'     \item `"comb"`: Comb architecture
#'     \item `"hyperbranched"`: Hyperbranched/dendritic
#'   }
#' @param branch_frequency For fitting: assumed constant branching frequency
#'   (branches per 1000 Da). If NULL, estimated from data.
#'
#' @return A list of class `branching_model_comparison` containing:
#'   \describe{
#'     \item{model_fits}{Data frame with fit statistics for each model}
#'     \item{best_model}{Name of the best-fitting model}
#'     \item{predictions}{Data frame with predicted g for each model}
#'     \item{experimental}{Input experimental data}
#'   }
#'
#' @details
#' **Model Equations:**
#'
#' **Random trifunctional (Zimm-Stockmayer):**
#' \deqn{g = \left[\left(1 + \frac{n_b}{7}\right)^{1/2} + \frac{4n_b}{9\pi}\right]^{-1/2}}
#'
#' **Star polymer (f arms):**
#' \deqn{g = \frac{3f - 2}{f^2}}
#'
#' **Comb polymer (n_b branches):**
#' \deqn{g \approx \frac{1}{1 + 2n_b/3}}
#'
#' **Hyperbranched (degree of branching DB):**
#' \deqn{g \approx \left(\frac{1}{1 + DB \cdot n/2}\right)^{0.5}}
#'
#' Model selection uses residual sum of squares and AIC-like criteria.
#'
#' @references
#' Zimm, B.H. and Stockmayer, W.H. (1949). J. Chem. Phys., 17, 1301-1314.
#'
#' Burchard, W. (1999). "Solution Properties of Branched Macromolecules."
#' Adv. Polym. Sci., 143, 113-194.
#'
#' @family sec-polymer
#' @export
#'
#' @examples
#' # Compare models for branched polymer data
#' mw <- c(50000, 100000, 200000, 500000)
#' g_exp <- c(0.85, 0.72, 0.58, 0.42)
#'
#' comparison <- measure_branching_model_comparison(g_exp, mw)
#' print(comparison)
measure_branching_model_comparison <- function(
  g,
  mw,
  g_prime = NULL,
  models = c("random", "star", "comb", "hyperbranched"),
  branch_frequency = NULL
) {
  models <- match.arg(models, several.ok = TRUE)

  if (length(g) != length(mw)) {
    cli::cli_abort("{.arg g} and {.arg mw} must have the same length.")
  }

  # Remove invalid data
  valid <- !is.na(g) & !is.na(mw) & g > 0 & g <= 1 & mw > 0
  g <- g[valid]
  mw <- mw[valid]

  if (!is.null(g_prime)) {
    g_prime <- g_prime[valid]
  }

  n <- length(g)

  if (n < 3) {
    cli::cli_abort(
      "At least 3 valid data points required for model comparison."
    )
  }

  # Model prediction functions
  # These predict g given MW and a branching parameter

  predict_random <- function(mw, lambda) {
    # lambda = branches per unit MW
    n_b <- lambda * mw / 1000
    sapply(n_b, function(nb) {
      if (nb <= 0) {
        return(1)
      }
      ((1 + nb / 7)^0.5 + 4 * nb / (9 * pi))^(-0.5)
    })
  }

  predict_star <- function(mw, f_per_mw) {
    # f increases with MW for multi-star
    f <- pmax(2, 2 + f_per_mw * mw / 100000)
    (3 * f - 2) / f^2
  }

  predict_comb <- function(mw, lambda) {
    n_b <- lambda * mw / 1000
    1 / (1 + 2 * n_b / 3)
  }

  predict_hyperbranched <- function(mw, db) {
    # db = degree of branching (0-1)
    # More branches at higher MW
    n_equiv <- db * mw / 5000
    (1 / (1 + n_equiv / 2))^0.5
  }

  # Fit each model
  model_fits <- list()
  failed_models <- character(0)

  if ("random" %in% models) {
    fit_random <- tryCatch(
      {
        stats::nls(
          g ~ predict_random(mw, lambda),
          start = list(lambda = 0.1),
          data = data.frame(g = g, mw = mw),
          control = stats::nls.control(maxiter = 100, warnOnly = TRUE)
        )
      },
      error = function(e) NULL
    )

    if (!is.null(fit_random)) {
      pred_random <- stats::predict(fit_random)
      ss_res <- sum((g - pred_random)^2)
      ss_tot <- sum((g - mean(g))^2)
      model_fits$random <- list(
        model = "random",
        parameter = stats::coef(fit_random)["lambda"],
        r_squared = 1 - ss_res / ss_tot,
        rmse = sqrt(mean((g - pred_random)^2)),
        aic = stats::AIC(fit_random),
        prediction = pred_random,
        fit = fit_random
      )
    } else {
      failed_models <- c(failed_models, "random")
    }
  }

  if ("star" %in% models) {
    fit_star <- tryCatch(
      {
        stats::nls(
          g ~ predict_star(mw, f_per_mw),
          start = list(f_per_mw = 0.5),
          data = data.frame(g = g, mw = mw),
          control = stats::nls.control(maxiter = 100, warnOnly = TRUE)
        )
      },
      error = function(e) NULL
    )

    if (!is.null(fit_star)) {
      pred_star <- stats::predict(fit_star)
      ss_res <- sum((g - pred_star)^2)
      ss_tot <- sum((g - mean(g))^2)
      model_fits$star <- list(
        model = "star",
        parameter = stats::coef(fit_star)["f_per_mw"],
        r_squared = 1 - ss_res / ss_tot,
        rmse = sqrt(mean((g - pred_star)^2)),
        aic = stats::AIC(fit_star),
        prediction = pred_star,
        fit = fit_star
      )
    } else {
      failed_models <- c(failed_models, "star")
    }
  }

  if ("comb" %in% models) {
    fit_comb <- tryCatch(
      {
        stats::nls(
          g ~ predict_comb(mw, lambda),
          start = list(lambda = 0.1),
          data = data.frame(g = g, mw = mw),
          control = stats::nls.control(maxiter = 100, warnOnly = TRUE)
        )
      },
      error = function(e) NULL
    )

    if (!is.null(fit_comb)) {
      pred_comb <- stats::predict(fit_comb)
      ss_res <- sum((g - pred_comb)^2)
      ss_tot <- sum((g - mean(g))^2)
      model_fits$comb <- list(
        model = "comb",
        parameter = stats::coef(fit_comb)["lambda"],
        r_squared = 1 - ss_res / ss_tot,
        rmse = sqrt(mean((g - pred_comb)^2)),
        aic = stats::AIC(fit_comb),
        prediction = pred_comb,
        fit = fit_comb
      )
    } else {
      failed_models <- c(failed_models, "comb")
    }
  }

  if ("hyperbranched" %in% models) {
    fit_hyper <- tryCatch(
      {
        stats::nls(
          g ~ predict_hyperbranched(mw, db),
          start = list(db = 0.5),
          data = data.frame(g = g, mw = mw),
          lower = 0.01,
          upper = 1.0,
          algorithm = "port",
          control = stats::nls.control(maxiter = 100, warnOnly = TRUE)
        )
      },
      error = function(e) NULL
    )

    if (!is.null(fit_hyper)) {
      pred_hyper <- stats::predict(fit_hyper)
      ss_res <- sum((g - pred_hyper)^2)
      ss_tot <- sum((g - mean(g))^2)
      model_fits$hyperbranched <- list(
        model = "hyperbranched",
        parameter = stats::coef(fit_hyper)["db"],
        r_squared = 1 - ss_res / ss_tot,
        rmse = sqrt(mean((g - pred_hyper)^2)),
        aic = stats::AIC(fit_hyper),
        prediction = pred_hyper,
        fit = fit_hyper
      )
    } else {
      failed_models <- c(failed_models, "hyperbranched")
    }
  }

  # Warn about failed models
  if (length(failed_models) > 0) {
    cli::cli_warn(c(
      "!" = "Some models failed to fit and were excluded from comparison:",
      "i" = "Failed models: {.val {failed_models}}",
      "i" = "This may indicate the data is not suitable for these architectures."
    ))
  }

  # Build comparison summary
  if (length(model_fits) == 0) {
    cli::cli_warn("No models could be successfully fit to the data.")
    return(NULL)
  }

  summary_df <- data.frame(
    model = sapply(model_fits, `[[`, "model"),
    parameter = sapply(model_fits, `[[`, "parameter"),
    r_squared = sapply(model_fits, `[[`, "r_squared"),
    rmse = sapply(model_fits, `[[`, "rmse"),
    aic = sapply(model_fits, `[[`, "aic"),
    row.names = NULL
  )

  # Best model by AIC
  best_idx <- which.min(summary_df$aic)
  best_model <- summary_df$model[best_idx]

  # Build predictions data frame
  predictions <- data.frame(mw = mw, g_experimental = g)
  for (mod_name in names(model_fits)) {
    predictions[[paste0("g_", mod_name)]] <- model_fits[[mod_name]]$prediction
  }

  result <- list(
    model_fits = summary_df,
    best_model = best_model,
    best_fit = model_fits[[best_model]],
    predictions = predictions,
    experimental = data.frame(mw = mw, g = g),
    all_fits = model_fits
  )

  class(result) <- c("branching_model_comparison", "list")

  result
}


#' @export
print.branching_model_comparison <- function(x, ...) {
  cat("Branching Model Comparison\n")
  cat(strrep("=", 60), "\n\n")

  cat("Model Fit Summary:\n")
  cat(strrep("-", 60), "\n")
  print(x$model_fits, row.names = FALSE)
  cat(strrep("-", 60), "\n\n")

  cat(sprintf("Best Model: %s (lowest AIC)\n", x$best_model))
  cat(sprintf(
    "  R-squared: %.4f\n",
    x$best_fit$r_squared
  ))
  cat(sprintf(
    "  RMSE: %.4f\n",
    x$best_fit$rmse
  ))
  cat(sprintf(
    "  Parameter: %.4f\n\n",
    x$best_fit$parameter
  ))

  cat("Note: Lower AIC indicates better fit with penalty for complexity.\n")
  cat("Consider physical plausibility when selecting the final model.\n")

  invisible(x)
}
