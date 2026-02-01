# ==============================================================================
# step_sec_band_broadening
#
# Axial dispersion (band broadening) correction for SEC
# ==============================================================================

#' Band Broadening Correction for SEC
#'
#' `step_sec_band_broadening()` creates a *specification* of a recipe step that
#' corrects for axial dispersion (band broadening) in SEC chromatograms. This
#' improves the accuracy of molecular weight distribution measurements.
#'
#' @param recipe A recipe object.
#' @param measures Character vector of measure column names to process.
#'   If `NULL`, all measure columns will be processed.
#' @param method Correction method. One of:
#'   - `"tung"` (default): Tung's linear correction
#'   - `"emg"`: Exponentially Modified Gaussian deconvolution
#' @param sigma Spreading parameter (standard deviation of the instrumental
#'   broadening function) in the same units as the location axis (typically
#'   minutes or mL). If `NULL`, must provide `calibration_peak`.
#' @param calibration_peak A `measure_tbl` or data frame with `location` and
#'   `value` columns representing a narrow standard peak used to estimate sigma.
#' @param tau Exponential time constant for EMG method. If `NULL` with EMG

#'   method, estimated from `calibration_peak`.
#' @param iterations Number of iterations for iterative correction. Default is 1
#'   for Tung's method (single pass). Higher values may improve correction but
#'   can introduce instability.
#' @param damping Damping factor (0-1) to prevent over-correction and
#'   instability. Default is 0.5. Lower values are more conservative.
#' @param role Role for generated columns.
#' @param trained Logical indicating if the step has been trained.
#' @param skip Logical. Should the step be skipped when baking?
#' @param id Unique step identifier.
#'
#' @return An updated recipe with the new step added.
#'
#' @details
#' Band broadening in SEC occurs due to:
#' \itemize{
#'   \item Axial diffusion during elution
#'   \item Non-ideal column packing
#'   \item Extra-column volume (tubing, connections, detector cell)
#' }
#'
#' This causes:
#' \itemize{
#'   \item Artificially broadened peaks
#'   \item Underestimated Mn (number-average MW)
#'   \item Overestimated dispersity (Mw/Mn)
#' }
#'
#' **Tung's Method** (default):
#'
#' The observed chromatogram F(V) is related to the true distribution W(V) by:
#' \deqn{F(V) = \int W(V') G(V - V') dV'}
#'
#' where G is a Gaussian spreading function with standard deviation sigma.
#' Tung's linear correction approximates:
#' \deqn{W(V) \approx F(V) - \sigma^2 \frac{d^2 F(V)}{dV^2}}
#'
#' **EMG Method**:
#'
#' Models band broadening as convolution with an Exponentially Modified
#' Gaussian, which better handles asymmetric peak shapes caused by tailing.
#'
#' **Sigma Determination**:
#'
#' The spreading parameter sigma should be determined from a narrow molecular
#' weight standard (e.g., polystyrene with PDI < 1.05). Use [estimate_sigma()]
#' to calculate sigma from such a standard.
#'
#' @note
#' \itemize{
#'   \item Correction is applied to the signal, not to molecular weight values
#'   \item Large corrections (> 50% change in peak width) may indicate
#'     unreliable sigma or poor chromatographic conditions
#'   \item This step preserves the area under the curve (mass conservation)
#' }
#'
#' @references
#' Tung, L. H. (1966). Method of calculating molecular weight distribution
#' function from gel permeation chromatograms. Journal of Applied Polymer
#' Science, 10(3), 375-385.
#'
#' @family sec-chromatography
#' @seealso [estimate_sigma()] for determining the spreading parameter
#' @export
#'
#' @examples
#' \dontrun{
#' library(recipes)
#' library(measure)
#'
#' # Using a known sigma value
#' rec <- recipe(~., data = sec_data) |>
#'   step_measure_input_long(signal, location = vars(time), col_name = "ri") |>
#'   step_sec_baseline() |>
#'   step_sec_band_broadening(sigma = 0.05) |>
#'   prep()
#'
#' # Estimating sigma from a narrow standard
#' narrow_std <- estimate_sigma(narrow_standard_peak)
#' rec <- recipe(~., data = sec_data) |>
#'   step_measure_input_long(signal, location = vars(time), col_name = "ri") |>
#'   step_sec_baseline() |>
#'   step_sec_band_broadening(sigma = narrow_std$sigma) |>
#'   prep()
#' }
step_sec_band_broadening <- function(
	recipe,
	measures = NULL,
	method = c("tung", "emg"),
	sigma = NULL,
	calibration_peak = NULL,
	tau = NULL,
	iterations = 1,
	damping = 0.5,
	role = NA,
	trained = FALSE,

	skip = FALSE,
	id = recipes::rand_id("sec_band_broadening")
) {
	method <- match.arg(method)

	# Validate inputs

	if (is.null(sigma) && is.null(calibration_peak)) {
		cli::cli_abort(
			c(
				"Either {.arg sigma} or {.arg calibration_peak} must be provided.",
				"i" = "Use {.fn estimate_sigma} on a narrow standard to determine sigma."
			)
		)
	}

	if (!is.null(sigma)) {
		if (!is.numeric(sigma) || length(sigma) != 1 || sigma <= 0) {
			cli::cli_abort("{.arg sigma} must be a positive numeric value.")
		}
	}

	if (!is.numeric(iterations) || iterations < 1) {
		cli::cli_abort("{.arg iterations} must be a positive integer.")
	}

	if (!is.numeric(damping) || damping <= 0 || damping > 1) {
		cli::cli_abort("{.arg damping} must be between 0 and 1.")
	}

	if (!is.null(tau)) {
		if (!is.numeric(tau) || length(tau) != 1 || tau <= 0) {
			cli::cli_abort("{.arg tau} must be a positive numeric value.")
		}
	}

	recipes::add_step(
		recipe,
		step_sec_band_broadening_new(
			measures = measures,
			method = method,
			sigma = sigma,
			calibration_peak = calibration_peak,
			tau = tau,
			iterations = as.integer(iterations),
			damping = damping,
			estimated_sigma = NULL,
			estimated_tau = NULL,
			role = role,
			trained = trained,
			skip = skip,
			id = id
		)
	)
}

step_sec_band_broadening_new <- function(
	measures,
	method,
	sigma,
	calibration_peak,
	tau,
	iterations,
	damping,
	estimated_sigma,
	estimated_tau,
	role,
	trained,
	skip,
	id
) {
	recipes::step(
		subclass = "sec_band_broadening",
		measures = measures,
		method = method,
		sigma = sigma,
		calibration_peak = calibration_peak,
		tau = tau,
		iterations = iterations,
		damping = damping,
		estimated_sigma = estimated_sigma,
		estimated_tau = estimated_tau,
		role = role,
		trained = trained,
		skip = skip,
		id = id
	)
}

#' @export
prep.step_sec_band_broadening <- function(x, training, info = NULL, ...) {
	check_for_measure(training)

	# Resolve which columns to process
	if (is.null(x$measures)) {
		measure_cols <- find_measure_cols(training)
	} else {
		measure_cols <- x$measures
	}

	# Estimate sigma from calibration peak if not provided
	if (is.null(x$sigma) && !is.null(x$calibration_peak)) {
		sigma_result <- estimate_sigma(x$calibration_peak)
		estimated_sigma <- sigma_result$sigma
		estimated_tau <- sigma_result$tau

		cli::cli_inform(
			c(
				"i" = "Estimated spreading parameters from calibration peak:",
				"*" = "sigma = {round(estimated_sigma, 4)}",
				"*" = "tau = {round(estimated_tau, 4)}"
			)
		)
	} else {
		estimated_sigma <- x$sigma
		estimated_tau <- x$tau
	}

	# Warn if sigma seems too large
	# Check against first sample to estimate relative correction
	if (length(measure_cols) > 0) {
		first_col <- training[[measure_cols[1]]]
		if (length(first_col) > 0) {
			first_peak <- first_col[[1]]
			peak_width <- .estimate_peak_width(first_peak$location, first_peak$value)

			if (!is.na(peak_width) && estimated_sigma > 0.5 * peak_width) {
				cli::cli_warn(
					c(
						"Sigma ({round(estimated_sigma, 4)}) is > 50% of peak width ({round(peak_width, 4)}).",
						"i" = "This may indicate poor chromatographic conditions or incorrect sigma.",
						"i" = "Consider verifying sigma with a narrow standard."
					)
				)
			}
		}
	}

	step_sec_band_broadening_new(
		measures = measure_cols,
		method = x$method,
		sigma = x$sigma,
		calibration_peak = x$calibration_peak,
		tau = x$tau,
		iterations = x$iterations,
		damping = x$damping,
		estimated_sigma = estimated_sigma,
		estimated_tau = estimated_tau,
		role = x$role,
		trained = TRUE,
		skip = x$skip,
		id = x$id
	)
}

#' @export
bake.step_sec_band_broadening <- function(object, new_data, ...) {
	sigma <- object$estimated_sigma
	tau <- object$estimated_tau

	method <- object$method
	iterations <- object$iterations
	damping <- object$damping

	# Validate sigma is available and valid
	if (is.null(sigma) || is.na(sigma)) {
		cli::cli_abort(
			c(
				"Band broadening correction failed: sigma not available.",
				"i" = "This indicates a problem during prep(). Check calibration peak data."
			)
		)
	}

	if (!is.numeric(sigma) || sigma <= 0) {
		cli::cli_abort(
			c(
				"Invalid sigma value: {sigma}",
				"i" = "Sigma must be a positive number."
			)
		)
	}

	for (col in object$measures) {
		new_data[[col]] <- new_measure_list(
			purrr::map(new_data[[col]], function(m) {
				corrected <- .apply_band_broadening_correction(
					location = m$location,
					value = m$value,
					sigma = sigma,
					tau = tau,
					method = method,
					iterations = iterations,
					damping = damping
				)
				new_measure_tbl(location = m$location, value = corrected)
			})
		)
	}

	tibble::as_tibble(new_data)
}

#' @export
print.step_sec_band_broadening <- function(
	x,
	width = max(20, options()$width - 30),
	...
) {
	title <- paste0("SEC band broadening correction (", x$method, ")")

	if (x$trained) {
		sigma_val <- x$estimated_sigma %||% x$sigma
		cat(title, ", sigma = ", round(sigma_val, 4), sep = "")
	} else {
		cat(title)
	}
	cat("\n")
	invisible(x)
}

#' @rdname tidy.step_sec
#' @export
#' @keywords internal
tidy.step_sec_band_broadening <- function(x, ...) {
	tibble::tibble(
		method = x$method,
		sigma = x$estimated_sigma %||% x$sigma %||% NA_real_,
		tau = x$estimated_tau %||% x$tau %||% NA_real_,
		iterations = x$iterations,
		damping = x$damping,
		id = x$id
	)
}

#' @rdname required_pkgs.step_sec
#' @export
#' @keywords internal
required_pkgs.step_sec_band_broadening <- function(x, ...) {
	c("measure.sec", "measure")
}

# ==============================================================================
# Sigma estimation
# ==============================================================================

#' Estimate Spreading Parameter from Narrow Standard
#'
#' Estimates the instrumental spreading parameter (sigma) from a narrow
#' molecular weight standard peak. This sigma value can then be used in
#' `step_sec_band_broadening()` to correct for band broadening.
#'
#' @param peak A `measure_tbl` or data frame with `location` and `value`
#'   columns representing a chromatographic peak from a narrow MW standard.
#' @param method Method for sigma estimation:
#'   - `"gaussian"` (default): Fit a Gaussian and extract sigma
#'   - `"fwhm"`: Calculate from full width at half maximum (sigma = FWHM / 2.355)
#'   - `"moments"`: Calculate from second moment of the peak
#'
#' @return A list with components:
#'   \describe{
#'     \item{sigma}{The estimated spreading parameter (Gaussian std dev)}
#'     \item{tau}{Exponential tail parameter (for EMG model), NA if not applicable}
#'     \item{fwhm}{Full width at half maximum}
#'     \item{asymmetry}{Peak asymmetry factor (> 1 indicates tailing)}
#'     \item{method}{The method used for estimation}
#'   }
#'
#' @details
#' The spreading parameter sigma represents the standard deviation of the
#' instrumental broadening function, assumed to be approximately Gaussian.
#'
#' For best results, use a narrow polydispersity standard (PDI < 1.05) run
#' under the same conditions as your samples.
#'
#' **Method Details:**
#' \itemize{
#'   \item **gaussian**: Fits a Gaussian function to the peak using nonlinear
#'     least squares. Most accurate for symmetric peaks.
#'   \item **fwhm**: Uses the relationship sigma = FWHM / 2.355 for a Gaussian.
#'     Fast and robust but assumes symmetric peak.
#'   \item **moments**: Calculates sigma from the second central moment.
#'     Accounts for asymmetry but sensitive to baseline.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # From a measure_tbl
#' sigma_result <- estimate_sigma(narrow_standard_peak)
#' sigma_result$sigma
#'
#' # Use in band broadening correction
#' rec <- recipe(~., data = sec_data) |>
#'   step_sec_band_broadening(sigma = sigma_result$sigma)
#' }
estimate_sigma <- function(peak, method = c("gaussian", "fwhm", "moments")) {
	method <- match.arg(method)

	# Extract location and value

	if (inherits(peak, "measure_tbl") || inherits(peak, "data.frame")) {
		if (!all(c("location", "value") %in% names(peak))) {
			cli::cli_abort(
				"{.arg peak} must have {.val location} and {.val value} columns."
			)
		}
		x <- peak$location
		y <- peak$value
	} else {
		cli::cli_abort(
			"{.arg peak} must be a measure_tbl or data frame."
		)
	}

	# Remove NAs
	valid <- !is.na(x) & !is.na(y)
	x <- x[valid]
	y <- y[valid]

	if (length(x) < 5) {
		cli::cli_abort("Peak must have at least 5 valid data points.")
	}

	# Normalize to positive values with baseline near zero
	y_min <- min(y)
	y <- y - y_min

	# Find peak maximum
	max_idx <- which.max(y)
	x_max <- x[max_idx]
	y_max <- y[max_idx]

	# Calculate FWHM
	half_max <- y_max / 2
	left_idx <- max(1, max_idx - 1)
	right_idx <- min(length(y), max_idx + 1)

	# Find left crossing
	for (i in seq(max_idx, 1)) {
		if (y[i] < half_max) {
			left_idx <- i
			break
		}
	}
	# Interpolate for more accurate position
	if (left_idx < max_idx && left_idx > 1) {
		denom_left <- y[left_idx + 1] - y[left_idx]
		if (abs(denom_left) > .Machine$double.eps) {
			x_left <- x[left_idx] +
				(half_max - y[left_idx]) /
					denom_left *
					(x[left_idx + 1] - x[left_idx])
		} else {
			x_left <- x[left_idx]
		}
	} else {
		x_left <- x[left_idx]
	}

	# Find right crossing
	for (i in seq(max_idx, length(y))) {
		if (y[i] < half_max) {
			right_idx <- i
			break
		}
	}
	# Interpolate
	if (right_idx > max_idx && right_idx > 1) {
		denom_right <- y[right_idx] - y[right_idx - 1]
		if (abs(denom_right) > .Machine$double.eps) {
			x_right <- x[right_idx - 1] +
				(half_max - y[right_idx - 1]) /
					denom_right *
					(x[right_idx] - x[right_idx - 1])
		} else {
			x_right <- x[right_idx]
		}
	} else {
		x_right <- x[right_idx]
	}

	fwhm <- abs(x_right - x_left)

	# Calculate asymmetry (right width / left width at 10% height)
	tenth_max <- y_max * 0.1
	left_10 <- x[max_idx]
	right_10 <- x[max_idx]

	for (i in seq(max_idx, 1)) {
		if (y[i] < tenth_max) {
			left_10 <- x[i]
			break
		}
	}
	for (i in seq(max_idx, length(y))) {
		if (y[i] < tenth_max) {
			right_10 <- x[i]
			break
		}
	}

	asymmetry <- if (abs(x_max - left_10) > .Machine$double.eps) {
		abs(right_10 - x_max) / abs(x_max - left_10)
	} else {
		1.0
	}

	# Estimate sigma based on method
	if (method == "fwhm") {
		sigma <- fwhm / 2.355 # FWHM = 2 * sqrt(2 * ln(2)) * sigma
		tau <- NA_real_
	} else if (method == "moments") {
		# Second central moment
		total_area <- sum(y, na.rm = TRUE)
		if (total_area > 0) {
			mean_x <- sum(x * y, na.rm = TRUE) / total_area
			var_x <- sum((x - mean_x)^2 * y, na.rm = TRUE) / total_area
			sigma <- sqrt(var_x)
		} else {
			sigma <- fwhm / 2.355
		}
		tau <- NA_real_
	} else {
		# Gaussian fit
		# Starting values
		start_sigma <- fwhm / 2.355

		fit <- tryCatch(
			{
				stats::nls(
					y ~ amplitude * exp(-(x - mu)^2 / (2 * sigma^2)),
					data = data.frame(x = x, y = y),
					start = list(amplitude = y_max, mu = x_max, sigma = start_sigma),
					control = stats::nls.control(maxiter = 100, warnOnly = TRUE)
				)
			},
			error = function(e) NULL
		)

		if (!is.null(fit)) {
			sigma <- abs(stats::coef(fit)["sigma"])
		} else {
			# Fall back to FWHM method
			cli::cli_warn(
				"Gaussian fit failed; using FWHM method instead."
			)
			sigma <- fwhm / 2.355
		}
		tau <- NA_real_
	}

	# For EMG, estimate tau from asymmetry if peak is tailed
	if (asymmetry > 1.2) {
		# Approximate: tau ~ (asymmetry - 1) * sigma for mild tailing
		tau <- (asymmetry - 1) * sigma * 0.5
	}

	list(
		sigma = as.numeric(sigma),
		tau = as.numeric(tau),
		fwhm = fwhm,
		asymmetry = asymmetry,
		method = method
	)
}

# ==============================================================================
# Internal correction functions
# ==============================================================================

#' Apply band broadening correction to a single chromatogram
#' @noRd
.apply_band_broadening_correction <- function(
	location,
	value,
	sigma,
	tau,
	method,
	iterations,
	damping
) {
	n <- length(value)

	if (n < 10) {
		cli::cli_warn(
			"Chromatogram has fewer than 10 points; skipping band broadening correction."
		)
		return(value)
	}

	# Store original for area conservation check
	original_area <- sum(value, na.rm = TRUE)

	# Calculate spacing (assume uniform)
	dx <- mean(diff(location), na.rm = TRUE)

	if (is.na(dx) || !is.finite(dx) || dx <= 0) {
		cli::cli_warn(
			c(
				"Invalid location spacing (dx = {dx}).",
				"i" = "Skipping band broadening correction."
			)
		)
		return(value)
	}

	if (method == "tung") {
		corrected <- .tung_correction(value, sigma, dx, iterations, damping)
	} else if (method == "emg") {
		corrected <- .emg_correction(value, sigma, tau, dx, iterations, damping)
	} else {
		corrected <- value
	}

	# Preserve area (mass conservation)
	corrected_area <- sum(corrected, na.rm = TRUE)
	if (abs(corrected_area) > .Machine$double.eps) {
		corrected <- corrected * (original_area / corrected_area)
	}

	# Ensure non-negative (physical constraint for concentration)
	n_negative <- sum(corrected < 0, na.rm = TRUE)
	if (n_negative > 0) {
		pct_negative <- round(100 * n_negative / length(corrected), 1)
		if (pct_negative > 5) {
			cli::cli_warn(
				c(
					"Clipped {n_negative} negative values ({pct_negative}% of points).",
					"i" = "This may indicate over-correction. Consider reducing sigma or damping."
				)
			)
		}
	}
	corrected[corrected < 0] <- 0

	corrected
}

#' Tung's linear correction method
#' @noRd
.tung_correction <- function(value, sigma, dx, iterations, damping) {
	corrected <- value

	for (iter in seq_len(iterations)) {
		# Calculate second derivative using central differences
		n <- length(corrected)
		d2y <- numeric(n)

		# Central difference for interior points
		for (i in seq(2, n - 1)) {
			d2y[i] <- (corrected[i + 1] - 2 * corrected[i] + corrected[i - 1]) / dx^2
		}

		# Forward/backward difference for endpoints
		if (n >= 3) {
			d2y[1] <- (corrected[3] - 2 * corrected[2] + corrected[1]) / dx^2
			d2y[n] <- (corrected[n] - 2 * corrected[n - 1] + corrected[n - 2]) / dx^2
		}

		# Tung's correction: W(V) = F(V) - sigma^2 * d2F/dV2
		correction <- sigma^2 * d2y

		# Apply with damping
		corrected <- corrected - damping * correction
	}

	corrected
}

#' EMG (Exponentially Modified Gaussian) correction
#' @noRd
.emg_correction <- function(value, sigma, tau, dx, iterations, damping) {
	# If tau is not available or very small, fall back to Tung's method
	if (is.null(tau) || is.na(tau) || tau < sigma * 0.1) {
		cli::cli_warn(
			c(
				"EMG method requested but tau is unavailable or too small.",
				"i" = "Falling back to Tung's method for band broadening correction."
			)
		)
		return(.tung_correction(value, sigma, dx, iterations, damping))
	}

	n <- length(value)

	# Build EMG kernel for deconvolution
	# EMG: convolution of Gaussian with exponential decay
	kernel_width <- ceiling(4 * (sigma + tau) / dx)
	kernel_x <- seq(-kernel_width, kernel_width) * dx

	# Gaussian part
	gaussian <- exp(-kernel_x^2 / (2 * sigma^2)) / (sigma * sqrt(2 * pi))

	# Exponential part (causal - only positive x contributes to tailing)
	exponential <- ifelse(kernel_x >= 0, exp(-kernel_x / tau) / tau, 0)

	# Convolve to get EMG kernel
	emg_kernel <- stats::convolve(gaussian, exponential, type = "open")
	emg_kernel <- emg_kernel / sum(emg_kernel) # Normalize

	# Simple iterative deconvolution (Richardson-Lucy style)
	corrected <- value

	for (iter in seq_len(iterations)) {
		# Reconvolve the estimate
		reconvolved <- stats::convolve(
			corrected,
			emg_kernel,
			type = "open"
		)

		# Trim to original length
		start_idx <- floor(length(reconvolved) / 2) - floor(n / 2) + 1
		end_idx <- start_idx + n - 1

		# Bounds check to prevent index errors
		start_idx <- max(1, start_idx)
		end_idx <- min(length(reconvolved), end_idx)

		if (end_idx - start_idx + 1 != n) {
			# If we can't extract the right length, pad with edge values
			extracted <- reconvolved[start_idx:end_idx]
			if (length(extracted) < n) {
				pad_length <- n - length(extracted)
				extracted <- c(
					rep(extracted[1], floor(pad_length / 2)),
					extracted,
					rep(extracted[length(extracted)], ceiling(pad_length / 2))
				)
			}
			reconvolved <- extracted[seq_len(n)]
		} else {
			reconvolved <- reconvolved[start_idx:end_idx]
		}

		# Update estimate with damping
		ratio <- value / pmax(reconvolved, .Machine$double.eps)
		ratio[!is.finite(ratio)] <- 1

		corrected <- corrected * (1 + damping * (ratio - 1))
	}

	corrected
}

#' Estimate peak width from chromatogram
#' @noRd
.estimate_peak_width <- function(location, value) {
	if (length(value) < 5 || all(is.na(value))) {
		return(NA_real_)
	}

	value <- value - min(value, na.rm = TRUE)
	max_val <- max(value, na.rm = TRUE)

	if (max_val <= 0) {
		return(NA_real_)
	}

	half_max <- max_val / 2
	max_idx <- which.max(value)

	# Find left and right crossings
	left_x <- location[1]
	right_x <- location[length(location)]

	for (i in seq(max_idx, 1)) {
		if (value[i] < half_max) {
			left_x <- location[i]
			break
		}
	}

	for (i in seq(max_idx, length(value))) {
		if (value[i] < half_max) {
			right_x <- location[i]
			break
		}
	}

	abs(right_x - left_x)
}
