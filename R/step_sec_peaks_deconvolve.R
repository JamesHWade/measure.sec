# ==============================================================================
# step_sec_peaks_deconvolve
#
# SEC-optimized wrapper for peak deconvolution
# ==============================================================================

#' SEC/GPC Peak Deconvolution
#'
#' `step_sec_peaks_deconvolve()` creates a *specification* of a recipe step
#' that resolves overlapping peaks in SEC/GPC chromatograms using curve fitting.
#' This is a thin wrapper around [measure::step_measure_peaks_deconvolve()]
#' with SEC-optimized defaults.
#'
#' @param recipe A recipe object. The step will be added to the sequence of
#'   operations for this recipe.
#' @param model Peak model to use:
#'   - `"emg"` (default): Exponentially Modified Gaussian, recommended for SEC
#'     chromatograms with peak tailing
#'   - `"gaussian"`: Symmetric Gaussian, for well-behaved symmetric peaks
#'   - `"bigaussian"`: Bi-Gaussian for flexible asymmetry
#'
#'   EMG is the default because chromatographic peaks typically exhibit tailing
#'   due to mass transfer kinetics and column effects.
#' @param optimizer Optimization method:
#'   - `"auto"` (default): Selects based on problem complexity
#'   - `"lbfgsb"`: L-BFGS-B (fast, local optimization)
#'   - `"multistart"`: Multiple starts for robustness
#'   - `"nelder_mead"`: Derivative-free simplex method
#' @param max_iter Maximum iterations for optimization. Default is `500`.
#' @param quality_threshold Minimum R-squared to accept fit. Default is `0.9`
#'   (stricter than base measure default of 0.8 for analytical quality).
#' @param smart_init Logical. Use smart initialization based on peak properties.
#'   Default is `TRUE`. Highly recommended for SEC data.
#' @param constrain_positions Logical. Enforce that peak centers maintain their
#'   relative ordering during optimization. Default is `TRUE`.
#' @param peaks_col Name of the peaks column. Default is `".peaks"`.
#' @param measures_col Name of the measures column containing the chromatogram.
#'   Default is `".measures"` but often `"ri"`, `"uv"`, etc. in SEC workflows.
#' @param role Not used by this step.
#' @param trained A logical to indicate if the step has been trained.
#' @param skip A logical. Should the step be skipped when baking?
#' @param id A character string that is unique to this step.
#'
#' @return An updated version of `recipe` with the new step added. The `.peaks`
#'   column will be updated with:
#'   - Refined peak parameters from curve fitting
#'   - `fit_r_squared`: R-squared of the overall fit
#'   - `area`: Integrated area under the fitted curve (analytical integration)
#'
#' @details
#' Peak deconvolution is essential for SEC/GPC analysis when peaks overlap,
#' which is common for:
#' - Bimodal molecular weight distributions
#' - Aggregate/monomer/fragment species in protein SEC
#' - Copolymer component separation
#'
#' **EMG Model (Default):**
#'
#' The Exponentially Modified Gaussian (EMG) is ideal for SEC because it
#' accounts for the characteristic tailing observed in chromatographic peaks:
#'
#' \deqn{EMG(x) = h \cdot \exp(\sigma^2/(2\tau^2) + (c-x)/\tau) \cdot
#'   \Phi((x-c)/\sigma - \sigma/\tau)}
#'

#' where \eqn{h} is height, \eqn{c} is center, \eqn{\sigma} is Gaussian width,
#' and \eqn{\tau} is the exponential decay parameter (tailing).
#'
#' **Quality Assessment:**
#'
#' The `quality_threshold` parameter sets the minimum acceptable R-squared
#' for fits. The default of 0.9 is stricter than the base measure default,
#' appropriate for quantitative analytical work.
#'
#' # Tidying
#'
#' When you [`tidy()`][recipes::tidy.recipe()] this step, a tibble with columns
#' `model`, `optimizer`, `quality_threshold`, and `id` is returned.
#'
#' @seealso [step_sec_peaks_detect()] for SEC-optimized peak detection,
#'   [measure::step_measure_peaks_deconvolve()] for the underlying implementation.
#'
#' @family sec-chromatography
#' @export
#'
#' @examples
#' \dontrun{
#' library(recipes)
#' library(measure)
#'
#' # SEC peak deconvolution with EMG model (default)
#' rec <- recipe(~., data = sec_data) |>
#'   step_measure_input_long(ri_signal, location = vars(elution_volume)) |>
#'   step_sec_peaks_detect() |>
#'   step_sec_peaks_deconvolve() |>
#'   prep()
#'
#' # Use Gaussian model for symmetric peaks
#' rec <- recipe(~., data = sec_data) |>
#'   step_measure_input_long(ri_signal, location = vars(elution_volume)) |>
#'   step_sec_peaks_detect() |>
#'   step_sec_peaks_deconvolve(model = "gaussian") |>
#'   prep()
#'
#' # Use multistart optimizer for complex multi-peak scenarios
#' rec <- recipe(~., data = sec_data) |>
#'   step_measure_input_long(ri_signal, location = vars(elution_volume)) |>
#'   step_sec_peaks_detect() |>
#'   step_sec_peaks_deconvolve(optimizer = "multistart") |>
#'   prep()
#' }
step_sec_peaks_deconvolve <- function(
	recipe,
	model = "emg",
	optimizer = "auto",
	max_iter = 500L,
	quality_threshold = 0.9,
	smart_init = TRUE,
	constrain_positions = TRUE,
	peaks_col = ".peaks",
	measures_col = ".measures",
	role = NA,
	trained = FALSE,
	skip = FALSE,
	id = recipes::rand_id("sec_peaks_deconvolve")
) {
	# Validate model - SEC supports gaussian, emg, bigaussian
	valid_models <- c("gaussian", "emg", "bigaussian")
	if (is.character(model)) {
		model <- rlang::arg_match(model, valid_models)
	}

	# Validate optimizer
	valid_optimizers <- c("auto", "lbfgsb", "multistart", "nelder_mead")
	optimizer <- rlang::arg_match(optimizer, valid_optimizers)

	# Validate other parameters
	if (!is.numeric(max_iter) || length(max_iter) != 1 || max_iter < 1) {
		cli::cli_abort("{.arg max_iter} must be a positive integer.")
	}
	if (
		!is.numeric(quality_threshold) ||
			length(quality_threshold) != 1 ||
			quality_threshold < 0 ||
			quality_threshold > 1
	) {
		cli::cli_abort("{.arg quality_threshold} must be between 0 and 1.")
	}
	if (!is.logical(smart_init) || length(smart_init) != 1) {
		cli::cli_abort("{.arg smart_init} must be TRUE or FALSE.")
	}
	if (!is.logical(constrain_positions) || length(constrain_positions) != 1) {
		cli::cli_abort("{.arg constrain_positions} must be TRUE or FALSE.")
	}

	recipes::add_step(
		recipe,
		step_sec_peaks_deconvolve_new(
			model = model,
			optimizer = optimizer,
			max_iter = as.integer(max_iter),
			quality_threshold = quality_threshold,
			smart_init = smart_init,
			constrain_positions = constrain_positions,
			peaks_col = peaks_col,
			measures_col = measures_col,
			role = role,
			trained = trained,
			skip = skip,
			id = id
		)
	)
}

step_sec_peaks_deconvolve_new <- function(
	model,
	optimizer,
	max_iter,
	quality_threshold,
	smart_init,
	constrain_positions,
	peaks_col,
	measures_col,
	role,
	trained,
	skip,
	id
) {
	recipes::step(
		subclass = "sec_peaks_deconvolve",
		model = model,
		optimizer = optimizer,
		max_iter = max_iter,
		quality_threshold = quality_threshold,
		smart_init = smart_init,
		constrain_positions = constrain_positions,
		peaks_col = peaks_col,
		measures_col = measures_col,
		role = role,
		trained = trained,
		skip = skip,
		id = id
	)
}

#' @export
prep.step_sec_peaks_deconvolve <- function(x, training, info = NULL, ...) {
	# Verify peaks column exists
	if (!x$peaks_col %in% names(training)) {
		cli::cli_abort(
			c(
				"Column {.val {x$peaks_col}} not found.",
				"i" = "Run {.fn step_sec_peaks_detect} first to detect peaks."
			)
		)
	}

	# Check for measures column - allow common SEC naming patterns
	measures_col <- x$measures_col
	if (!measures_col %in% names(training)) {
		# Check for common SEC measure column names
		measure_cols <- find_measure_cols(training)
		if (length(measure_cols) > 0) {
			measures_col <- measure_cols[1]
			cli::cli_inform(
				c(
					"i" = "Column {.val {x$measures_col}} not found.",
					"i" = "Using {.val {measures_col}} instead."
				)
			)
		} else {
			cli::cli_abort("No measure columns found in data.")
		}
	}

	step_sec_peaks_deconvolve_new(
		model = x$model,
		optimizer = x$optimizer,
		max_iter = x$max_iter,
		quality_threshold = x$quality_threshold,
		smart_init = x$smart_init,
		constrain_positions = x$constrain_positions,
		peaks_col = x$peaks_col,
		measures_col = measures_col,
		role = x$role,
		trained = TRUE,
		skip = x$skip,
		id = x$id
	)
}

#' @export
bake.step_sec_peaks_deconvolve <- function(object, new_data, ...) {
	# Delegate to measure's deconvolution step
	# We create a temporary recipe step and bake it
	model <- object$model
	optimizer <- object$optimizer
	max_iter <- object$max_iter
	quality_threshold <- object$quality_threshold
	smart_init <- object$smart_init
	constrain_positions <- object$constrain_positions
	peaks_col <- object$peaks_col
	measures_col <- object$measures_col

	# Call measure's bake directly with our parameters
	# Create a temporary step object with measure's class
	temp_step <- list(
		model = model,
		optimizer = optimizer,
		max_iter = max_iter,
		tol = 1e-6, # Use measure's default
		n_starts = 5L, # Use measure's default for multistart
		quality_threshold = quality_threshold,
		constrain_positions = constrain_positions,
		store_components = FALSE,
		smart_init = smart_init,
		peaks_col = peaks_col,
		measures_col = measures_col,
		trained = TRUE,
		skip = FALSE,
		id = object$id
	)
	class(temp_step) <- c("step_measure_peaks_deconvolve", "step")

	# Use measure's bake method
	measure:::bake.step_measure_peaks_deconvolve(temp_step, new_data, ...)
}

#' @export
print.step_sec_peaks_deconvolve <- function(
	x,
	width = max(20, options()$width - 30),
	...
) {
	model_name <- if (is.character(x$model)) x$model else x$model$name
	cat(
		glue::glue(
			"SEC peak deconvolution ({model_name}, {x$optimizer})"
		)
	)
	cat("\n")
	invisible(x)
}

#' @rdname tidy.step_sec
#' @export
#' @keywords internal
tidy.step_sec_peaks_deconvolve <- function(x, ...) {
	model_name <- if (is.character(x$model)) x$model else "custom"

	tibble::tibble(
		model = model_name,
		optimizer = x$optimizer,
		max_iter = x$max_iter,
		quality_threshold = x$quality_threshold,
		smart_init = x$smart_init,
		id = x$id
	)
}

#' @rdname required_pkgs.step_sec
#' @export
#' @keywords internal
required_pkgs.step_sec_peaks_deconvolve <- function(x, ...) {
	c("measure.sec", "measure")
}
