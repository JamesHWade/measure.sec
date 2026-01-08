# ==============================================================================
# Tests for standards-helpers.R
# ==============================================================================

# -- get_recognized_standards tests --------------------------------------------

test_that("get_recognized_standards returns a named list", {
  standards <- get_recognized_standards()

  expect_type(standards, "list")
  expect_true(length(standards) > 0)
  expect_true(all(nchar(names(standards)) > 0))
})

test_that("get_recognized_standards contains expected standard types", {
  standards <- get_recognized_standards()
  names <- names(standards)

  # Check for PS standards
  expect_true(any(grepl("^PS", names)))

  # Check for PMMA standards
  expect_true(any(grepl("^PMMA", names)))

  # Check for generic standards
  expect_true(any(grepl("^Standard", names)))
})

# -- get_standard_names tests --------------------------------------------------

test_that("get_standard_names returns character vector", {
  names <- get_standard_names()

  expect_type(names, "character")
  expect_true(length(names) > 0)
})

test_that("get_standard_names matches get_recognized_standards", {
  standards <- get_recognized_standards()
  names <- get_standard_names()

  expect_equal(names, names(standards))
})

# -- get_standard_type tests ---------------------------------------------------

test_that("get_standard_type correctly classifies PS standards", {
  expect_equal(get_standard_type("PS A"), "PS")
  expect_equal(get_standard_type("PS B"), "PS")
  expect_equal(get_standard_type("PS 1683"), "PS")
})

test_that("get_standard_type correctly classifies PMMA standards", {
  expect_equal(get_standard_type("PMMA A"), "PMMA")
  expect_equal(get_standard_type("PMMA D"), "PMMA")
})

test_that("get_standard_type correctly classifies PEG standards", {
  expect_equal(get_standard_type("PEG A"), "PEG")
  expect_equal(get_standard_type("PEG C"), "PEG")
})

test_that("get_standard_type correctly classifies Pullulan standards", {
  expect_equal(get_standard_type("Pullulan A"), "Pullulan")
  expect_equal(get_standard_type("Pullulan B"), "Pullulan")
})

test_that("get_standard_type correctly classifies Generic standards", {
  expect_equal(get_standard_type("Standard 1"), "Generic")
  expect_equal(get_standard_type("Standard 4"), "Generic")
})

test_that("get_standard_type returns Other for unknown types", {
  expect_equal(get_standard_type("Unknown"), "Other")
  expect_equal(get_standard_type("My Custom Standard"), "Other")
})

# -- normalize_standard_name tests ---------------------------------------------

test_that("normalize_standard_name converts to lowercase", {
  expect_equal(normalize_standard_name("PS A"), "ps a")
  expect_equal(normalize_standard_name("PMMA B"), "pmma b")
})

test_that("normalize_standard_name replaces separators with spaces", {
  expect_equal(normalize_standard_name("PS-A"), "ps a")
  expect_equal(normalize_standard_name("PS_A"), "ps a")
  expect_equal(normalize_standard_name("PS_A_Test"), "ps a test")
})

test_that("normalize_standard_name converts std to standard", {
  expect_equal(normalize_standard_name("Std A"), "standard a")
  expect_equal(normalize_standard_name("StdA"), "standard a")
  expect_equal(normalize_standard_name("Std 1"), "standard 1")
  expect_equal(normalize_standard_name("Std1"), "standard 1")
  expect_equal(normalize_standard_name("std B"), "standard b")
})

test_that("normalize_standard_name removes date prefixes", {
  # YYYYMMDD format
  result <- normalize_standard_name("20231215_PS_A")
  expect_true(grepl("ps a", result))

  # YYYY-MM-DD format (converted to YYYY MM DD with space)
  result <- normalize_standard_name("2023-12-15 PS A")
  expect_true(grepl("ps a", result))
})

test_that("normalize_standard_name removes vendor identifiers", {
  result <- normalize_standard_name("Phillips PS A")
  expect_true(grepl("ps a", result))

  result <- normalize_standard_name("DOW PMMA B")
  expect_true(grepl("pmma b", result))
})

test_that("normalize_standard_name handles empty input", {
  expect_equal(normalize_standard_name(character(0)), character(0))
})

test_that("normalize_standard_name handles vector input", {
  input <- c("PS A", "Std B", "PMMA-C")
  result <- normalize_standard_name(input)

  expect_length(result, 3)
  expect_true(grepl("ps a", result[1]))
  expect_true(grepl("standard b", result[2]))
  expect_true(grepl("pmma c", result[3]))
})

# -- match_standards tests -----------------------------------------------------

test_that("match_standards returns tibble with expected columns", {
  samples <- c("PS A", "Unknown", "PMMA B")
  result <- match_standards(samples, return_all = TRUE)

  expect_s3_class(result, "tbl_df")
  expect_true("sample_name" %in% names(result))
  expect_true("normalized_name" %in% names(result))
  expect_true("matched_standard" %in% names(result))
  expect_true("standard_type" %in% names(result))
})

test_that("match_standards correctly matches PS standards", {
  samples <- c("PS A", "PS B", "PS-C", "PS_D", "PS 1683")
  result <- match_standards(samples)

  expect_equal(nrow(result), 5)
  expect_equal(
    result$matched_standard,
    c("PS A", "PS B", "PS C", "PS D", "PS 1683")
  )
})

test_that("match_standards correctly matches PMMA standards", {
  samples <- c("PMMA A", "PMMA-B", "PMMA_C")
  result <- match_standards(samples)

  expect_equal(nrow(result), 3)
  expect_equal(result$matched_standard, c("PMMA A", "PMMA B", "PMMA C"))
})

test_that("match_standards correctly matches generic standards", {
  samples <- c("Standard 1", "Std 2", "Std3", "StdD")
  result <- match_standards(samples)

  # Standard 1, Standard 2, Standard 3 should match
  # StdD -> standard d, which doesn't match Standard 1-4
  expect_true("Standard 1" %in% result$matched_standard)
  expect_true("Standard 2" %in% result$matched_standard)
  expect_true("Standard 3" %in% result$matched_standard)
})

test_that("match_standards handles prefixed names", {
  samples <- c("20231215_PS_A", "Phillips PMMA B", "DOW Std 1")
  result <- match_standards(samples)

  expect_equal(nrow(result), 3)
  expect_equal(result$matched_standard, c("PS A", "PMMA B", "Standard 1"))
})

test_that("match_standards excludes non-matches by default", {
  samples <- c("PS A", "Unknown Sample", "PMMA B", "Random Name")
  result <- match_standards(samples, return_all = FALSE)

  expect_equal(nrow(result), 2)
  expect_false("Unknown Sample" %in% result$sample_name)
})

test_that("match_standards includes all with return_all = TRUE", {
  samples <- c("PS A", "Unknown Sample", "PMMA B")
  result <- match_standards(samples, return_all = TRUE)

  expect_equal(nrow(result), 3)
  expect_true(is.na(result$matched_standard[
    result$sample_name == "Unknown Sample"
  ]))
})

test_that("match_standards includes standard_type", {
  samples <- c("PS A", "PMMA B", "Standard 1")
  result <- match_standards(samples)

  expect_equal(result$standard_type, c("PS", "PMMA", "Generic"))
})

# -- is_standard tests ---------------------------------------------------------

test_that("is_standard returns logical vector", {
  samples <- c("PS A", "Unknown", "PMMA B")
  result <- is_standard(samples)

  expect_type(result, "logical")
  expect_length(result, 3)
})

test_that("is_standard correctly identifies standards", {
  samples <- c("PS A", "Unknown Sample", "PMMA-B", "Random")
  result <- is_standard(samples)

  expect_equal(result, c(TRUE, FALSE, TRUE, FALSE))
})

# -- auto_detect_standards tests -----------------------------------------------

test_that("auto_detect_standards finds sample_name column", {
  data <- data.frame(
    sample_name = c("PS A", "Unknown", "PMMA B"),
    signal = c(100, 200, 150)
  )

  result <- auto_detect_standards(data, add_columns = TRUE)

  expect_true("matched_standard" %in% names(result))
  expect_true("standard_type" %in% names(result))
  expect_equal(result$matched_standard, c("PS A", NA, "PMMA B"))
})

test_that("auto_detect_standards uses specified name_col", {
  data <- data.frame(
    my_names = c("PS A", "PMMA B"),
    signal = c(100, 150)
  )

  result <- auto_detect_standards(
    data,
    name_col = "my_names",
    add_columns = TRUE
  )

  expect_equal(result$matched_standard, c("PS A", "PMMA B"))
})

test_that("auto_detect_standards errors on missing column", {
  data <- data.frame(
    signal = c(100, 200)
  )

  expect_error(
    auto_detect_standards(data),
    "auto-detect"
  )
})

test_that("auto_detect_standards with add_columns = FALSE returns match tibble", {
  data <- data.frame(
    sample_name = c("PS A", "PMMA B"),
    signal = c(100, 150)
  )

  result <- auto_detect_standards(data, add_columns = FALSE)

  expect_s3_class(result, "tbl_df")
  expect_true("sample_name" %in% names(result))
  expect_true("matched_standard" %in% names(result))
  expect_false("signal" %in% names(result))
})

# -- Edge cases ----------------------------------------------------------------

test_that("normalize_standard_name handles NA values gracefully", {
  # This might produce unexpected results, but shouldn't error
  result <- normalize_standard_name(c("PS A", NA, "PMMA B"))
  expect_length(result, 3)
})

test_that("match_standards handles empty input", {
  result <- match_standards(character(0))
  expect_equal(nrow(result), 0)
})

test_that("match_standards handles single input", {
  result <- match_standards("PS A")
  expect_equal(nrow(result), 1)
  expect_equal(result$matched_standard, "PS A")
})
