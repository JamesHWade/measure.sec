# baseline-only exclusions are ignored and malformed metadata fails

    Code
      recipes::bake(prepped, new_data = malformed)
    Condition
      Error in `purrr::map()`:
      i In index: 1.
      Caused by error in `.sec_validate_excluded_regions()`:
      ! `.excluded_regions` has invalid column types for start, end, or purpose.

# missing calibration is tolerated only inside excluded support

    Code
      recipes::bake(recipes::prep(step_sec_mw_distribution(recipes::recipe(
        signal_data), measures = "ri", mw_col = "mw", mw_scale = "absolute", type = "both",
      n_points = 10, normalize = FALSE)), new_data = NULL)
    Condition
      Error in `step_sec_mw_distribution()`:
      Caused by error in `purrr::map()`:
      i In index: 1.
      Caused by error in `.sec_mw_slices()`:
      ! Positive signal includes slices without a valid molecular-weight calibration; adjust the integration range or calibration.
