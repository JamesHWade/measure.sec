# positive MW support with a calibration gap errors

    Code
      .sec_mw_slices(signal = list(location = 1:4, value = rep(1, 4)), mw_curve = list(
        location = 1:4, value = c(100, NA_real_, 300, 400)), mw_scale = "absolute")
    Condition
      Error in `.sec_mw_slices()`:
      ! Positive signal includes slices without a valid molecular-weight calibration; adjust the integration range or calibration.

# off-grid support outside valid calibration endpoints errors

    Code
      .sec_mw_slices(signal = list(location = 0:2, value = rep(1, 3)), mw_curve = list(
        location = 0:2, value = c(NA_real_, 1, 2)), mw_scale = "absolute",
      integration_range = c(0.5, 2))
    Condition
      Error in `.sec_mw_slices()`:
      ! Positive signal includes slices without a valid molecular-weight calibration; adjust the integration range or calibration.

---

    Code
      .sec_mw_slices(signal = list(location = 0:2, value = rep(1, 3)), mw_curve = list(
        location = 0:2, value = c(1, 2, NA_real_)), mw_scale = "absolute",
      integration_range = c(0, 1.5))
    Condition
      Error in `.sec_mw_slices()`:
      ! Positive signal includes slices without a valid molecular-weight calibration; adjust the integration range or calibration.

# row integration metadata uses overrides and validates metadata

    Code
      .sec_row_integration_metadata(malformed, 1)
    Condition
      Error in `.sec_validate_integration_range()`:
      ! `integration_range` must contain two increasing finite values.

---

    Code
      .sec_row_integration_metadata(malformed, 1)
    Condition
      Error in `.sec_validate_excluded_regions()`:
      ! `.excluded_regions` contains malformed region metadata.

# paired MW axes validate before support selection

    Code
      .sec_mw_slices(signal = list(location = 1:3, value = rep(1, 3)), mw_curve = list(
        location = c(1, 2, 4), value = c(100, 200, 400)), mw_scale = "absolute",
      integration_range = c(10, 11))
    Condition
      Error in `.sec_mw_slices()`:
      ! Molecular-weight and concentration curve locations must align.

