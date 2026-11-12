# mass normalization reports sampled positive and negative areas

    Code
      cat(conditionMessage(warning), "\n", sep = "")
    Output
      Negative detector values were excluded during mass normalization.
      i 2 of 9 signal points (22.2%) were set to zero.
      i Column "ri": maximum negative share of absolute sampled signal area was 76.9% in row 2.

# normalization failures are reported before clipping warnings

    Code
      for (name in names(errors)) {
        cat(name, ":\n", sep = "")
        cat(errors[[name]], "\n", sep = "")
      }
    Output
      all_negative:
      i In index: 1.
      Caused by error in `.sec_normalize_concentration_measure()`:
      ! Cannot normalize a concentration signal with nonpositive integrated area.
      all_zero:
      i In index: 1.
      Caused by error in `.sec_normalize_concentration_measure()`:
      ! Cannot normalize a concentration signal with nonpositive integrated area.
      late_row_all_negative:
      i In index: 2.
      Caused by error in `.sec_normalize_concentration_measure()`:
      ! Cannot normalize a concentration signal with nonpositive integrated area.
      positive_area_overflow:
      i In index: 1.
      Caused by error in `.sec_normalize_concentration_measure()`:
      ! Cannot normalize a concentration signal with non-finite positive integrated area.
      negative_area_overflow:
      i In index: 1.
      Caused by error in `.sec_normalize_concentration_measure()`:
      ! Concentration normalization requires a finite, nonnegative negative-signal area.
      axis_width_overflow:
      i In index: 1.
      Caused by error in `.sec_normalize_concentration_measure()`:
      ! Concentration normalization requires finite, positive axis widths.
      denominator_overflow:
      i In index: 1.
      Caused by error in `.sec_normalize_concentration_measure()`:
      ! Concentration normalization requires a finite, positive integrated area-flow denominator.
      denominator_underflow:
      i In index: 1.
      Caused by error in `.sec_normalize_concentration_measure()`:
      ! Concentration normalization requires a finite, positive integrated area-flow denominator.
      factor_underflow:
      i In index: 1.
      Caused by error in `.sec_normalize_concentration_measure()`:
      ! Concentration normalization requires a finite, positive normalization factor.
      factor_overflow:
      i In index: 1.
      Caused by error in `.sec_normalize_concentration_measure()`:
      ! Concentration normalization requires a finite, positive normalization factor.
      output_overflow:
      i In index: 1.
      Caused by error in `.sec_normalize_concentration_measure()`:
      ! Concentration normalization produced non-finite output values.
      output_underflow:
      i In index: 1.
      Caused by error in `.sec_normalize_concentration_measure()`:
      ! Concentration normalization underflowed all positive signal values to zero.

