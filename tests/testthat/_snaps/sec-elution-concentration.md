# clipping warnings identify raw axis units without changing table shape

    Code
      for (warning in warnings) {
        cat(warning$elution_units, ": ", sep = "")
        cat(warning$diagnostics$positive_area, "/", warning$diagnostics$negative_area,
        "\n", sep = "")
        cat(conditionMessage(warning), "\n", sep = "")
      }
    Output
      min: 3/1
      Negative detector values were excluded during mass normalization.
      i 1 of 3 signal points (33.3%) were set to zero.
      i Column "ri": maximum negative share of absolute sampled signal area was 25% in row 1.
      s: 180/60
      Negative detector values were excluded during mass normalization.
      i 1 of 3 signal points (33.3%) were set to zero.
      i Column "ri": maximum negative share of absolute sampled signal area was 25% in row 1.
      mL: 1.5/0.5
      Negative detector values were excluded during mass normalization.
      i 1 of 3 signal points (33.3%) were set to zero.
      i Column "ri": maximum negative share of absolute sampled signal area was 25% in row 1.
      uL: 1500/500
      Negative detector values were excluded during mass normalization.
      i 1 of 3 signal points (33.3%) were set to zero.
      i Column "ri": maximum negative share of absolute sampled signal area was 25% in row 1.
