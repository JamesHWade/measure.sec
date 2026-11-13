# measure.sec 0.0.0.9000

* `step_sec_detector_delay()` accepts explicit minute, second, mL, and uL
  axes. Positive delays move late target signals earlier with both correction
  methods, fixing the reversed interpolation direction. Index correction
  requires a uniform grid; interpolation preserves missing gaps and rejects
  unrepresentable results. Current reference/target grids must match (#82).

* `step_sec_concentration()` accepts explicit minute, second, mL, and uL
  elution axes without changing coordinates. Mass normalization converts the
  sampled area to mL; invalid declarations and malformed curves fail before
  normalization. Named concentration units now convert consistently across
  detector calculations, fixing a factor-of-1,000 error (#82).

* `step_sec_baseline()` retains the current input column's recorded
  concentration unit after correction, including unchanged-signal fallbacks.
  Untagged inputs stay untagged and invalid recorded units are rejected;
  no other cached column metadata is copied (#60).

* The new detector qualification guide uses fixed synthetic traces to check
  signed RI conversion, retained raw outputs, mg/mL versus g/mL mass balance,
  and the difference between negative-point counts and sampled-area shares.
  It documents minute-axis assumptions and physical qualification limits (#60).

* `measure.sec` now requires R 4.5.0 or later (#70).

* `step_sec_mw_distribution()` now respects row integration windows and
  exclusions, conserves mass within each retained interval, and keeps the CDF
  flat across excluded support. Differential curves encode discontinuities
  with repeated boundary locations and zero density across gaps; output grids
  include support boundaries when needed. A new guide demonstrates the full
  calibrated selection workflow (#57).

* `step_sec_mw_averages()` and `step_sec_mw_fractions()` now use each baked
  sample's integration window and integration-purpose exclusions. An explicit
  `integration_range` overrides the stored window. Clipped boundaries are
  interpolated in the declared MW scale, and trapezoidal weights keep
  disconnected support separate without extrapolating calibration (#57).

* `step_sec_intrinsic_visc()` now converts recorded or explicit concentration
  units before calculating dL/g or mL/g. Composition uses calibrated two-channel
  responses to recover component mass fractions; UV/RI ratios distinguish raw
  detector responses from normalized concentrations. Inferred unit tags must
  persist when baking new data, and ratio smoothing preserves invalid points
  (#60).

* `step_sec_concentration()` reports clipped detector values with a structured
  `measure_sec_negative_signal_clipped` warning. Its per-curve diagnostics
  include positive and negative sampled areas and the negative share of
  absolute sampled signal area, which can differ greatly from the percentage
  of negative points. Normalization rejects nonfinite or fully underflowed
  output and preserves ordinary numerical results (#60).
* `step_sec_lals()`, `step_sec_rals()`, and `step_sec_mals()` accept finite,
  nonzero signed dn/dc values because the scattering optical constant squares
  dn/dc. LALS/RALS also require explicit or recorded concentration units,
  convert mg/mL to g/mL before single-angle calculations, and reject
  misaligned detector curves (#60).
* `step_sec_ri()` now respects `output_col`. A distinct output name preserves
  the raw source curve and its role; a named output requires one selected RI
  measure and cannot overwrite another column visible at that recipe stage.
  `NULL` keeps existing in-place processing, including multiple measures (#75).
* `step_sec_ri()` accepts finite, nonzero signed scalar and sample-specific
  refractive-index increments. Its documentation now distinguishes the
  positive raw-to-RI conversion gain from the signed increment, matching the
  existing calibration arithmetic (#60).
* `step_sec_ri()` and `step_sec_uv()` record the concentration units produced
  by detector normalization and allow mg/mL or g/mL output (#60).
