# measure.sec 0.0.0.9000

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

* `step_sec_concentration()` reports the count and percentage of negative
  detector values excluded during mass normalization, and sets those values to
  zero in the normalized curve (#60).
* `step_sec_lals()`, `step_sec_rals()`, and `step_sec_mals()` accept finite,
  nonzero signed dn/dc values because the scattering optical constant squares
  dn/dc. LALS/RALS also require explicit or recorded concentration units,
  convert mg/mL to g/mL before single-angle calculations, and reject
  misaligned detector curves (#60).
* `step_sec_ri()` and `step_sec_uv()` record the concentration units produced
  by detector normalization and allow mg/mL or g/mL output (#60).
