# measure.sec 0.0.0.9000

* `step_sec_concentration()` reports the count and percentage of negative
  detector values excluded during mass normalization, and sets those values to
  zero in the normalized curve (#60).
* `step_sec_lals()` and `step_sec_rals()` require explicit or recorded
  concentration units, convert mg/mL to g/mL before single-angle scattering
  calculations, and reject misaligned detector curves (#60).
* `step_sec_ri()` and `step_sec_uv()` record the concentration units produced
  by detector normalization and allow mg/mL or g/mL output (#60).
