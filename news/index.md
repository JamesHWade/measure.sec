# Changelog

## measure.sec 0.0.0.9000

- [`step_sec_concentration()`](https://jameshwade.github.io/measure.sec/reference/step_sec_concentration.md)
  reports the count and percentage of negative detector values excluded
  during mass normalization, and sets those values to zero in the
  normalized curve
  ([\#60](https://github.com/JamesHWade/measure.sec/issues/60)).
- [`step_sec_lals()`](https://jameshwade.github.io/measure.sec/reference/step_sec_lals.md),
  [`step_sec_rals()`](https://jameshwade.github.io/measure.sec/reference/step_sec_rals.md),
  and
  [`step_sec_mals()`](https://jameshwade.github.io/measure.sec/reference/step_sec_mals.md)
  accept finite, nonzero signed dn/dc values because the scattering
  optical constant squares dn/dc. LALS/RALS also require explicit or
  recorded concentration units, convert mg/mL to g/mL before
  single-angle calculations, and reject misaligned detector curves
  ([\#60](https://github.com/JamesHWade/measure.sec/issues/60)).
- [`step_sec_ri()`](https://jameshwade.github.io/measure.sec/reference/step_sec_ri.md)
  and
  [`step_sec_uv()`](https://jameshwade.github.io/measure.sec/reference/step_sec_uv.md)
  record the concentration units produced by detector normalization and
  allow mg/mL or g/mL output
  ([\#60](https://github.com/JamesHWade/measure.sec/issues/60)).
