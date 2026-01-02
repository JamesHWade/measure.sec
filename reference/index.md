# Package index

## Visualization

Publication-ready plots for SEC/GPC results

- [`sec_results()`](https://jameshwade.github.io/measure-sec/reference/sec_results.md)
  : Create SEC Results Object
- [`autoplot(`*`<sec_results>`*`)`](https://jameshwade.github.io/measure-sec/reference/autoplot.sec_results.md)
  : Automatic Plot for SEC Results
- [`plot_sec()`](https://jameshwade.github.io/measure-sec/reference/plot_sec.md)
  : Quick SEC Plot
- [`plot_sec_chromatogram()`](https://jameshwade.github.io/measure-sec/reference/plot_sec_chromatogram.md)
  : Plot SEC Chromatogram
- [`plot_sec_mwd()`](https://jameshwade.github.io/measure-sec/reference/plot_sec_mwd.md)
  : Plot Molecular Weight Distribution
- [`plot_sec_multidetector()`](https://jameshwade.github.io/measure-sec/reference/plot_sec_multidetector.md)
  : Plot Multi-Detector SEC Overlay
- [`plot_sec_conformation()`](https://jameshwade.github.io/measure-sec/reference/plot_sec_conformation.md)
  : Plot SEC Conformation Data
- [`plot_sec_calibration()`](https://jameshwade.github.io/measure-sec/reference/plot_sec_calibration.md)
  : Plot SEC Calibration Curve
- [`plot_sec_composition()`](https://jameshwade.github.io/measure-sec/reference/plot_sec_composition.md)
  : Plot SEC Composition Distribution

## Detector Processing

Recipe steps for SEC detector signal processing

- [`step_sec_ri()`](https://jameshwade.github.io/measure-sec/reference/step_sec_ri.md)
  : RI Detector Processing for SEC
- [`step_sec_uv()`](https://jameshwade.github.io/measure-sec/reference/step_sec_uv.md)
  : UV Detector Processing for SEC
- [`step_sec_dad()`](https://jameshwade.github.io/measure-sec/reference/step_sec_dad.md)
  : Diode Array Detector Processing for SEC
- [`step_sec_mals()`](https://jameshwade.github.io/measure-sec/reference/step_sec_mals.md)
  : Multi-Angle Light Scattering Processing for SEC
- [`step_sec_lals()`](https://jameshwade.github.io/measure-sec/reference/step_sec_lals.md)
  : Low-Angle Light Scattering Processing for SEC
- [`step_sec_rals()`](https://jameshwade.github.io/measure-sec/reference/step_sec_rals.md)
  : Right-Angle Light Scattering Processing for SEC
- [`step_sec_dls()`](https://jameshwade.github.io/measure-sec/reference/step_sec_dls.md)
  : Dynamic Light Scattering Processing for SEC
- [`step_sec_viscometer()`](https://jameshwade.github.io/measure-sec/reference/step_sec_viscometer.md)
  : Differential Viscometer Processing for SEC
- [`step_sec_detector_delay()`](https://jameshwade.github.io/measure-sec/reference/step_sec_detector_delay.md)
  : Correct Inter-Detector Volume Delays

## Signal Processing

Baseline correction and signal enhancement

- [`step_sec_baseline()`](https://jameshwade.github.io/measure-sec/reference/step_sec_baseline.md)
  : SEC/GPC Baseline Correction
- [`step_sec_band_broadening()`](https://jameshwade.github.io/measure-sec/reference/step_sec_band_broadening.md)
  : Band Broadening Correction for SEC
- [`step_sec_concentration()`](https://jameshwade.github.io/measure-sec/reference/step_sec_concentration.md)
  : Convert Detector Signal to Concentration
- [`estimate_sigma()`](https://jameshwade.github.io/measure-sec/reference/estimate_sigma.md)
  : Estimate Spreading Parameter from Narrow Standard

## Calibration

Molecular weight calibration methods

- [`step_sec_conventional_cal()`](https://jameshwade.github.io/measure-sec/reference/step_sec_conventional_cal.md)
  : Conventional Calibration for SEC Using Narrow Standards
- [`step_sec_broad_standard()`](https://jameshwade.github.io/measure-sec/reference/step_sec_broad_standard.md)
  : Broad Standard Calibration for SEC/GPC
- [`step_sec_universal_cal()`](https://jameshwade.github.io/measure-sec/reference/step_sec_universal_cal.md)
  : Universal Calibration for SEC

## Molecular Weight Calculations

MW averages, distributions, and fractions

- [`step_sec_mw_averages()`](https://jameshwade.github.io/measure-sec/reference/step_sec_mw_averages.md)
  : Calculate Molecular Weight Averages for SEC/GPC
- [`step_sec_mw_distribution()`](https://jameshwade.github.io/measure-sec/reference/step_sec_mw_distribution.md)
  : Generate Molecular Weight Distribution Curve
- [`step_sec_mw_fractions()`](https://jameshwade.github.io/measure-sec/reference/step_sec_mw_fractions.md)
  : Calculate Molecular Weight Fractions for SEC/GPC

## Composition Analysis

Copolymer and composition analysis

- [`step_sec_composition()`](https://jameshwade.github.io/measure-sec/reference/step_sec_composition.md)
  : Calculate Copolymer Composition from Detector Signals
- [`step_sec_uv_ri_ratio()`](https://jameshwade.github.io/measure-sec/reference/step_sec_uv_ri_ratio.md)
  : Calculate UV/RI Ratio for Composition Analysis
- [`step_sec_intrinsic_visc()`](https://jameshwade.github.io/measure-sec/reference/step_sec_intrinsic_visc.md)
  : Calculate Intrinsic Viscosity for SEC

## Protein SEC

Biopharmaceutical SEC analysis

- [`step_sec_protein()`](https://jameshwade.github.io/measure-sec/reference/step_sec_protein.md)
  : Protein SEC Analysis Workflow
- [`step_sec_aggregates()`](https://jameshwade.github.io/measure-sec/reference/step_sec_aggregates.md)
  : Quantify Protein Aggregates and Fragments in SEC
- [`step_sec_oligomer()`](https://jameshwade.github.io/measure-sec/reference/step_sec_oligomer.md)
  : Oligomeric Species Analysis for Protein SEC

## Quality Control Functions

Functions for system suitability testing and QC

- [`measure_sec_resolution()`](https://jameshwade.github.io/measure-sec/reference/measure_sec_resolution.md)
  : Calculate Peak Resolution
- [`measure_sec_plate_count()`](https://jameshwade.github.io/measure-sec/reference/measure_sec_plate_count.md)
  : Calculate Theoretical Plate Count
- [`measure_sec_asymmetry()`](https://jameshwade.github.io/measure-sec/reference/measure_sec_asymmetry.md)
  : Calculate Peak Asymmetry Factor
- [`measure_sec_recovery()`](https://jameshwade.github.io/measure-sec/reference/measure_sec_recovery.md)
  : Calculate Mass Recovery
- [`measure_sec_suitability()`](https://jameshwade.github.io/measure-sec/reference/measure_sec_suitability.md)
  : System Suitability Test for SEC
- [`measure_sec_column_performance()`](https://jameshwade.github.io/measure-sec/reference/measure_sec_column_performance.md)
  : Column Performance Metrics for SEC

## Polymer Analysis Functions

Functions for polymer characterization

- [`measure_mh_parameters()`](https://jameshwade.github.io/measure-sec/reference/measure_mh_parameters.md)
  : Estimate Mark-Houwink Parameters
- [`measure_branching_index()`](https://jameshwade.github.io/measure-sec/reference/measure_branching_index.md)
  : Calculate Branching Index
- [`measure_branching_frequency()`](https://jameshwade.github.io/measure-sec/reference/measure_branching_frequency.md)
  : Calculate Branching Frequency from Zimm-Stockmayer Theory
- [`measure_branching_model_comparison()`](https://jameshwade.github.io/measure-sec/reference/measure_branching_model_comparison.md)
  : Compare Branching Models
- [`measure_rg_mw_scaling()`](https://jameshwade.github.io/measure-sec/reference/measure_rg_mw_scaling.md)
  : Fit Rg-MW Scaling Relationship
- [`measure_conformation_data()`](https://jameshwade.github.io/measure-sec/reference/measure_conformation_data.md)
  : Create Conformation Plot Data

## Export Functions

Data extraction and export utilities

- [`measure_sec_slice_table()`](https://jameshwade.github.io/measure-sec/reference/measure_sec_slice_table.md)
  : Extract Slice-by-Slice SEC Data
- [`measure_sec_summary_table()`](https://jameshwade.github.io/measure-sec/reference/measure_sec_summary_table.md)
  : Generate SEC Summary Table
- [`measure_sec_compare()`](https://jameshwade.github.io/measure-sec/reference/measure_sec_compare.md)
  : Compare Multiple SEC Samples

## Datasets

Example datasets for SEC analysis

- [`sec_triple_detect`](https://jameshwade.github.io/measure-sec/reference/sec_triple_detect.md)
  : Multi-Detector SEC Data
- [`sec_calibration_standards`](https://jameshwade.github.io/measure-sec/reference/sec_calibration_standards.md)
  : SEC Calibration Standards - Complete Dataset
- [`sec_ps_standards`](https://jameshwade.github.io/measure-sec/reference/sec_ps_standards.md)
  : SEC Polystyrene Calibration Standards
- [`sec_pmma_standards`](https://jameshwade.github.io/measure-sec/reference/sec_pmma_standards.md)
  : SEC PMMA Calibration Standards
- [`sec_copolymer`](https://jameshwade.github.io/measure-sec/reference/sec_copolymer.md)
  : Copolymer SEC Data for Composition Analysis
- [`sec_protein`](https://jameshwade.github.io/measure-sec/reference/sec_protein.md)
  : Protein SEC Data with Aggregates
- [`sec_branched`](https://jameshwade.github.io/measure-sec/reference/sec_branched.md)
  : Branched Polymer SEC Data
- [`sec_system_suitability`](https://jameshwade.github.io/measure-sec/reference/sec_system_suitability.md)
  : System Suitability Test Data
