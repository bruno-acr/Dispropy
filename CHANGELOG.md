# Changelog

## 0.3.0

- Added `GPSFitWarning` for GPS fits based on fewer than 50 valid pairs.
- Added warnings when fitted GPS parameters are within 1% of optimization
  bounds, indicating possible weak identification.
- Added GPS reliability diagnostics to `DataFrame.attrs["gps_model"]`:
  `parameters_near_bounds`, `near_bound_parameters`, `valid_pair_count`, and
  `recommended_min_valid_pairs`.
- Added reliability guidance and recommended responses to GPS warnings.
- Added tests for small samples, boundary estimates, and healthy fits.
- Added an executable notebook demonstrating the complete public API and GPS
  reliability diagnostics.

## 0.2.0

- Added the Empirical Bayes Geometric Mean (`EBGM`) using the GPS model.
- Added posterior outputs `qn`, `ebgm`, `eb05`, and `eb95`.
- Added the `signal_ebgm` screening flag based on `EB05 > 2` and `A >= 3`.
- Added `scipy` as a runtime dependency.
- Added GPS model parameters and diagnostics to `DataFrame.attrs["gps_model"]`.
- Retained vectorized ROR, PRR, IC, IC025, and IC975 calculations.
- Expanded validation and automated tests.
