# Changelog

## 0.5.0

- **Breaking:** `calculate_ror`, `calculate_prr` and
  `calculate_disproportionality` now default to `correction=0.0` (no
  continuity correction) instead of `correction=0.5`. Code that relied on
  the implicit `0.5` correction must now pass `correction=0.5` explicitly
  to get the same numbers as before. With the new default, ROR/PRR (and
  their confidence intervals) will be `NaN`/infinite for any row where a
  contingency cell is 0, unless `correction` is set to a positive value.
  `calculate_disproportionality` forwards `correction` to both metrics
  unchanged, so this applies there too.
- `correction` now consistently rejects negative, non-numeric, `NaN` and
  infinite values with a clear `ValueError`, for all three functions. The
  original A, B, C and D columns are never mutated by the correction.
- Expanded the `calculate_ror`/`calculate_prr`/`calculate_disproportionality`
  docstrings and the README to explain the continuity correction, its
  default, that the same corrected counts drive the point estimate,
  standard error and confidence interval, and the implications of leaving
  it disabled on sparse tables.
- Clarified in the README that the ROR/PRR validation against
  `PhViD::ROR()`/`PhViD::PRR()` was performed with `correction=0.5`
  explicitly, since that is no longer the library's default.

## 0.4.0

- Fixed a bias in `calculate_ebgm`'s GPS hyperparameter fit: the likelihood
  now truncates at `n_star`, the smallest observed count in the fitting
  sample, following DuMouchel and Pregibon (2001). Real disproportionality
  tables normally list only pairs that were actually reported, so counts
  below that minimum (typically zero) are already excluded before reaching
  `calculate_ebgm`; fitting the untruncated mixture to such a sample biased
  all five hyperparameters and, through them, every EBGM/EB05/EB95/Qn value.
  Tables that already include zero-count rows are unaffected (`n_star` is 0
  and the correction is a no-op).
- Added `n_star` to `result.attrs["gps_model"]` for transparency.
- Cross-validated ROR, PRR, IC and EBGM against the R reference
  implementations `PhViD` and `openEBGM` on both a large simulated dataset
  and the real FDA CAERS dataset; see the README's Validation section.
- Clarified in the README that IC is the shrinkage approximation of Norén
  et al. (2013), not the full Dirichlet-based BCPNN posterior of Bate et
  al. (1998); the two correlate strongly but are not numerically identical,
  particularly for pairs with very few reports.
- Translated the README to English and added a References section citing
  the source of every implemented formula.

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
