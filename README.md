# dispropy

`dispropy` calculates disproportionality metrics for pharmacovigilance from
2x2 contingency tables stored in a `pandas.DataFrame`.

|                     | Event of interest | Other events |
|---------------------|--------------------|---------------|
| Drug of interest     | A                  | B             |
| Other drugs          | C                  | D             |

## Metrics

### ROR (Reporting Odds Ratio)

```
ROR = (A * D) / (B * C)
```

The log-scale 95% confidence interval uses the standard delta-method
variance `1/A + 1/B + 1/C + 1/D` [[1]](#references).

### PRR (Proportional Reporting Ratio)

```
PRR = [A / (A + B)] / [C / (C + D)]
```

The log-scale 95% confidence interval uses the delta-method variance
`1/A - 1/(A+B) + 1/C - 1/(C+D)` [[2]](#references).

### Continuity correction (ROR and PRR)

`calculate_ror` and `calculate_prr` both take a `correction` argument. A
continuity correction is a small constant added to sparse contingency
tables to avoid division by zero and undefined logarithms, so that ROR,
PRR and their confidence intervals can still be computed when a cell is 0.

- The library's default behavior is `correction=0`, i.e. no correction: the
  raw A, B, C and D counts are used as-is.
- If you omit `correction`, `0` is used automatically — calling
  `calculate_ror(data)` is equivalent to `calculate_ror(data,
  correction=0)`, and likewise for `calculate_prr`.
- You can apply the conventional continuity correction by passing
  `correction=0.5` explicitly.
- Whatever value you choose is added to all four contingency cells (A, B,
  C and D) before ROR/PRR and their respective confidence intervals are
  computed, so both the point estimate and the interval use the same
  corrected counts.
- `correction` only accepts nonnegative numbers; a negative value raises a
  `ValueError` with a clear message.

**Warning:** with the default `correction=0`, a table that has any cell
equal to 0 will produce an undefined ROR, PRR, their logarithms, or their
confidence intervals (infinite or `NaN`). Pass `correction=0.5` (or another
positive value) if your data may contain empty cells and you want finite
results.

#### Examples

```python
# Default behavior (no correction applied)
calculate_ror(data, "A", "B", "C", "D")
calculate_prr(data, "A", "B", "C", "D")

# Equivalent, explicit form
calculate_ror(data, "A", "B", "C", "D", correction=0)
calculate_prr(data, "A", "B", "C", "D", correction=0)

# With the conventional continuity correction
calculate_ror(data, "A", "B", "C", "D", correction=0.5)
calculate_prr(data, "A", "B", "C", "D", correction=0.5)
```

### IC (Information Component)

```
IC = log2((Obs + shrinkage) / (Exp + shrinkage))
Obs = A
Exp = ((A + B) * (A + C)) / (A + B + C + D)
```

with `shrinkage` defaulting to `0.5`. `IC025` and `IC975` (the 95%
credibility interval) use the simplified normal approximation of Norén,
Hopstadius and Bate (2013) [[5]](#references):

```
IC025 = IC - 3.3 * (Obs + 0.5)^(-0.5) - 2.0 * (Obs + 0.5)^(-1.5)
IC975 = IC + 2.4 * (Obs + 0.5)^(-0.5) - 0.5 * (Obs + 0.5)^(-1.5)
```

The IC concept itself was introduced as part of the Bayesian Confidence
Propagation Neural Network (BCPNN) by Bate et al. (1998) [[3]](#references),
which places a Dirichlet(1,1,1,1) prior on the four joint cell probabilities
and derives IC as a posterior expectation over the *full* 2x2 table. What
`dispropy` implements is the later, computationally simpler shrinkage
observed-to-expected approximation of [[5]](#references), which only needs
`Obs` and `Exp` (not the full Dirichlet posterior). The two correlate
strongly but are **not numerically identical** — see
[Validation](#validation-against-reference-implementations).

### EBGM (Empirical Bayes Geometric Mean)

DuMouchel's Gamma-Poisson Shrinker (GPS) method [[4]](#references). A
mixture of two Gamma distributions is fitted jointly to all rows by maximum
likelihood, giving five hyperparameters (`alpha1`, `beta1`, `alpha2`,
`beta2`, `weight`). For each pair, the output includes the posterior
mixture weight `qn`, `ebgm` (the posterior geometric mean of the relative
reporting rate), and `eb05`/`eb95` (the 5th and 95th percentiles of its
posterior distribution).

The hyperparameter likelihood is truncated at `n_star`, the smallest
observed count in the fitting sample, following DuMouchel and Pregibon
(2001) [[8]](#references). This matters in practice: a real
disproportionality table almost never lists the pairs with zero reports
(enumerating the full drug x event grid is normally infeasible), so the
fitting sample implicitly excludes counts below `n_star` before it ever
reaches `dispropy`. Fitting the untruncated GPS likelihood to such a sample
biases all five hyperparameters — this was caught during validation (see
below) and fixed in 0.4.0. If your table does include zero-count rows,
`n_star` is 0 and the correction has no effect. `n_star` is recorded in
`result.attrs["gps_model"]` for transparency.

## Validation against reference implementations

Formulas were checked against independent R reference implementations, not
just against `dispropy`'s own test suite, to confirm the numbers a
researcher gets are not just internally consistent but actually correct.
Two datasets were used: `PhViD`'s built-in simulated dataset (102,483
drug-event pairs) and the real FDA CAERS dataset (17,189 dietary-supplement
product-event pairs) that `openEBGM` ships with.

- **ROR and PRR** were compared against `PhViD::ROR()` and `PhViD::PRR()` on
  both datasets, calling `calculate_ror`/`calculate_prr` with
  `correction=0.5` to match `PhViD`'s continuity correction. Results matched
  to floating-point precision on the simulated data (correlation 1.000000,
  maximum relative difference ~1e-14) and on the real CAERS data
  (correlation 1.000000, maximum relative difference ~5e-15). Note that
  `calculate_ror`/`calculate_prr` default to `correction=0` (no correction)
  — this validation applies to the `correction=0.5` case specifically, not
  to the library's current default; see
  [Continuity correction](#continuity-correction-ror-and-prr).
- **IC** was compared against a reimplementation of the full Dirichlet-based
  BCPNN posterior of Bate et al. (1998) (as used by `PhViD::BCPNN`) on both
  datasets. The two correlate strongly but are not numerically identical,
  confirming the distinction described above: r ≈ 0.98 on the simulated
  data (99.3% agreement on the `ic025 > 0` signal call) and r ≈ 0.92 on the
  real CAERS data (99.4% signal agreement). The divergence is larger on
  CAERS because its counts are mostly small (`A` between 1 and 54), and
  that is exactly where the two approximations diverge most: mean |ΔIC| ≈
  0.47 for `A` between 1 and 2 versus ≈ 0.01 for `A` above 100 in the
  simulated data.
- **EBGM** was compared against the `openEBGM` package on the real CAERS
  data. This comparison is what surfaced the `n_star` truncation bug
  described above: before the fix, EBGM correlated at only 0.20 with the
  reference, and 23.5% of pairs were flagged as `EB05 > 2` signals that
  should not have been. After adding the truncation, `qn`, `ebgm`, `eb05`
  and `eb95` all correlate above 0.9999 with `openEBGM`, and agreement on
  the `EBGM > 2` signal call is 99.9%. Remaining differences are consistent
  with ordinary optimizer precision, not a systematic bias.

## Installation

```bash
pip install -e .
```

For development:

```bash
pip install -e ".[dev]"
pytest
```

## Usage

A runnable tutorial covering every metric, validation rule, signal flag and
GPS diagnostic is available at
[`examples/dispropy_tutorial.ipynb`](examples/dispropy_tutorial.ipynb).

```python
import pandas as pd
from dispropy import calculate_disproportionality

df = pd.DataFrame({
    "drug": ["Drug A", "Drug A", "Drug A", "Drug A"],
    "event": ["Event W", "Event X", "Event Y", "Event Z"],
    "A": [2, 5, 10, 20],
    "B": [98, 95, 90, 80],
    "C": [20, 25, 20, 30],
    "D": [880, 875, 880, 870],
})

result = calculate_disproportionality(
    df,
    a_col="A",
    b_col="B",
    c_col="C",
    d_col="D",
    metrics=("ror", "prr", "ic", "ebgm"),
    correction=0.5,
    shrinkage=0.5,
    add_signal_flags=True,
)

print(result)
print(result.attrs["gps_model"])
```

ROR, PRR and IC are calculated by default. EBGM must be requested explicitly
because it involves a global numerical fit and requires at least two pairs
with a positive expected count. The fitted GPS parameters are stored in
`result.attrs["gps_model"]`.

`calculate_ror`, `calculate_prr`, `calculate_ic` and `calculate_ebgm` can
also be called individually. Column names for A, B, C and D are free-form.

## Signal flags and interpretation

With `add_signal_flags=True`, the screening criteria are:

- `signal_ror`: `ror_lower_95 > 1`;
- `signal_prr`: `prr_lower_95 > 1`;
- `signal_ic`: `ic025 > 0`;
- `signal_ebgm`: `eb05 > 2` and at least three observed cases (`A >= 3`).

All four use a strict `>` rather than `>=`, because the flag marks whether
the lower bound of the interval excludes the null value (ROR/PRR = 1, IC =
0, EB05 = 2): a bound exactly at the null value does not reject it. This
mirrors the "point estimate minus 1.96 SE" screening approach compared
against the IC in [[2]](#references), rather than the original combined
criterion from [[6]](#references) (`PRR >= 2`, `chi-squared >= 4`, `N >= 3`),
which `dispropy` does not implement. If you need the latter, compute it from
the columns and `prr` output yourself.

ROR or PRR above 1 and IC above 0 suggest reporting higher than the
comparator or than expected. These metrics and flags indicate only
statistical disproportionality. They do not demonstrate causality, and
clinical and pharmacological assessment is still required.

## Limitations

- Does not correct for biases inherent to spontaneous reporting.
- Does not estimate incidence or absolute risk.
- The GPS fit depends on the number and composition of the analyzed pairs.
- Does not replace clinical assessment.
- Does not implement stratification by age, sex, country or period.

## When results may not be reliable

The library's validation confirms that the columns exist, are numeric,
contain no missing or negative values, and that no row is entirely empty.
This does not prove that `A`, `B`, `C` and `D` were constructed correctly.
Users must check the counting unit, deduplication, drug and event
definitions, comparator population and consistency of totals.

ROR, PRR and IC can be statistically unstable when counts are small. IC's
shrinkage is applied by default, and ROR/PRR's continuity correction can be
enabled with `correction=0.5`; both avoid undefined operations and reduce
part of that instability, but they do not create information or replace an
assessment of precision and clinical relevance. IC additionally trades some
accuracy for speed relative to the full BCPNN posterior (see
[Validation](#validation-against-reference-implementations)); for pairs
with very few reports, treat IC as a fast screening approximation rather
than a substitute for the full Bayesian posterior.

EBGM fits five hyperparameters using all valid pairs. The library emits
`GPSFitWarning` when there are fewer than 50 valid pairs or when a parameter
ends up within 1% of an optimizer bound. The 50-pair threshold is a
conservative operational diagnostic, equivalent to 10 pairs per
hyperparameter, and not a formal cutoff established in the literature. GPS
was developed for large frequency tables, and reference implementations also
check convergence within the parameter space and stability across solutions.
See [[4]](#references) and [[7]](#references).

When a warning appears, inspect `result.attrs["gps_model"]`. The fields
`parameters_near_bounds` and `near_bound_parameters` indicate proximity to
the bounds; `valid_pair_count` reports how many pairs supported the fit;
`n_star` reports the zero-truncation floor actually used (see
[EBGM](#ebgm-empirical-bayes-geometric-mean) above). Do not interpret EBGM
in isolation: review how the table was built, increase and diversify the
set of pairs when possible, and run a sensitivity analysis or independent
statistical validation before using the result for a decision.

## References

1. Rothman KJ, Lanes S, Sacks ST. The reporting odds ratio and its
   advantages over the proportional reporting ratio. *Pharmacoepidemiol
   Drug Saf.* 2004;13(8):519-523.
   [doi:10.1002/pds.1001](https://doi.org/10.1002/pds.1001)
2. van Puijenbroek EP, Bate A, Leufkens HGM, Lindquist M, Orre R, Egberts
   ACG. A comparison of measures of disproportionality for signal detection
   in spontaneous reporting systems for adverse drug reactions.
   *Pharmacoepidemiol Drug Saf.* 2002;11(1):3-10.
   [doi:10.1002/pds.668](https://doi.org/10.1002/pds.668)
3. Bate A, Lindquist M, Edwards IR, Olsson S, Orre R, Lansner A, De Freitas
   RM. A Bayesian neural network method for adverse drug reaction signal
   generation. *Eur J Clin Pharmacol.* 1998;54(4):315-321.
   [doi:10.1007/s002280050466](https://doi.org/10.1007/s002280050466)
4. DuMouchel W. Bayesian data mining in large frequency tables, with an
   application to the FDA spontaneous reporting system. *Am Stat.*
   1999;53(3):177-190.
   [doi:10.1080/00031305.1999.10474456](https://doi.org/10.1080/00031305.1999.10474456)
5. Norén GN, Hopstadius J, Bate A. Shrinkage observed-to-expected ratios for
   robust and transparent large-scale pattern discovery. *Stat Methods Med
   Res.* 2013;22(1):57-69.
   [doi:10.1177/0962280211403604](https://doi.org/10.1177/0962280211403604)
6. Evans SJ, Waller PC, Davis S. Use of proportional reporting ratios (PRRs)
   for signal generation from spontaneous adverse drug reaction reports.
   *Pharmacoepidemiol Drug Saf.* 2001;10(6):483-486.
   [doi:10.1002/pds.677](https://doi.org/10.1002/pds.677)
7. Canida T, Ihrie J. openEBGM: An R implementation of the Gamma-Poisson
   shrinker data mining model. *R J.* 2017;9(2):84-97.
   [journal.r-project.org/articles/RJ-2017-063](https://journal.r-project.org/articles/RJ-2017-063/)
8. DuMouchel W, Pregibon D. Empirical Bayes screening for multi-item
   associations. In *Proceedings of the Seventh ACM SIGKDD International
   Conference on Knowledge Discovery and Data Mining* (KDD '01), 2001,
   pp. 67-76.
   [doi:10.1145/502512.502526](https://doi.org/10.1145/502512.502526)
