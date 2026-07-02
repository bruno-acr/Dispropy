"""Vectorized disproportionality metrics for pharmacovigilance."""

from collections.abc import Iterable, Sequence
from typing import Hashable
import warnings

import numpy as np
import pandas as pd
from scipy import optimize, special, stats

from .diagnostics import GPSFitWarning
from .validation import get_contingency_arrays

_VALID_METRICS = frozenset({"ror", "prr", "ic", "ebgm"})
_GPS_INITIAL = (0.2, 0.1, 2.0, 4.0, 1.0 / 3.0)
_GPS_BOUNDS = (
    (1e-4, 100.0),
    (1e-4, 100.0),
    (1e-4, 100.0),
    (1e-4, 100.0),
    (1e-4, 1.0 - 1e-4),
)
_GPS_PARAMETER_NAMES = ("alpha1", "beta1", "alpha2", "beta2", "weight")
_GPS_MIN_VALID_PAIRS = 50
_GPS_BOUND_PROXIMITY = 0.01


def _output_frame(df: pd.DataFrame, inplace: bool) -> pd.DataFrame:
    return df if inplace else df.copy()


def _validate_positive_parameter(value: float, name: str, *, allow_zero: bool) -> float:
    if not isinstance(value, (int, float, np.number)) or not np.isfinite(value):
        raise ValueError(f"{name} must be a finite number.")
    if value < 0 or (not allow_zero and value == 0):
        qualifier = "nonnegative" if allow_zero else "greater than zero"
        raise ValueError(f"{name} must be {qualifier}.")
    return float(value)


def _expected_counts(
    a: pd.Series, b: pd.Series, c: pd.Series, d: pd.Series
) -> pd.Series:
    total = a + b + c + d
    if (total == 0).any():
        raise ValueError("A+B+C+D must be greater than zero for every row.")
    return ((a + b) * (a + c)) / total


def calculate_ror(
    df: pd.DataFrame,
    a_col: Hashable,
    b_col: Hashable,
    c_col: Hashable,
    d_col: Hashable,
    correction: float = 0.5,
    inplace: bool = False,
) -> pd.DataFrame:
    """Calculate ROR and its log-scale 95% confidence interval."""
    correction = _validate_positive_parameter(correction, "correction", allow_zero=True)
    a, b, c, d = get_contingency_arrays(df, a_col, b_col, c_col, d_col)
    ac, bc, cc, dc = (value + correction for value in (a, b, c, d))
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        ror = (ac * dc) / (bc * cc)
        log_ror = np.log(ror)
        se = np.sqrt(1.0 / ac + 1.0 / bc + 1.0 / cc + 1.0 / dc)
    result = _output_frame(df, inplace)
    result["ror"] = ror
    result["log_ror"] = log_ror
    result["se_log_ror"] = se
    result["ror_lower_95"] = np.exp(log_ror - 1.96 * se)
    result["ror_upper_95"] = np.exp(log_ror + 1.96 * se)
    return result


def calculate_prr(
    df: pd.DataFrame,
    a_col: Hashable,
    b_col: Hashable,
    c_col: Hashable,
    d_col: Hashable,
    correction: float = 0.5,
    inplace: bool = False,
) -> pd.DataFrame:
    """Calculate PRR and its approximate log-scale 95% confidence interval."""
    correction = _validate_positive_parameter(correction, "correction", allow_zero=True)
    a, b, c, d = get_contingency_arrays(df, a_col, b_col, c_col, d_col)
    ac, bc, cc, dc = (value + correction for value in (a, b, c, d))
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        prr = (ac / (ac + bc)) / (cc / (cc + dc))
        log_prr = np.log(prr)
        se = np.sqrt(1.0 / ac - 1.0 / (ac + bc) + 1.0 / cc - 1.0 / (cc + dc))
    result = _output_frame(df, inplace)
    result["prr"] = prr
    result["log_prr"] = log_prr
    result["se_log_prr"] = se
    result["prr_lower_95"] = np.exp(log_prr - 1.96 * se)
    result["prr_upper_95"] = np.exp(log_prr + 1.96 * se)
    return result


def calculate_ic(
    df: pd.DataFrame,
    a_col: Hashable,
    b_col: Hashable,
    c_col: Hashable,
    d_col: Hashable,
    shrinkage: float = 0.5,
    inplace: bool = False,
) -> pd.DataFrame:
    """Calculate IC, IC025 and IC975 from observed and expected counts."""
    shrinkage = _validate_positive_parameter(shrinkage, "shrinkage", allow_zero=False)
    a, b, c, d = get_contingency_arrays(df, a_col, b_col, c_col, d_col)
    expected = _expected_counts(a, b, c, d)
    adjusted = a + shrinkage
    ic = np.log2(adjusted / (expected + shrinkage))
    result = _output_frame(df, inplace)
    result["observed_count"] = a
    result["expected_count"] = expected
    result["ic"] = ic
    result["ic025"] = ic - 3.3 * np.power(adjusted, -0.5) - 2.0 * np.power(adjusted, -1.5)
    result["ic975"] = ic + 2.4 * np.power(adjusted, -0.5) - 0.5 * np.power(adjusted, -1.5)
    return result


def _negative_log_likelihood(
    parameters: np.ndarray, observed: np.ndarray, expected: np.ndarray
) -> float:
    alpha1, beta1, alpha2, beta2, weight = parameters
    p1 = beta1 / (beta1 + expected)
    p2 = beta2 / (beta2 + expected)
    components = np.vstack(
        (
            np.log(weight) + stats.nbinom.logpmf(observed, alpha1, p1),
            np.log1p(-weight) + stats.nbinom.logpmf(observed, alpha2, p2),
        )
    )
    value = -float(np.sum(special.logsumexp(components, axis=0)))
    return value if np.isfinite(value) else 1e300


def _parameters_near_bounds(parameters: Sequence[float]) -> list[str]:
    """Return parameter names within 1% of either optimization bound."""
    affected = []
    for name, value, (lower, upper) in zip(
        _GPS_PARAMETER_NAMES, parameters, _GPS_BOUNDS
    ):
        tolerance = _GPS_BOUND_PROXIMITY * (upper - lower)
        if value - lower <= tolerance or upper - value <= tolerance:
            affected.append(name)
    return affected


def _posterior_quantile(
    quantile: float,
    observed: float,
    expected: float,
    parameters: Sequence[float],
    posterior_weight: float,
) -> float:
    alpha1, beta1, alpha2, beta2, _ = parameters

    def objective(value: float) -> float:
        first = stats.gamma.cdf(value, a=alpha1 + observed, scale=1.0 / (beta1 + expected))
        second = stats.gamma.cdf(value, a=alpha2 + observed, scale=1.0 / (beta2 + expected))
        return posterior_weight * first + (1.0 - posterior_weight) * second - quantile

    upper = 1.0
    while objective(upper) < 0 and upper < 1e12:
        upper *= 10.0
    return float(optimize.brentq(objective, 0.0, upper, xtol=1e-10, rtol=1e-10))


def calculate_ebgm(
    df: pd.DataFrame,
    a_col: Hashable,
    b_col: Hashable,
    c_col: Hashable,
    d_col: Hashable,
    inplace: bool = False,
    initial_parameters: Sequence[float] = _GPS_INITIAL,
) -> pd.DataFrame:
    """Fit the two-gamma GPS model and calculate Qn, EBGM, EB05 and EB95.

    The prior is estimated jointly from all rows. At least two rows with a
    positive expected count are required. Fitted parameters and optimizer
    diagnostics are stored in ``result.attrs["gps_model"]``.
    """
    a, b, c, d = get_contingency_arrays(df, a_col, b_col, c_col, d_col)
    expected = _expected_counts(a, b, c, d)
    valid = np.isfinite(a.to_numpy()) & np.isfinite(expected.to_numpy()) & (expected.to_numpy() > 0)
    if valid.sum() < 2:
        raise ValueError("EBGM requires at least two rows with expected_count greater than zero.")
    valid_pair_count = int(valid.sum())
    if valid_pair_count < _GPS_MIN_VALID_PAIRS:
        warnings.warn(
            "The GPS model is being fitted with "
            f"{valid_pair_count} valid pairs; fewer than the recommended "
            f"operational minimum of {_GPS_MIN_VALID_PAIRS}. The result is "
            "still calculated, but the five mixture parameters may be weakly "
            "identified and uncertainty may be greater.",
            GPSFitWarning,
            stacklevel=2,
        )
    initial = np.asarray(initial_parameters, dtype=float)
    if initial.shape != (5,):
        raise ValueError("initial_parameters must contain alpha1, beta1, alpha2, beta2, and weight.")

    observed_fit = a.to_numpy()[valid]
    expected_fit = expected.to_numpy()[valid]
    fit = optimize.minimize(
        _negative_log_likelihood,
        x0=initial,
        args=(observed_fit, expected_fit),
        method="L-BFGS-B",
        bounds=_GPS_BOUNDS,
    )
    if not fit.success or not np.all(np.isfinite(fit.x)):
        raise RuntimeError(f"GPS model optimization failed: {fit.message}")

    near_bound_parameters = _parameters_near_bounds(fit.x)
    if near_bound_parameters:
        warnings.warn(
            "GPS optimization converged, but parameter(s) are within 1% of "
            "an optimization bound: "
            f"{', '.join(near_bound_parameters)}. The fit may be poorly "
            "identified; inspect data heterogeneity and consider a larger, "
            "more representative set of drug-event pairs.",
            GPSFitWarning,
            stacklevel=2,
        )

    alpha1, beta1, alpha2, beta2, weight = fit.x
    observed = a.to_numpy()
    expected_array = expected.to_numpy()
    log_component1 = np.log(weight) + stats.nbinom.logpmf(
        observed, alpha1, beta1 / (beta1 + expected_array)
    )
    log_component2 = np.log1p(-weight) + stats.nbinom.logpmf(
        observed, alpha2, beta2 / (beta2 + expected_array)
    )
    qn = np.exp(log_component1 - np.logaddexp(log_component1, log_component2))
    log_lambda = qn * (special.digamma(alpha1 + observed) - np.log(beta1 + expected_array))
    log_lambda += (1.0 - qn) * (
        special.digamma(alpha2 + observed) - np.log(beta2 + expected_array)
    )
    ebgm = np.exp(log_lambda)
    eb05 = np.fromiter(
        (_posterior_quantile(0.05, n, e, fit.x, q) for n, e, q in zip(observed, expected_array, qn)),
        dtype=float,
        count=len(df),
    )
    eb95 = np.fromiter(
        (_posterior_quantile(0.95, n, e, fit.x, q) for n, e, q in zip(observed, expected_array, qn)),
        dtype=float,
        count=len(df),
    )

    result = _output_frame(df, inplace)
    result["observed_count"] = a
    result["expected_count"] = expected
    result["qn"] = qn
    result["ebgm"] = ebgm
    result["eb05"] = eb05
    result["eb95"] = eb95
    result.attrs["gps_model"] = {
        "alpha1": float(alpha1),
        "beta1": float(beta1),
        "alpha2": float(alpha2),
        "beta2": float(beta2),
        "weight": float(weight),
        "converged": bool(fit.success),
        "parameters_near_bounds": bool(near_bound_parameters),
        "near_bound_parameters": near_bound_parameters,
        "valid_pair_count": valid_pair_count,
        "recommended_min_valid_pairs": _GPS_MIN_VALID_PAIRS,
        "log_likelihood": float(-fit.fun),
    }
    return result


def calculate_disproportionality(
    df: pd.DataFrame,
    a_col: Hashable,
    b_col: Hashable,
    c_col: Hashable,
    d_col: Hashable,
    correction: float = 0.5,
    shrinkage: float = 0.5,
    metrics: Iterable[str] = ("ror", "prr", "ic"),
    add_signal_flags: bool = False,
    inplace: bool = False,
) -> pd.DataFrame:
    """Calculate selected disproportionality metrics for a DataFrame.

    Supported metrics are ``ror``, ``prr``, ``ic`` and ``ebgm``. EBGM is not
    enabled by default because it fits a model jointly and is more expensive.
    """
    get_contingency_arrays(df, a_col, b_col, c_col, d_col)
    if isinstance(metrics, str):
        raise TypeError("metrics must be an iterable of metric names, not a string.")
    selected = tuple(dict.fromkeys(str(metric).lower() for metric in metrics))
    invalid = sorted(set(selected) - _VALID_METRICS)
    if invalid:
        raise ValueError(f"Invalid metric(s): {invalid}. Valid metrics are: {sorted(_VALID_METRICS)}.")

    result = _output_frame(df, inplace)
    calculators = {
        "ror": lambda: calculate_ror(result, a_col, b_col, c_col, d_col, correction, True),
        "prr": lambda: calculate_prr(result, a_col, b_col, c_col, d_col, correction, True),
        "ic": lambda: calculate_ic(result, a_col, b_col, c_col, d_col, shrinkage, True),
        "ebgm": lambda: calculate_ebgm(result, a_col, b_col, c_col, d_col, True),
    }
    for metric in selected:
        calculators[metric]()

    if add_signal_flags:
        if "ror" in selected:
            result["signal_ror"] = result["ror_lower_95"] > 1.0
        if "prr" in selected:
            result["signal_prr"] = result["prr_lower_95"] > 1.0
        if "ic" in selected:
            result["signal_ic"] = result["ic025"] > 0.0
        if "ebgm" in selected:
            result["signal_ebgm"] = (result["eb05"] > 2.0) & (result[a_col] >= 3)
    return result
