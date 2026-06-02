"""Vectorized disproportionality metrics for pharmacovigilance."""

from __future__ import annotations

from collections.abc import Iterable
from typing import Hashable, Literal

import numpy as np
import pandas as pd

from .validation import get_contingency_arrays


ColumnName = Hashable
MetricName = Literal["ror", "prr", "ic"]
VALID_METRICS = {"ror", "prr", "ic"}


def _prepare_result_frame(df: pd.DataFrame, inplace: bool) -> pd.DataFrame:
    return df if inplace else df.copy()


def _validate_positive_parameter(value: float, name: str) -> None:
    if value <= 0:
        raise ValueError(f"{name} must be greater than zero.")


def _validate_non_negative_parameter(value: float, name: str) -> None:
    if value < 0:
        raise ValueError(f"{name} must be greater than or equal to zero.")


def _validate_total_count(
    a: pd.Series,
    b: pd.Series,
    c: pd.Series,
    d: pd.Series,
) -> None:
    if ((a + b + c + d) == 0).any():
        raise ValueError("The sum A+B+C+D cannot be zero for any row.")


def calculate_ror(
    df: pd.DataFrame,
    a_col: ColumnName,
    b_col: ColumnName,
    c_col: ColumnName,
    d_col: ColumnName,
    correction: float = 0.5,
    inplace: bool = False,
) -> pd.DataFrame:
    """Calculate Reporting Odds Ratio (ROR) and its 95% confidence interval.

    ROR is calculated as ``(A * D) / (B * C)`` after applying a continuity
    correction to A, B, C, and D.
    """
    _validate_non_negative_parameter(correction, "correction")
    result = _prepare_result_frame(df, inplace)
    a, b, c, d = get_contingency_arrays(result, a_col, b_col, c_col, d_col)

    a_calc = a + correction
    b_calc = b + correction
    c_calc = c + correction
    d_calc = d + correction

    ror = (a_calc * d_calc) / (b_calc * c_calc)
    log_ror = np.log(ror)
    se_log_ror = np.sqrt(
        1 / a_calc + 1 / b_calc + 1 / c_calc + 1 / d_calc
    )

    result["ror"] = ror
    result["log_ror"] = log_ror
    result["se_log_ror"] = se_log_ror
    result["ror_lower_95"] = np.exp(log_ror - 1.96 * se_log_ror)
    result["ror_upper_95"] = np.exp(log_ror + 1.96 * se_log_ror)
    return result


def calculate_prr(
    df: pd.DataFrame,
    a_col: ColumnName,
    b_col: ColumnName,
    c_col: ColumnName,
    d_col: ColumnName,
    correction: float = 0.5,
    inplace: bool = False,
) -> pd.DataFrame:
    """Calculate Proportional Reporting Ratio (PRR) and its 95% interval.

    PRR is calculated as ``[A / (A + B)] / [C / (C + D)]`` after applying a
    continuity correction to A, B, C, and D.
    """
    _validate_non_negative_parameter(correction, "correction")
    result = _prepare_result_frame(df, inplace)
    a, b, c, d = get_contingency_arrays(result, a_col, b_col, c_col, d_col)

    a_calc = a + correction
    b_calc = b + correction
    c_calc = c + correction
    d_calc = d + correction

    prr = (a_calc / (a_calc + b_calc)) / (
        c_calc / (c_calc + d_calc)
    )
    log_prr = np.log(prr)
    se_log_prr = np.sqrt(
        (1 / a_calc)
        - (1 / (a_calc + b_calc))
        + (1 / c_calc)
        - (1 / (c_calc + d_calc))
    )

    result["prr"] = prr
    result["log_prr"] = log_prr
    result["se_log_prr"] = se_log_prr
    result["prr_lower_95"] = np.exp(log_prr - 1.96 * se_log_prr)
    result["prr_upper_95"] = np.exp(log_prr + 1.96 * se_log_prr)
    return result


def calculate_ic(
    df: pd.DataFrame,
    a_col: ColumnName,
    b_col: ColumnName,
    c_col: ColumnName,
    d_col: ColumnName,
    shrinkage: float = 0.5,
    inplace: bool = False,
) -> pd.DataFrame:
    """Calculate Information Component (IC), IC025, and IC975.

    The expected count is calculated from the original 2x2 table values:
    ``((A + B) * (A + C)) / (A + B + C + D)``. Shrinkage is applied only to
    the observed-to-expected ratio and interval formulas.
    """
    _validate_positive_parameter(shrinkage, "shrinkage")
    result = _prepare_result_frame(df, inplace)
    a, b, c, d = get_contingency_arrays(result, a_col, b_col, c_col, d_col)
    _validate_total_count(a, b, c, d)

    observed = a
    total = a + b + c + d
    expected = ((a + b) * (a + c)) / total
    observed_shrunk = observed + shrinkage

    ic = np.log2(observed_shrunk / (expected + shrinkage))
    ic025 = (
        ic
        - 3.3 * np.power(observed_shrunk, -0.5)
        - 2.0 * np.power(observed_shrunk, -1.5)
    )
    ic975 = (
        ic
        + 2.4 * np.power(observed_shrunk, -0.5)
        - 0.5 * np.power(observed_shrunk, -1.5)
    )

    result["observed_count"] = observed
    result["expected_count"] = expected
    result["ic"] = ic
    result["ic025"] = ic025
    result["ic975"] = ic975
    return result


def calculate_disproportionality(
    df: pd.DataFrame,
    a_col: ColumnName,
    b_col: ColumnName,
    c_col: ColumnName,
    d_col: ColumnName,
    correction: float = 0.5,
    shrinkage: float = 0.5,
    metrics: Iterable[MetricName] = ("ror", "prr", "ic"),
    add_signal_flags: bool = False,
    inplace: bool = False,
) -> pd.DataFrame:
    """Calculate selected disproportionality metrics for a DataFrame.

    Parameters
    ----------
    df:
        DataFrame containing A, B, C, and D contingency table columns.
    a_col, b_col, c_col, d_col:
        Column names corresponding to A, B, C, and D.
    correction:
        Continuity correction used for ROR and PRR.
    shrinkage:
        Shrinkage value used for IC, IC025, and IC975.
    metrics:
        Metrics to calculate. Valid values are ``"ror"``, ``"prr"``, and
        ``"ic"``.
    add_signal_flags:
        If True, add boolean screening flags for calculated metrics.
    inplace:
        If True, modify the input DataFrame and return it. If False, return a
        copy with the calculated columns.
    """
    selected_metrics = tuple(metrics)
    invalid_metrics = set(selected_metrics) - VALID_METRICS
    if invalid_metrics:
        invalid = ", ".join(sorted(str(metric) for metric in invalid_metrics))
        valid = ", ".join(sorted(VALID_METRICS))
        raise ValueError(f"Invalid metric(s): {invalid}. Valid metrics: {valid}.")

    result = _prepare_result_frame(df, inplace)
    get_contingency_arrays(result, a_col, b_col, c_col, d_col)

    if "ror" in selected_metrics:
        calculate_ror(
            result,
            a_col,
            b_col,
            c_col,
            d_col,
            correction=correction,
            inplace=True,
        )

    if "prr" in selected_metrics:
        calculate_prr(
            result,
            a_col,
            b_col,
            c_col,
            d_col,
            correction=correction,
            inplace=True,
        )

    if "ic" in selected_metrics:
        calculate_ic(
            result,
            a_col,
            b_col,
            c_col,
            d_col,
            shrinkage=shrinkage,
            inplace=True,
        )

    if add_signal_flags:
        if "ror" in selected_metrics:
            result["signal_ror"] = result["ror_lower_95"] > 1
        if "prr" in selected_metrics:
            result["signal_prr"] = result["prr_lower_95"] > 1
        if "ic" in selected_metrics:
            result["signal_ic"] = result["ic025"] > 0

    return result
