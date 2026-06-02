"""Validation helpers for disproportionality contingency table columns."""

from __future__ import annotations

from typing import Hashable

import pandas as pd
from pandas.api.types import is_numeric_dtype


ColumnName = Hashable


def validate_contingency_columns(
    df: pd.DataFrame,
    a_col: ColumnName,
    b_col: ColumnName,
    c_col: ColumnName,
    d_col: ColumnName,
) -> None:
    """Validate 2x2 contingency table columns in a pandas DataFrame.

    Parameters
    ----------
    df:
        DataFrame containing the contingency table columns.
    a_col, b_col, c_col, d_col:
        Column names corresponding to A, B, C, and D.

    Raises
    ------
    TypeError
        If ``df`` is not a DataFrame or any selected column is not numeric.
    ValueError
        If any selected column is missing, contains missing values, or contains
        negative values.
    """
    if not isinstance(df, pd.DataFrame):
        raise TypeError("df must be a pandas DataFrame.")

    columns = (a_col, b_col, c_col, d_col)
    missing_columns = [column for column in columns if column not in df.columns]
    if missing_columns:
        missing = ", ".join(str(column) for column in missing_columns)
        raise ValueError(f"Missing contingency column(s): {missing}.")

    columns_with_missing = [
        column for column in columns if df[column].isna().any()
    ]
    if columns_with_missing:
        with_missing = ", ".join(str(column) for column in columns_with_missing)
        raise ValueError(
            "Contingency columns cannot contain missing values. "
            f"Column(s) with missing values: {with_missing}."
        )

    non_numeric_columns = [
        column for column in columns if not is_numeric_dtype(df[column])
    ]
    if non_numeric_columns:
        non_numeric = ", ".join(str(column) for column in non_numeric_columns)
        raise TypeError(
            "Contingency columns must be numeric. "
            f"Non-numeric column(s): {non_numeric}."
        )

    columns_with_negative_values = [
        column for column in columns if (df[column] < 0).any()
    ]
    if columns_with_negative_values:
        negative = ", ".join(
            str(column) for column in columns_with_negative_values
        )
        raise ValueError(
            "Contingency columns must contain values greater than or equal to "
            f"zero. Column(s) with negative values: {negative}."
        )


def get_contingency_arrays(
    df: pd.DataFrame,
    a_col: ColumnName,
    b_col: ColumnName,
    c_col: ColumnName,
    d_col: ColumnName,
) -> tuple[pd.Series, pd.Series, pd.Series, pd.Series]:
    """Return A, B, C, and D columns as float Series after validation."""
    validate_contingency_columns(df, a_col, b_col, c_col, d_col)
    return (
        df[a_col].astype(float),
        df[b_col].astype(float),
        df[c_col].astype(float),
        df[d_col].astype(float),
    )
