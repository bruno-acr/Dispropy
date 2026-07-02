"""Validation helpers for 2x2 contingency tables."""

from typing import Hashable

import pandas as pd
from pandas.api.types import is_numeric_dtype


def validate_contingency_columns(
    df: pd.DataFrame,
    a_col: Hashable,
    b_col: Hashable,
    c_col: Hashable,
    d_col: Hashable,
) -> None:
    """Validate columns representing A, B, C and D of a 2x2 table.

    Raises clear errors for an invalid DataFrame, missing or duplicate column
    references, nonnumeric values, missing values, and negative counts.
    """
    if not isinstance(df, pd.DataFrame):
        raise TypeError("df must be a pandas DataFrame.")

    columns = (a_col, b_col, c_col, d_col)
    if len(set(columns)) != 4:
        raise ValueError("a_col, b_col, c_col, and d_col must be distinct columns.")

    missing = [column for column in columns if column not in df.columns]
    if missing:
        raise ValueError(f"Contingency column(s) not found: {missing}.")

    nonnumeric = [column for column in columns if not is_numeric_dtype(df[column])]
    if nonnumeric:
        raise TypeError(f"Contingency column(s) must be numeric: {nonnumeric}.")

    values = df.loc[:, list(columns)]
    if values.isna().to_numpy().any():
        raise ValueError("Contingency columns must not contain missing values.")
    if (values < 0).to_numpy().any():
        raise ValueError("Contingency columns must contain values greater than or equal to zero.")
    if (values.sum(axis=1) == 0).any():
        raise ValueError("A+B+C+D must be greater than zero for every row.")


def get_contingency_arrays(
    df: pd.DataFrame,
    a_col: Hashable,
    b_col: Hashable,
    c_col: Hashable,
    d_col: Hashable,
) -> tuple[pd.Series, pd.Series, pd.Series, pd.Series]:
    """Return validated A, B, C and D columns as float Series."""
    validate_contingency_columns(df, a_col, b_col, c_col, d_col)
    return tuple(df[column].astype(float) for column in (a_col, b_col, c_col, d_col))  # type: ignore[return-value]
