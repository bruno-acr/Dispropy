import numpy as np
import pandas as pd
import pytest

from dispropy.validation import validate_contingency_columns


def table(**overrides):
    values = {"A": [10], "B": [90], "C": [20], "D": [880]}
    values.update(overrides)
    return pd.DataFrame(values)


def test_requires_dataframe():
    with pytest.raises(TypeError, match="DataFrame"):
        validate_contingency_columns({}, "A", "B", "C", "D")


def test_rejects_missing_column():
    with pytest.raises(ValueError, match="not found"):
        validate_contingency_columns(table().drop(columns="D"), "A", "B", "C", "D")


def test_rejects_nonnumeric_column():
    with pytest.raises(TypeError, match="numeric"):
        validate_contingency_columns(table(A=["10"]), "A", "B", "C", "D")


def test_rejects_negative_values():
    with pytest.raises(ValueError, match="greater than or equal to zero"):
        validate_contingency_columns(table(A=[-1]), "A", "B", "C", "D")


def test_rejects_missing_values():
    with pytest.raises(ValueError, match="missing"):
        validate_contingency_columns(table(A=[np.nan]), "A", "B", "C", "D")


def test_rejects_all_zero_row():
    with pytest.raises(ValueError, match=r"A\+B\+C\+D"):
        validate_contingency_columns(
            pd.DataFrame({"A": [0], "B": [0], "C": [0], "D": [0]}),
            "A",
            "B",
            "C",
            "D",
        )
