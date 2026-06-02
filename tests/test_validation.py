import pandas as pd
import pytest

from dispropy.validation import validate_contingency_columns


def test_validate_rejects_non_dataframe():
    with pytest.raises(TypeError, match="pandas DataFrame"):
        validate_contingency_columns({"A": [1]}, "A", "B", "C", "D")


def test_validate_rejects_missing_column():
    df = pd.DataFrame({"A": [1], "B": [2], "C": [3]})

    with pytest.raises(ValueError, match="Missing contingency column"):
        validate_contingency_columns(df, "A", "B", "C", "D")


def test_validate_rejects_non_numeric_column():
    df = pd.DataFrame({"A": [1], "B": ["x"], "C": [3], "D": [4]})

    with pytest.raises(TypeError, match="must be numeric"):
        validate_contingency_columns(df, "A", "B", "C", "D")


def test_validate_rejects_negative_values():
    df = pd.DataFrame({"A": [1], "B": [-2], "C": [3], "D": [4]})

    with pytest.raises(ValueError, match="negative values"):
        validate_contingency_columns(df, "A", "B", "C", "D")


def test_validate_rejects_missing_values():
    df = pd.DataFrame({"A": [1], "B": [None], "C": [3], "D": [4]})

    with pytest.raises(ValueError, match="missing values"):
        validate_contingency_columns(df, "A", "B", "C", "D")
