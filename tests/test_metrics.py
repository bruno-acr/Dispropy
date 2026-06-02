import numpy as np
import pandas as pd
import pytest

import dispropy as disp
from dispropy import (
    calculate_disproportionality,
    calculate_ic,
    calculate_prr,
    calculate_ror,
)


def example_df():
    return pd.DataFrame({"A": [10], "B": [90], "C": [20], "D": [880]})


def test_calculate_ror_with_known_values():
    result = calculate_ror(example_df(), "A", "B", "C", "D", correction=0)
    expected_ror = (10 * 880) / (90 * 20)

    assert result.loc[0, "ror"] == pytest.approx(expected_ror)
    assert result.loc[0, "log_ror"] == pytest.approx(np.log(expected_ror))


def test_calculate_prr_with_known_values():
    result = calculate_prr(example_df(), "A", "B", "C", "D", correction=0)
    expected_prr = (10 / (10 + 90)) / (20 / (20 + 880))

    assert result.loc[0, "prr"] == pytest.approx(expected_prr)
    assert result.loc[0, "log_prr"] == pytest.approx(np.log(expected_prr))


def test_calculate_ic_with_known_values():
    result = calculate_ic(example_df(), "A", "B", "C", "D", shrinkage=0.5)
    expected_count = ((10 + 90) * (10 + 20)) / (10 + 90 + 20 + 880)
    expected_ic = np.log2((10 + 0.5) / (expected_count + 0.5))
    expected_ic025 = (
        expected_ic
        - 3.3 * np.power(10 + 0.5, -0.5)
        - 2.0 * np.power(10 + 0.5, -1.5)
    )

    assert result.loc[0, "expected_count"] == pytest.approx(expected_count)
    assert result.loc[0, "ic"] == pytest.approx(expected_ic)
    assert result.loc[0, "ic025"] == pytest.approx(expected_ic025)


def test_calculate_disproportionality_with_zeros():
    df = pd.DataFrame({"A": [0], "B": [0], "C": [1], "D": [10]})

    result = calculate_disproportionality(
        df,
        "A",
        "B",
        "C",
        "D",
        correction=0.5,
        shrinkage=0.5,
        add_signal_flags=True,
    )

    assert np.isfinite(result.loc[0, "ror"])
    assert np.isfinite(result.loc[0, "prr"])
    assert np.isfinite(result.loc[0, "ic"])
    assert "signal_ic" in result.columns


def test_calculate_disproportionality_rejects_zero_total_count():
    df = pd.DataFrame({"A": [0], "B": [0], "C": [0], "D": [0]})

    with pytest.raises(ValueError, match="A\\+B\\+C\\+D cannot be zero"):
        calculate_disproportionality(df, "A", "B", "C", "D")


def test_inplace_false_does_not_modify_original_dataframe():
    df = example_df()

    result = calculate_disproportionality(df, "A", "B", "C", "D", inplace=False)

    assert "ror" not in df.columns
    assert "ror" in result.columns


def test_inplace_true_modifies_original_dataframe():
    df = example_df()

    result = calculate_disproportionality(df, "A", "B", "C", "D", inplace=True)

    assert result is df
    assert "ror" in df.columns
    assert "prr" in df.columns
    assert "ic" in df.columns


def test_selected_metrics_calculates_only_ror():
    result = calculate_disproportionality(
        example_df(),
        "A",
        "B",
        "C",
        "D",
        metrics=("ror",),
    )

    assert "ror" in result.columns
    assert "prr" not in result.columns
    assert "ic" not in result.columns


def test_invalid_metric_raises_value_error():
    with pytest.raises(ValueError, match="Invalid metric"):
        calculate_disproportionality(
            example_df(),
            "A",
            "B",
            "C",
            "D",
            metrics=("ror", "ebgm"),
        )


def test_short_api_aliases_are_available():
    df = example_df()

    ror_result = disp.ror(df, "A", "B", "C", "D", correction=0)
    prr_result = disp.prr(df, "A", "B", "C", "D", correction=0)
    ic_result = disp.ic(df, "A", "B", "C", "D")

    assert ror_result.loc[0, "ror"] == pytest.approx((10 * 880) / (90 * 20))
    assert prr_result.loc[0, "prr"] == pytest.approx(
        (10 / (10 + 90)) / (20 / (20 + 880))
    )
    assert "ic025" in ic_result.columns
