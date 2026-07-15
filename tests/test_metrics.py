import numpy as np
import pandas as pd
import pytest
from scipy import special, stats

from dispropy import (
    calculate_disproportionality,
    calculate_ebgm,
    calculate_ic,
    calculate_prr,
    calculate_ror,
)


def table():
    return pd.DataFrame({"A": [10], "B": [90], "C": [20], "D": [880]})


def ebgm_table():
    return pd.DataFrame(
        {
            "A": [0, 1, 2, 5, 10, 20, 40, 3],
            "B": [100, 99, 98, 95, 90, 80, 60, 97],
            "C": [5, 10, 20, 25, 20, 30, 40, 50],
            "D": [895, 890, 880, 875, 880, 870, 860, 850],
        }
    )


def test_ror_known_value():
    result = calculate_ror(table(), "A", "B", "C", "D", correction=0)
    assert result.loc[0, "ror"] == pytest.approx((10 * 880) / (90 * 20))


def test_prr_known_value():
    result = calculate_prr(table(), "A", "B", "C", "D", correction=0)
    expected = (10 / 100) / (20 / 900)
    assert result.loc[0, "prr"] == pytest.approx(expected)


def test_ic_known_value():
    result = calculate_ic(table(), "A", "B", "C", "D")
    expected_count = (100 * 30) / 1000
    expected_ic = np.log2(10.5 / (expected_count + 0.5))
    expected_ic025 = expected_ic - 3.3 * 10.5**-0.5 - 2 * 10.5**-1.5
    assert result.loc[0, "expected_count"] == pytest.approx(expected_count)
    assert result.loc[0, "ic"] == pytest.approx(expected_ic)
    assert result.loc[0, "ic025"] == pytest.approx(expected_ic025)


def test_zero_counts_with_defaults_are_finite():
    df = pd.DataFrame({"A": [0], "B": [0], "C": [0], "D": [1]})
    result = calculate_disproportionality(df, "A", "B", "C", "D")
    assert np.isfinite(result[["ror", "prr", "ic", "ic025", "ic975"]]).all().all()


def test_all_zero_table_is_rejected_for_ic():
    df = pd.DataFrame({"A": [0], "B": [0], "C": [0], "D": [0]})
    with pytest.raises(ValueError, match=r"A\+B\+C\+D"):
        calculate_ic(df, "A", "B", "C", "D")


def test_inplace_behavior():
    original = table()
    copied = calculate_ror(original, "A", "B", "C", "D")
    assert "ror" not in original
    assert "ror" in copied
    returned = calculate_ror(original, "A", "B", "C", "D", inplace=True)
    assert returned is original
    assert "ror" in original


def test_metric_selection_and_invalid_metric():
    result = calculate_disproportionality(table(), "A", "B", "C", "D", metrics=("ror",))
    assert "ror" in result
    assert "prr" not in result
    assert "ic" not in result
    with pytest.raises(ValueError, match="Invalid metric"):
        calculate_disproportionality(table(), "A", "B", "C", "D", metrics=("unknown",))


def test_ebgm_outputs_and_model_metadata():
    result = calculate_ebgm(ebgm_table(), "A", "B", "C", "D")
    assert {"observed_count", "expected_count", "qn", "ebgm", "eb05", "eb95"} <= set(result)
    assert np.isfinite(result[["qn", "ebgm", "eb05", "eb95"]]).all().all()
    assert result["qn"].between(0, 1).all()
    assert (result["eb05"] <= result["ebgm"]).all()
    assert (result["ebgm"] <= result["eb95"]).all()
    assert result.attrs["gps_model"]["converged"] is True


def test_ebgm_in_main_api_and_signal_flag():
    result = calculate_disproportionality(
        ebgm_table(), "A", "B", "C", "D", metrics=("ebgm",), add_signal_flags=True
    )
    assert "signal_ebgm" in result
    assert result["signal_ebgm"].equals((result["eb05"] > 2) & (result["A"] >= 3))


def test_ebgm_requires_multiple_valid_rows():
    with pytest.raises(ValueError, match="at least two rows"):
        calculate_ebgm(table(), "A", "B", "C", "D")


def test_ebgm_n_star_matches_min_observed_count_in_fitting_sample():
    with_zero = ebgm_table()
    result_with_zero = calculate_ebgm(with_zero, "A", "B", "C", "D")
    assert result_with_zero.attrs["gps_model"]["n_star"] == 0

    without_zero = with_zero.copy()
    without_zero["A"] = without_zero["A"] + 1
    result_without_zero = calculate_ebgm(without_zero, "A", "B", "C", "D")
    assert result_without_zero.attrs["gps_model"]["n_star"] == 1


def test_ebgm_log_likelihood_matches_manual_zero_truncated_formula():
    df = ebgm_table().copy()
    df["A"] = df["A"] + 1  # no zero counts, so n_star must be 1
    result = calculate_ebgm(df, "A", "B", "C", "D")
    model = result.attrs["gps_model"]
    assert model["n_star"] == 1

    observed = result["observed_count"].to_numpy()
    expected = result["expected_count"].to_numpy()
    alpha1, beta1 = model["alpha1"], model["beta1"]
    alpha2, beta2, weight = model["alpha2"], model["beta2"], model["weight"]
    p1 = beta1 / (beta1 + expected)
    p2 = beta2 / (beta2 + expected)
    log_f1 = stats.nbinom.logpmf(observed, alpha1, p1) - stats.nbinom.logsf(0, alpha1, p1)
    log_f2 = stats.nbinom.logpmf(observed, alpha2, p2) - stats.nbinom.logsf(0, alpha2, p2)
    components = np.vstack((np.log(weight) + log_f1, np.log1p(-weight) + log_f2))
    manual_log_likelihood = np.sum(special.logsumexp(components, axis=0))

    assert manual_log_likelihood == pytest.approx(model["log_likelihood"], rel=1e-6)
