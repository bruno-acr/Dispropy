import warnings

import numpy as np
import pandas as pd

from dispropy import GPSFitWarning, calculate_ebgm


def _small_homogeneous_table() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "A": [0, 1, 2, 5, 10, 20, 40, 3],
            "B": [100, 99, 98, 95, 90, 80, 60, 97],
            "C": [5, 10, 20, 25, 20, 30, 40, 50],
            "D": [895, 890, 880, 875, 880, 870, 860, 850],
        }
    )


def _healthy_table() -> pd.DataFrame:
    rng = np.random.default_rng(20260702)
    size = 300
    expected = rng.uniform(0.5, 25.0, size=size)
    component = rng.random(size) < 0.65
    relative_rate = np.empty(size)
    relative_rate[component] = rng.gamma(1.2, 1 / 1.5, component.sum())
    relative_rate[~component] = rng.gamma(4.0, 1 / 1.2, (~component).sum())
    observed = rng.poisson(expected * relative_rate)
    total = np.full(size, 100_000, dtype=int)
    drug_total = rng.integers(500, 5_000, size=size)
    event_total = np.maximum(
        observed + 1, np.rint(expected * total / drug_total).astype(int)
    )
    a = observed.astype(int)
    b = drug_total - a
    c = event_total - a
    d = total - a - b - c
    return pd.DataFrame({"A": a, "B": b, "C": c, "D": d})


def test_warns_and_records_parameters_near_bounds():
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", GPSFitWarning)
        result = calculate_ebgm(
            _small_homogeneous_table(), "A", "B", "C", "D"
        )

    messages = [str(item.message) for item in caught]
    assert any("optimization bound" in message for message in messages)
    model = result.attrs["gps_model"]
    assert model["parameters_near_bounds"] is True
    assert model["near_bound_parameters"]


def test_warns_for_small_valid_sample():
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", GPSFitWarning)
        result = calculate_ebgm(
            _small_homogeneous_table(), "A", "B", "C", "D"
        )

    messages = [str(item.message) for item in caught]
    assert any("fewer than" in message and "50" in message for message in messages)
    assert result.attrs["gps_model"]["valid_pair_count"] == 8


def test_healthy_fit_emits_no_gps_warning():
    with warnings.catch_warnings():
        warnings.simplefilter("error", GPSFitWarning)
        result = calculate_ebgm(_healthy_table(), "A", "B", "C", "D")

    model = result.attrs["gps_model"]
    assert model["parameters_near_bounds"] is False
    assert model["near_bound_parameters"] == []
    assert model["valid_pair_count"] >= 50
