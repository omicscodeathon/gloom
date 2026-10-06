"""Unit tests for the zero / missing / invalid value handling of step 2 (sanitize_expression)."""
import logging

import numpy as np
import pandas as pd

from gloom_utils import sanitize_expression


def test_true_zeros_are_kept_and_log_transform_is_zero():
    df = pd.DataFrame({"s1": [0.0, 5.0], "s2": [3.0, 0.0]}, index=["g1", "g2"])
    clean, report = sanitize_expression(df, "unit")
    assert np.array_equal(clean.to_numpy(), df.to_numpy())          # nothing was replaced
    assert report["n_zero_kept"] == 2
    assert report["n_invalid"] == 0
    logged = np.log2(clean + 1)
    assert logged.loc["g1", "s1"] == 0.0                            # log2(0 + 1) = 0, not NaN
    assert logged.loc["g2", "s2"] == 0.0
    assert np.isclose(logged.loc["g1", "s2"], 2.0)                  # log2(3 + 1)


def test_zeros_enter_group_means():
    """A gene with values [0, 4] must have mean 2 (the old NaN trick gave 4)."""
    df = pd.DataFrame({"s1": [0.0], "s2": [4.0]}, index=["g"])
    clean, _ = sanitize_expression(df)
    assert clean.mean(axis=1).iloc[0] == 2.0


def test_nan_is_kept_as_missing_and_counted():
    df = pd.DataFrame({"s1": [1.0, np.nan], "s2": [np.nan, 2.0]}, index=["g1", "g2"])
    clean, report = sanitize_expression(df)
    assert clean.isna().to_numpy().tolist() == [[False, True], [True, False]]
    assert report["n_nan_kept"] == 2
    assert report["n_invalid"] == 0


def test_negative_and_infinite_values_are_invalid_and_warned(caplog):
    df = pd.DataFrame({"s1": [-1.0, np.inf], "s2": [-np.inf, 2.0]}, index=["g1", "g2"])
    with caplog.at_level(logging.WARNING):
        clean, report = sanitize_expression(df, "unit")
    assert report["n_negative_invalid"] == 1
    assert report["n_inf_invalid"] == 2
    assert report["n_invalid"] == 3
    assert clean.isna().sum().sum() == 3
    assert clean.loc["g2", "s2"] == 2.0                              # the valid value survives
    assert any("INVALID" in rec.getMessage() for rec in caplog.records)


def test_input_is_not_modified():
    df = pd.DataFrame({"s1": [-1.0, 0.0]}, index=["g1", "g2"])
    before = df.copy()
    sanitize_expression(df)
    pd.testing.assert_frame_equal(df, before)
