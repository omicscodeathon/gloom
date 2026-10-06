"""Unit tests for the ranking metrics and bootstrap confidence intervals (hand-computed examples)."""
import numpy as np
import pytest

from gloom_utils import (auprc_trapezoid, auroc_score, average_precision_score_np,
                         bootstrap_ranking_metrics, compute_ranking_metrics,
                         paired_bootstrap_difference, stratified_bootstrap_indices,
                         tie_broken_score)

# Ranking: scores 0.9 (pos), 0.8 (neg), 0.7 (pos), 0.6 (neg), 0.5 (neg)
Y = np.array([1, 0, 1, 0, 0])
S = np.array([0.9, 0.8, 0.7, 0.6, 0.5])


def test_auroc_by_hand():
    # positive pairs: 0.9 beats 3 negatives, 0.7 beats 2 of 3 -> 5 / 6
    assert auroc_score(Y, S) == pytest.approx(5 / 6)


def test_auroc_ties_count_half():
    assert auroc_score([1, 0], [0.5, 0.5]) == pytest.approx(0.5)


def test_average_precision_by_hand():
    # hits at ranks 1 and 3: precision 1 and 2/3 -> mean = 5/6
    assert average_precision_score_np(Y, S) == pytest.approx(5 / 6)


def test_trapezoidal_auprc_by_hand():
    # PR points (recall, precision): (0,1) (.5,1) (.5,.5) (1,2/3) (1,.5) (1,.4)
    # area = .5*1 + 0 + .5*(0.5 + 2/3)/2 = 0.5 + 0.2916667
    assert auprc_trapezoid(Y, S) == pytest.approx(0.5 + 0.5 * (0.5 + 2 / 3) / 2)


def test_perfect_ranking():
    y = [1, 1, 0, 0]
    s = [0.9, 0.8, 0.2, 0.1]
    m = compute_ranking_metrics(y, s, ks=(2,))
    assert m["auroc"] == pytest.approx(1.0)
    assert m["average_precision"] == pytest.approx(1.0)
    assert m["auprc"] == pytest.approx(1.0)
    assert m["precision_at_2"] == pytest.approx(1.0)
    assert m["recall_at_2"] == pytest.approx(1.0)
    assert m["ef_at_2"] == pytest.approx(2.0)                         # (2/2) / (2/4)


def test_precision_recall_enrichment_at_k_by_hand():
    m = compute_ranking_metrics(Y, S, ks=(1, 2, 3))
    assert m["precision_at_1"] == pytest.approx(1.0)
    assert m["recall_at_1"] == pytest.approx(0.5)
    assert m["ef_at_1"] == pytest.approx(1.0 / 0.4)                   # base rate = 2/5
    assert m["precision_at_2"] == pytest.approx(0.5)
    assert m["recall_at_2"] == pytest.approx(0.5)
    assert m["ef_at_2"] == pytest.approx(0.5 / 0.4)
    assert m["precision_at_3"] == pytest.approx(2 / 3)
    assert m["recall_at_3"] == pytest.approx(1.0)


def test_k_larger_than_n_is_capped():
    m = compute_ranking_metrics(Y, S, ks=(10,))
    assert m["precision_at_10"] == pytest.approx(2 / 5)
    assert m["recall_at_10"] == pytest.approx(1.0)


def test_metrics_do_not_depend_on_input_order():
    perm = np.array([3, 0, 4, 2, 1])
    a = compute_ranking_metrics(Y, S, ks=(2,))
    b = compute_ranking_metrics(Y[perm], S[perm], ks=(2,))
    for key in a:
        assert a[key] == pytest.approx(b[key])


def test_tie_broken_score_orders_by_primary_then_secondary():
    score = tie_broken_score([1, 1, 2], [5, 3, 0])
    assert score.tolist() == [1.0, 0.0, 2.0]


def test_stratified_bootstrap_keeps_class_counts():
    y = np.array([1] * 5 + [0] * 20)
    idx = stratified_bootstrap_indices(y, np.random.default_rng(0))
    assert len(idx) == 25
    assert y[idx].sum() == 5


def _synthetic(n_pos=30, n_neg=170, shift=1.0, seed=0):
    rng = np.random.default_rng(seed)
    y = np.array([1] * n_pos + [0] * n_neg)
    s = np.concatenate([rng.normal(shift, 1, n_pos), rng.normal(0, 1, n_neg)])
    return y, s


def test_bootstrap_ci_shape_columns_and_bounds():
    y, s = _synthetic()
    res = bootstrap_ranking_metrics(y, s, ks=(10, 50), n_boot=60, seed=1)
    assert list(res.columns) == ["metric", "estimate", "ci_low", "ci_high", "n_boot"]
    expected = {"auroc", "auprc", "average_precision",
                "precision_at_10", "recall_at_10", "ef_at_10",
                "precision_at_50", "recall_at_50", "ef_at_50"}
    assert set(res["metric"]) == expected
    assert len(res) == len(expected)
    assert (res["ci_low"] <= res["ci_high"]).all()
    assert np.isfinite(res[["estimate", "ci_low", "ci_high"]].to_numpy()).all()
    assert (res["n_boot"] == 60).all()


def test_bootstrap_is_reproducible_with_fixed_seed():
    y, s = _synthetic()
    a = bootstrap_ranking_metrics(y, s, ks=(10,), n_boot=40, seed=5)
    b = bootstrap_ranking_metrics(y, s, ks=(10,), n_boot=40, seed=5)
    assert a.equals(b)


def test_bootstrap_ci_brackets_the_point_estimate_for_auroc():
    y, s = _synthetic(shift=1.5)
    res = bootstrap_ranking_metrics(y, s, ks=(10,), n_boot=200, seed=2).set_index("metric")
    row = res.loc["auroc"]
    assert row["ci_low"] <= row["estimate"] <= row["ci_high"]
    assert row["ci_high"] - row["ci_low"] > 0


def test_paired_bootstrap_identical_scorers_have_zero_difference():
    y, s = _synthetic()
    res = paired_bootstrap_difference(y, s, s, ks=(10,), n_boot=50, seed=3)
    assert (res["diff"] == 0).all()
    assert (res["ci_low"] == 0).all() and (res["ci_high"] == 0).all()
    assert (res["p_value"] == 1.0).all()


def test_paired_bootstrap_detects_a_clearly_better_scorer():
    rng = np.random.default_rng(4)
    y = np.array([1] * 20 + [0] * 80)
    s_good = y + rng.normal(0, 0.01, size=100)                        # perfect separation
    s_random = rng.random(100)
    res = paired_bootstrap_difference(y, s_good, s_random, ks=(10,), n_boot=200, seed=4)
    row = res.set_index("metric").loc["auroc"]
    assert row["diff"] > 0.3
    assert row["ci_low"] > 0
    assert row["p_value"] <= 0.05
