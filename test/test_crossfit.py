"""Unit tests for stratified cross-fitting (no gene is scored by a model that saw it)."""
import numpy as np
import pytest

from gloom_utils import crossfit_oof_scores, make_stratified_folds


def test_folds_partition_all_genes_and_are_stratified():
    y = np.array([1] * 10 + [0] * 90)
    folds = make_stratified_folds(y, n_splits=5, seed=1)
    assert folds.shape == (100,)
    assert sorted(set(folds.tolist())) == [0, 1, 2, 3, 4]            # every gene has exactly one fold
    for f in range(5):
        assert (folds == f).sum() == 20                               # 100 genes / 5 folds
        assert y[folds == f].sum() == 2                               # 10 positives / 5 folds


def test_folds_with_uneven_class_sizes_still_partition():
    y = np.array([1] * 7 + [0] * 43)
    folds = make_stratified_folds(y, n_splits=5, seed=3)
    assert set(folds.tolist()) == {0, 1, 2, 3, 4}
    pos_per_fold = [int(y[folds == f].sum()) for f in range(5)]
    assert sum(pos_per_fold) == 7 and max(pos_per_fold) - min(pos_per_fold) <= 1


def test_folds_are_reproducible_and_seed_dependent():
    y = np.array([1] * 10 + [0] * 90)
    assert np.array_equal(make_stratified_folds(y, 5, 7), make_stratified_folds(y, 5, 7))
    assert not np.array_equal(make_stratified_folds(y, 5, 7), make_stratified_folds(y, 5, 8))


def test_no_gene_is_scored_by_a_model_that_saw_it():
    n = 60
    rng = np.random.default_rng(0)
    X = np.column_stack([np.arange(n), rng.normal(size=n)])           # column 0 = gene id
    y = np.array([1] * 12 + [0] * 48)
    held_out_sets = []

    def score_fn(X_train, y_train, X_test, seed):
        train_ids = set(X_train[:, 0].astype(int).tolist())
        test_ids = set(X_test[:, 0].astype(int).tolist())
        assert not (train_ids & test_ids), "a held-out gene was part of the training data"
        held_out_sets.append(test_ids)
        return np.full(len(X_test), float(y_train.mean()))

    oof_mean, oof_rep, fold_table = crossfit_oof_scores(X, y, score_fn, n_splits=5, n_repeats=2, seed=0)
    assert len(held_out_sets) == 10                                   # 5 folds x 2 repeats
    for rep in range(2):
        sets = held_out_sets[rep * 5:(rep + 1) * 5]
        assert set().union(*sets) == set(range(n))                    # folds cover every gene ...
        assert sum(len(s) for s in sets) == n                         # ... exactly once
    assert oof_rep.shape == (n, 2)
    assert not np.isnan(oof_rep.to_numpy()).any()
    assert fold_table.shape == (n, 2)


def test_memorising_model_cannot_leak_labels_through_crossfitting():
    """A scorer that returns the training label of any gene it has seen would be perfect
    on resubstitution; under cross-fitting every gene is unseen, so all scores stay neutral."""
    n = 40
    X = np.column_stack([np.arange(n), np.zeros(n)])
    y = np.array([1] * 8 + [0] * 32)

    def memoriser(X_train, y_train, X_test, seed):
        lookup = {int(i): int(l) for i, l in zip(X_train[:, 0], y_train)}
        return np.array([lookup.get(int(i), 0.5) for i in X_test[:, 0]])

    oof_mean, _, _ = crossfit_oof_scores(X, y, memoriser, n_splits=4, n_repeats=3, seed=5)
    assert np.allclose(oof_mean.to_numpy(), 0.5)


def test_output_index_and_mean_over_repeats():
    X = np.arange(30, dtype=float).reshape(-1, 1)
    y = np.array([1] * 6 + [0] * 24)
    genes = [f"G{i}" for i in range(30)]
    calls = {"n": 0}

    def score_fn(X_train, y_train, X_test, seed):
        calls["n"] += 1
        return np.full(len(X_test), float(calls["n"]))

    oof_mean, oof_rep, _ = crossfit_oof_scores(X, y, score_fn, n_splits=3, n_repeats=2, seed=0, index=genes)
    assert list(oof_mean.index) == genes
    assert np.allclose(oof_mean.to_numpy(), oof_rep.mean(axis=1).to_numpy())


def test_pu_bagging_crossfit_recovers_signal_on_synthetic_data():
    pytest.importorskip("sklearn")
    from gloom_utils import auroc_score, pu_bagging_fit_predict

    rng = np.random.default_rng(1)
    n_pos, n_unl = 40, 160
    X = np.vstack([rng.normal(2.0, 1.0, size=(n_pos, 4)), rng.normal(0.0, 1.0, size=(n_unl, 4))])
    y = np.array([1] * n_pos + [0] * n_unl)

    def score_fn(X_train, y_train, X_test, seed):
        return pu_bagging_fit_predict(
            X_train[y_train == 1], X_train[y_train == 0], X_test,
            n_estimators=5, base_n_trees=10, seed=seed, n_jobs=1,
        )

    oof_mean, _, _ = crossfit_oof_scores(X, y, score_fn, n_splits=4, n_repeats=1, seed=0)
    assert auroc_score(y, oof_mean.to_numpy()) > 0.8


def test_pu_bagging_is_reproducible_for_a_fixed_seed():
    pytest.importorskip("sklearn")
    from gloom_utils import pu_bagging_fit_predict

    rng = np.random.default_rng(2)
    X_pos = rng.normal(1.5, 1.0, size=(15, 3))
    X_unl = rng.normal(0.0, 1.0, size=(60, 3))
    X_new = rng.normal(0.0, 1.0, size=(10, 3))
    kwargs = dict(n_estimators=4, base_n_trees=8, n_jobs=1)
    a = pu_bagging_fit_predict(X_pos, X_unl, X_new, seed=11, **kwargs)
    b = pu_bagging_fit_predict(X_pos, X_unl, X_new, seed=11, **kwargs)
    assert np.array_equal(a, b)
