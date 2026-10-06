"""Unit tests for the enrichment helpers: BH correction and the explicit-background hypergeometric test."""
from math import comb

import numpy as np
import pytest

pytest.importorskip("scipy")

from gloom_utils import benjamini_hochberg, hypergeom_enrichment


def test_bh_by_hand():
    p = [0.01, 0.04, 0.03, 0.005]
    # sorted: .005, .01, .03, .04 -> scaled .02, .02, .04, .04 -> monotone already
    assert benjamini_hochberg(p).tolist() == pytest.approx([0.02, 0.04, 0.04, 0.02])


def test_bh_is_monotone_and_capped_at_one():
    p = np.array([0.9, 0.5, 0.2, 0.001, 0.8])
    adj = benjamini_hochberg(p)
    assert (adj <= 1.0).all() and (adj >= p).all()
    order = np.argsort(p)
    assert (np.diff(adj[order]) >= -1e-12).all()


def test_bh_ignores_nan_when_counting_tests():
    adj = benjamini_hochberg([0.01, np.nan, 0.02])
    assert np.isnan(adj[1])
    assert adj[0] == pytest.approx(0.02) and adj[2] == pytest.approx(0.02)   # m = 2, not 3


def test_hypergeometric_p_value_by_hand():
    bg = list("ABCDEFGHIJ")                                   # N = 10
    gene_sets = {"T": ["A", "B", "C", "X", "Y"]}              # X, Y are outside the background
    res = hypergeom_enrichment(["A", "B", "Z"], gene_sets, bg)  # Z is outside the background
    row = res.iloc[0]
    assert row["n_query"] == 2                                # Z ignored
    assert row["n_term_in_background"] == 3                   # X, Y ignored
    assert row["n_term_total"] == 5
    assert row["n_background"] == 10
    assert row["n_overlap"] == 2
    # P(X >= 2) = C(3,2) * C(7,0) / C(10,2) = 3 / 45
    assert row["p_raw"] == pytest.approx(3 / 45)
    assert row["genes"] == "A;B"


def test_larger_background_gives_smaller_p_value():
    members = ["A", "B", "C"]
    small_bg = list("ABCDEFGHIJ")
    big_bg = small_bg + [f"G{i}" for i in range(90)]          # N = 100
    p_small = hypergeom_enrichment(["A", "B"], {"T": members}, small_bg).iloc[0]["p_raw"]
    p_big = hypergeom_enrichment(["A", "B"], {"T": members}, big_bg).iloc[0]["p_raw"]
    assert p_big == pytest.approx(comb(3, 2) / comb(100, 2))
    assert p_big < p_small


def test_perfect_overlap_p_value():
    bg = [f"G{i}" for i in range(20)]
    members = bg[:5]
    res = hypergeom_enrichment(members, {"T": members}, bg)
    assert res.iloc[0]["p_raw"] == pytest.approx(1 / comb(20, 5), rel=1e-6)


def test_all_terms_are_returned_unfiltered_with_bh_over_all_tests():
    bg = [f"G{i}" for i in range(30)]
    gene_sets = {
        "Cell cycle": bg[:5],
        "Random pathway": bg[10:15],
        "No overlap": bg[20:25],
    }
    res = hypergeom_enrichment(bg[:5], gene_sets, bg)
    assert set(res["term"]) == set(gene_sets)                  # nothing filtered out
    assert {"term", "n_query", "n_term_in_background", "n_overlap", "p_raw", "p_adj"} <= set(res.columns)
    no_overlap = res.set_index("term").loc["No overlap"]
    assert no_overlap["n_overlap"] == 0 and no_overlap["p_raw"] == pytest.approx(1.0)
    # BH over the 3 tested terms: adjusted p of the best term = p_raw * 3 / 1
    top = res.iloc[0]
    assert top["term"] == "Cell cycle"
    assert top["p_adj"] == pytest.approx(min(1.0, top["p_raw"] * 3))
    assert res["p_adj"].between(0, 1).all()


def test_terms_without_background_genes_are_skipped():
    bg = ["A", "B", "C", "D"]
    res = hypergeom_enrichment(["A"], {"in": ["A", "B"], "out": ["X", "Y"]}, bg)
    assert list(res["term"]) == ["in"]


def test_case_and_whitespace_are_normalised():
    res = hypergeom_enrichment([" a "], {"T": ["A", "b"]}, ["a", "B", "c"])
    assert res.iloc[0]["n_overlap"] == 1 and res.iloc[0]["n_background"] == 3
