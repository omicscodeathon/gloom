"""Unit tests for the remaining pure helpers: cohort check, QC verdicts, pairing, network summary, evidence classes."""
import numpy as np
import pandas as pd
import pytest

from gloom_utils import (adjacency_from_edges, classify_candidate_evidence,
                         detect_cohort_collinearity, evaluate_qc_warnings,
                         infer_patient_ids, network_summary, pair_samples_by_patient,
                         resolve_feature_groups, update_report_section)


# ── cohort design ─────────────────────────────────────────────────────────────────────────────────

def _series(values, prefix="s"):
    return pd.Series(values, index=[f"{prefix}{i}" for i in range(len(values))])


def test_perfect_collinearity_between_cohort_and_group_is_detected():
    groups = _series(["tumor"] * 3 + ["normal"] * 3)
    cohorts = _series(["TCGA"] * 3 + ["GTEx"] * 3)
    res = detect_cohort_collinearity(groups, cohorts)
    assert res["confounded"] is True
    assert res["table"].shape == (2, 2)


def test_mixed_cohorts_are_not_confounded():
    groups = _series(["tumor", "tumor", "normal", "normal"])
    cohorts = _series(["A", "B", "A", "B"])
    assert detect_cohort_collinearity(groups, cohorts)["confounded"] is False


def test_single_shared_cohort_is_not_confounded():
    groups = _series(["tumor", "tumor", "normal", "normal"])
    cohorts = _series(["TCGA"] * 4)
    assert detect_cohort_collinearity(groups, cohorts)["confounded"] is False


def test_partial_overlap_is_not_perfect_collinearity():
    groups = _series(["tumor", "tumor", "tumor", "normal", "normal", "normal"])
    cohorts = _series(["A", "A", "B", "B", "B", "B"])      # cohort B holds both groups
    assert detect_cohort_collinearity(groups, cohorts)["confounded"] is False


# ── QC verdicts ───────────────────────────────────────────────────────────────────────────────────

def test_qc_warnings_flag_the_original_tcga_gtex_pattern():
    msgs = evaluate_qc_warnings(median_log2fc=4.393, fraction_de=9715 / 10986, n_up=9700, n_down=15)
    assert len(msgs) == 3


def test_qc_warnings_are_silent_for_a_plausible_contrast():
    assert evaluate_qc_warnings(median_log2fc=0.05, fraction_de=0.30, n_up=1200, n_down=900) == []


def test_qc_warns_when_there_is_no_downregulation_at_all():
    msgs = evaluate_qc_warnings(0.0, 0.1, n_up=50, n_down=0)
    assert any("ratio" in m for m in msgs)


def test_qc_warns_on_extreme_down_asymmetry():
    msgs = evaluate_qc_warnings(0.0, 0.1, n_up=5, n_down=100)
    assert any("ratio" in m for m in msgs)


# ── report sections ───────────────────────────────────────────────────────────────────────────────

def test_update_report_section_replaces_and_removes(tmp_path):
    path = tmp_path / "qc_cohort_warning.txt"
    update_report_section(path, "A", ["first"])
    update_report_section(path, "B", ["second"])
    update_report_section(path, "A", ["first v2"])
    text = path.read_text(encoding="utf-8")
    assert "first v2" in text and "second" in text and "first\n" not in text
    update_report_section(path, "A", None)
    assert "first v2" not in path.read_text(encoding="utf-8")
    update_report_section(path, "B", None)
    assert not path.exists()


# ── patient pairing ───────────────────────────────────────────────────────────────────────────────

def test_infer_patient_ids_from_tcga_barcodes_and_map():
    ids = infer_patient_ids(["TCGA-AA-1234-01A-11R-A", "TCGA-AA-1234-11A-01R-A", "weird"], None)
    assert ids["TCGA-AA-1234-01A-11R-A"] == "TCGA-AA-1234"
    assert ids["TCGA-AA-1234-11A-01R-A"] == "TCGA-AA-1234"
    assert pd.isna(ids["weird"])
    mapped = infer_patient_ids(["weird"], {"weird": "P9"})
    assert mapped["weird"] == "P9"


def test_pairing_keeps_only_patients_with_both_samples():
    tumor = pd.Series({"T1": "P1", "T2": "P2", "T3": "P3"})
    normal = pd.Series({"N1": "P1", "N3": "P3", "N4": "P4"})
    t_ids, n_ids = pair_samples_by_patient(tumor, normal)
    assert t_ids == ["T1", "T3"] and n_ids == ["N1", "N3"]


def test_pairing_uses_first_sample_when_a_patient_has_duplicates():
    tumor = pd.Series({"T1a": "P1", "T1b": "P1"})
    normal = pd.Series({"N1": "P1"})
    t_ids, n_ids = pair_samples_by_patient(tumor, normal)
    assert t_ids == ["T1a"] and n_ids == ["N1"]


# ── network summary ───────────────────────────────────────────────────────────────────────────────

def test_network_summary_triangle_plus_isolated_node():
    pytest.importorskip("scipy")
    A = adjacency_from_edges([0, 1, 0], [1, 2, 2], n_nodes=4)     # triangle on nodes 0,1,2; node 3 isolated
    s = network_summary(A)
    assert s["n_nodes"] == 4 and s["n_edges"] == 3
    assert s["density"] == pytest.approx(2 * 3 / (4 * 3))
    assert s["mean_degree"] == pytest.approx(1.5)
    assert s["largest_component"] == 3
    assert s["isolated_nodes"] == 1
    assert s["clustering_coefficient"] == pytest.approx(3 / 4)    # three nodes with C = 1, one with 0
    assert s["degree"].tolist() == [2, 2, 2, 0]


def test_adjacency_ignores_duplicates_and_self_loops():
    pytest.importorskip("scipy")
    A = adjacency_from_edges([0, 1, 0, 2], [1, 0, 1, 2], n_nodes=3)
    assert network_summary(A)["n_edges"] == 1


# ── feature groups ────────────────────────────────────────────────────────────────────────────────

def test_feature_groups_resolve_against_existing_columns_only():
    cols = ["tumor_mean", "abs_log2fc", "cohens_d", "degree", "tumor_degree", "mystery"]
    g = resolve_feature_groups(
        cols,
        expression=["tumor_mean", "abs_log2fc", "cohens_d", "normal_mean"],   # normal_mean absent
        network=["degree", "tumor_degree", "delta_degree"],                   # delta_degree absent
        de_stats=["abs_log2fc", "cohens_d"],
    )
    assert g["expression_only"] == ["tumor_mean", "abs_log2fc", "cohens_d"]
    assert g["network_only"] == ["degree", "tumor_degree"]
    assert g["expression_plus_network"] == ["tumor_mean", "abs_log2fc", "cohens_d", "degree", "tumor_degree"]
    assert g["combined_without_de_stats"] == ["tumor_mean", "degree", "tumor_degree"]
    assert g["full"] == cols
    assert g["unassigned"] == ["mystery"]


# ── candidate evidence ────────────────────────────────────────────────────────────────────────────

@pytest.mark.parametrize(
    "args, expected",
    [
        ((0.4, 0.5, 0.6, 20, 100, 1000), "already_LUAD"),
        ((0.0, 0.3, 0.6, 0, 40, 1000), "lung_cancer_unspecified"),
        ((0.0, 0.0, 0.2, 0, 0, 50), "other_cancer"),
        ((0.0, 0.0, 0.0, 1, 0, 2), "indirect_only"),
        ((0.0, 0.0, 0.0, 0, 0, 0), "no_association_found"),
        ((None, None, None, None, None, None), "not_assessed"),
        ((None, None, None, 0, 0, 0), "not_assessed"),        # one source failed: no claim of novelty
        ((None, None, None, 5, 0, 0), "already_LUAD"),         # positive evidence from the other source stands
    ],
)
def test_candidate_evidence_classes(args, expected):
    assert classify_candidate_evidence(*args) == expected
