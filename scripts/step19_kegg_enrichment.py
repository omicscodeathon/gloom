"""
step19_kegg_enrichment.py - KEGG Pathway Enrichment of non-LCGene candidates
=============================================================================
Over-representation analysis of the non-LCGene candidate genes (step14) against KEGG pathways.

PRIMARY OUTPUT (complete, unfiltered)
  * Explicit background = the genes actually analysed by GLOOM (the analysis universe,
    gene_labels.csv) — NOT all genes of KEGG.  Pathway members outside the universe are ignored.
  * One-sided hypergeometric test for EVERY KEGG term with at least one background gene,
    Benjamini-Hochberg correction across all tested terms.  No filtering on overlap size or on
    keywords.
  * Columns: term, n_query, n_term_in_background, n_overlap, p_raw, p_adj (+ n_term_total,
    n_background, fold_enrichment, genes).  For compatibility with the dashboards the legacy
    aliases pathway / pvalue / padj / overlap_genes / pathway_size / odds_ratio are also written.
  * Files: kegg_all_candidates.csv (all high-confidence candidates), kegg_upregulated.csv,
    kegg_downregulated.csv, kegg_summary.csv.

SECONDARY OUTPUT (interpretive only, config.ENRICHMENT_THEMATIC_SUBSET)
  kegg_lung_cancer_subset.csv lists the rows of the primary tables whose pathway name contains a
  predefined lung- or cancer-related keyword (LUNG_RELEVANT_PATTERNS).  It is a reading aid
  ONLY: it is not used to determine statistical significance (p_adj is taken from the complete
  BH correction, not recomputed on the subset) nor to claim biological validity.

Gene sets: config.ENRICHMENT_GENESETS_FILE (local GMT) or, when None, the Enrichr library
config.ENRICHMENT_LIBRARY downloaded with gseapy.get_library (needs gseapy + internet).
"""

import logging
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config
from gloom_utils import hypergeom_enrichment, read_gmt
config.create_output_dirs()

logging.basicConfig(
    level=getattr(logging, config.LOG_LEVEL),
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.FileHandler(config.LOG_FILE, encoding="utf-8"), logging.StreamHandler(sys.stdout)],
)
log = logging.getLogger(__name__)

# ── Parameters ────────────────────────────────────────────────────────────────
GENE_SET_LIBRARY = getattr(config, "ENRICHMENT_LIBRARY", "KEGG_2021_Human")
PROB_THRESHOLD   = 0.80    # min predicted_prob to include
TOP_N_PLOT       = 20      # max pathways shown per figure
PADJ_CUTOFF      = 0.05    # significance threshold (display / summary only)
MIN_GENE_OVERLAP = 3       # min overlapping genes for a pathway to be PLOTTED / counted as hit
MIN_QUERY_GENES  = 3       # smaller query sets are not tested

# ── Lung/cancer keyword list — INTERPRETIVE DISPLAY ONLY ─────────────────────
# Substring patterns matched (case-insensitive) against KEGG pathway names to build the
# secondary file kegg_lung_cancer_subset.csv.  The list is NOT used to decide significance and
# is NOT applied to the primary result tables.
LUNG_RELEVANT_PATTERNS = [
    # ── 1. Direct lung cancer pathways ──────────────────────────────────────
    "non-small cell lung",
    "small cell lung",
    "lung cancer",

    # ── 2. Core LUAD oncogenic signaling ────────────────────────────────────
    "erbb signaling",           # EGFR family
    "mapk signaling",           # KRAS/BRAF downstream
    "pi3k-akt signaling",
    "mtor signaling",
    "ras signaling",
    "jak-stat signaling",
    "vegf signaling",
    "hippo signaling",
    "wnt signaling",
    "tgf-beta signaling",
    "notch signaling",
    "hedgehog signaling",
    "foxo signaling",
    "ampk signaling",           # STK11/LKB1 downstream

    # ── 3. Cancer hallmark processes ─────────────────────────────────────────
    "pathways in cancer",
    "p53 signaling",
    "cell cycle",
    "apoptosis",
    "dna replication",
    "nucleotide excision repair",
    "mismatch repair",
    "homologous recombination",
    "base excision repair",
    "ubiquitin mediated proteolysis",
    "senescence",
    "proteoglycans in cancer",
    "micrornas in cancer",
    "transcriptional misregulation in cancer",
    "central carbon metabolism in cancer",

    # ── 4. Tumor microenvironment & invasion ─────────────────────────────────
    "ecm-receptor interaction",
    "focal adhesion",
    "adherens junction",
    "tight junction",
    "regulation of actin cytoskeleton",
    "cell adhesion molecules",
    "cytokine-cytokine receptor interaction",
    "chemokine signaling",
    "nf-kappa b signaling",
    "tnf signaling",
    "il-17 signaling",
    "pd-l1 expression",
    "natural killer cell",
    "t cell receptor signaling",
    "b cell receptor signaling",
    "toll-like receptor signaling",
    "complement and coagulation",
    "oxidative phosphorylation",
    "glycolysis",
    "hypoxia",
]


def _matches_lung_pattern(pathway_name: str) -> bool:
    """Return True if pathway_name contains any predefined lung/cancer keyword (display only)."""
    name_lower = pathway_name.lower()
    return any(pat in name_lower for pat in LUNG_RELEVANT_PATTERNS)


# ── Gene sets and background ──────────────────────────────────────────────────

def _load_gene_sets() -> dict:
    """KEGG gene sets from a local GMT file or from Enrichr via gseapy.get_library."""
    gmt = getattr(config, "ENRICHMENT_GENESETS_FILE", None)
    if gmt and Path(gmt).exists():
        log.info(f"  Reading gene sets from {gmt}")
        return read_gmt(gmt)
    try:
        import gseapy as gp
    except ImportError:
        log.error("gseapy is not installed. Run: pip install gseapy>=0.10.8 "
                  "(or set config.ENRICHMENT_GENESETS_FILE to a local GMT file).")
        raise
    log.info(f"  Downloading gene-set library {GENE_SET_LIBRARY} (gseapy.get_library) …")
    lib = gp.get_library(name=GENE_SET_LIBRARY, organism="Human")
    return {term: [str(g).upper() for g in genes] for term, genes in lib.items()}


def _analysis_universe() -> list:
    """Genes actually analysed by GLOOM = genes of gene_labels.csv (the analysis universe)."""
    labels = pd.read_csv(config.LABELS_FILE)
    return labels["gene"].astype(str).str.strip().str.upper().tolist()


# ── Core enrichment ───────────────────────────────────────────────────────────

def _odds_ratio(k, n, K, N):
    """Haldane-corrected odds ratio of the 2x2 table (overlap / query-only / term-only / rest)."""
    a = k + 0.5
    b = (n - k) + 0.5
    c = (K - k) + 0.5
    d = (N - K - n + k) + 0.5
    return (a * d) / (b * c)


def _run_enrichment(gene_list: list, label: str, gene_sets: dict, background: list):
    """Complete unfiltered table for one query set. Returns (DataFrame, reason-if-empty)."""
    query = [g for g in {str(x).strip().upper() for x in gene_list}]
    bg = set(background)
    n_in_bg = len([g for g in query if g in bg])
    if n_in_bg < MIN_QUERY_GENES:
        reason = (f"Only {n_in_bg} {label.lower()} genes are in the analysis universe; "
                  f"below the minimum of {MIN_QUERY_GENES} genes needed to test.")
        log.warning(f"  [{label}] {reason}")
        return pd.DataFrame(), reason

    log.info(f"  [{label}] Hypergeometric test: {n_in_bg} query genes, background = {len(bg):,} genes …")
    df = hypergeom_enrichment(
        query, gene_sets, bg,
        min_term_size=int(getattr(config, "ENRICHMENT_MIN_TERM_SIZE", 1)),
    )
    if df.empty:
        reason = f"No KEGG term has a member in the analysis universe for the {label.lower()} set."
        return pd.DataFrame(), reason

    N = int(df["n_background"].iloc[0])
    df["odds_ratio"] = [
        _odds_ratio(k, n, K, N)
        for k, n, K in zip(df["n_overlap"], df["n_query"], df["n_term_in_background"])
    ]
    df["neg_log10_padj"] = -np.log10(df["p_adj"].clip(lower=1e-300))
    # legacy column aliases (dashboards in step17 / output.py read these names)
    df["pathway"] = df["term"]
    df["pvalue"] = df["p_raw"]
    df["padj"] = df["p_adj"]
    df["overlap_genes"] = df["n_overlap"]
    df["pathway_size"] = df["n_term_in_background"]
    df["subset"] = label
    log.info(f"  [{label}] {len(df)} terms tested (all kept); "
             f"{int(((df['p_adj'] < PADJ_CUTOFF) & (df['n_overlap'] >= MIN_GENE_OVERLAP)).sum())} "
             f"with BH-adjusted p < {PADJ_CUTOFF} and overlap >= {MIN_GENE_OVERLAP}")
    return df, ""


def _hits(df: pd.DataFrame) -> pd.DataFrame:
    """Rows reported as significant hits (display / summary rule only; the table itself is not filtered)."""
    return df[(df["p_adj"] < PADJ_CUTOFF) & (df["n_overlap"] >= MIN_GENE_OVERLAP)]


# ── Visualisations (use the complete tables; hits chosen by p_adj, not by keywords) ─────

def _barplot(df: pd.DataFrame, label: str, out_path: Path) -> None:
    sig = _hits(df).head(TOP_N_PLOT).sort_values("neg_log10_padj")
    if sig.empty:
        log.warning(f"  [{label}] No significant pathways — bar-plot skipped.")
        return

    cmap   = matplotlib.colormaps["RdYlGn"]
    colors = cmap(sig["neg_log10_padj"] / sig["neg_log10_padj"].max())

    fig, ax = plt.subplots(figsize=(11, max(4, len(sig) * 0.40)))
    ax.barh(sig["term"], sig["neg_log10_padj"], color=colors, edgecolor="white", height=0.7)
    ax.axvline(-np.log10(PADJ_CUTOFF), color="grey", linestyle="--", lw=0.9,
               label=f"padj = {PADJ_CUTOFF}")
    ax.set_xlabel("-log10(BH-adjusted p-value)", fontsize=10)
    ax.set_title(
        f"KEGG enrichment (unfiltered, background = analysis universe) — {label}\n"
        f"(top {len(sig)} by padj, padj < {PADJ_CUTOFF})",
        fontsize=11,
    )
    ax.legend(fontsize=8)
    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    log.info(f"  [{label}] Bar-plot -> {out_path.name}")


def _dotplot(df: pd.DataFrame, label: str, out_path: Path) -> None:
    sig = _hits(df).head(TOP_N_PLOT).copy()
    if sig.empty:
        log.warning(f"  [{label}] Dot-plot skipped.")
        return

    sig["gene_ratio"] = sig["n_overlap"] / sig["n_term_in_background"]
    sig = sig.sort_values("gene_ratio")

    fig, ax = plt.subplots(figsize=(11, max(4, len(sig) * 0.42)))
    sc = ax.scatter(
        sig["gene_ratio"],
        sig["term"],
        s=sig["n_overlap"] * 14,
        c=sig["neg_log10_padj"],
        cmap="RdYlGn",
        edgecolors="grey",
        linewidths=0.4,
        alpha=0.85,
    )
    plt.colorbar(sc, ax=ax, label="-log10(BH-adjusted p-value)", shrink=0.6)
    ax.set_xlabel("Gene ratio  (overlap / pathway genes in background)", fontsize=10)
    ax.set_title(
        f"KEGG enrichment (unfiltered) — {label}\n"
        f"(dot size = overlapping genes)",
        fontsize=11,
    )
    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    log.info(f"  [{label}] Dot-plot -> {out_path.name}")


def _combined_overview(results: dict) -> None:
    subsets = {k: v for k, v in results.items() if not v.empty}
    if not subsets:
        return

    n   = len(subsets)
    fig, axes = plt.subplots(1, n, figsize=(11 * n, 10), squeeze=False)

    for ax, (label, df) in zip(axes[0], subsets.items()):
        sig = _hits(df).head(15).sort_values("neg_log10_padj")
        if sig.empty:
            ax.text(0.5, 0.5, "No significant\npathways",
                    ha="center", va="center", transform=ax.transAxes, fontsize=13)
            ax.set_title(label, fontsize=14)
            continue
        cmap   = matplotlib.colormaps["RdYlGn"]
        colors = cmap(sig["neg_log10_padj"] / sig["neg_log10_padj"].max())
        ax.barh(sig["term"], sig["neg_log10_padj"],
                color=colors, edgecolor="white", height=0.7)
        ax.axvline(-np.log10(PADJ_CUTOFF), color="grey", linestyle="--", lw=1.2,
                   label=f"padj = {PADJ_CUTOFF}")
        ax.set_xlabel("-log10(padj)", fontsize=13)
        ax.set_title(f"KEGG — {label}", fontsize=14, fontweight="bold")
        ax.tick_params(axis="x", labelsize=12)
        ax.tick_params(axis="y", labelsize=11)
        ax.legend(fontsize=11)

    plt.suptitle(
        f"KEGG Pathway Enrichment of Non-LCGene Candidates (unfiltered)\n"
        f"(prob ≥ {PROB_THRESHOLD}, library: {GENE_SET_LIBRARY}, background = analysis universe)",
        fontsize=16, fontweight="bold", y=1.02,
    )
    plt.tight_layout()
    out = config.FIGURES_DIR / "kegg_overview.png"
    fig.savefig(out, dpi=220, bbox_inches="tight")
    plt.close(fig)
    log.info(f"  Combined overview -> {out.name}")


# ── Main ──────────────────────────────────────────────────────────────────────

def run_kegg_enrichment() -> dict:
    log.info("=" * 60)
    log.info("STEP 19 — KEGG PATHWAY ENRICHMENT (complete, unfiltered)")
    log.info("=" * 60)

    nc_path = config.RESULTS_DIR / "non_lcgene_candidates.csv"
    if not nc_path.exists():
        raise FileNotFoundError(f"non_lcgene_candidates.csv not found: {nc_path}\nRun Step 14 first.")

    nc = pd.read_csv(nc_path, index_col=0)
    log.info(f"  Loaded {len(nc):,} non-LCGene candidates")

    high_conf  = nc[nc["predicted_prob"] >= PROB_THRESHOLD]
    all_genes  = high_conf.index.tolist()
    up_genes   = high_conf[
        (high_conf["is_de_significant"] == True) & (high_conf["direction"] == "up")
    ].index.tolist()
    down_genes = high_conf[
        (high_conf["is_de_significant"] == True) & (high_conf["direction"] == "down")
    ].index.tolist()

    log.info(f"  High-confidence (prob >= {PROB_THRESHOLD}): {len(high_conf):,} genes")
    log.info(f"  Gene lists — all: {len(all_genes)}  up: {len(up_genes)}  down: {len(down_genes)}")

    background = _analysis_universe()
    log.info(f"  Background (analysis universe): {len(background):,} genes")
    gene_sets = _load_gene_sets()
    log.info(f"  Gene sets loaded: {len(gene_sets)} terms")

    analyses = [
        (all_genes,  "All candidates", config.KEGG_ALL_FILE),
        (up_genes,   "Upregulated",    config.KEGG_UP_FILE),
        (down_genes, "Downregulated",  config.KEGG_DOWN_FILE),
    ]

    results = {}
    reasons = {}
    gene_counts = {}
    for gene_list, label, out_csv in analyses:
        gene_counts[label] = len(gene_list)
        df, reason = _run_enrichment(gene_list, label, gene_sets, background)
        if not df.empty:
            df.to_csv(out_csv, index=False)
            log.info(f"  [{label}] Saved COMPLETE table ({len(df)} terms) -> {out_csv.name}")
            slug = label.lower().replace(" ", "_")
            _barplot(df, label, config.FIGURES_DIR / f"kegg_barplot_{slug}.png")
            _dotplot(df, label, config.FIGURES_DIR / f"kegg_dotplot_{slug}.png")
        results[label] = df
        reasons[label] = reason

    # ── Secondary output: thematic lung/cancer subset (interpretive only) ─────────────────────
    if getattr(config, "ENRICHMENT_THEMATIC_SUBSET", True):
        parts = []
        for label, df in results.items():
            if df.empty:
                continue
            sub = df[df["term"].apply(_matches_lung_pattern)].copy()
            sub.insert(0, "query_set", label)
            parts.append(sub)
        subset_path = Path(config.ENRICHMENT_DIR) / "kegg_lung_cancer_subset.csv"
        subset_df = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
        subset_df.to_csv(subset_path, index=False)
        log.info(f"  Thematic lung/cancer subset (interpretive only; p_adj NOT recomputed) "
                 f"-> {subset_path.name}  ({len(subset_df)} rows)")

    # ── Summary table ──────────────────────────────────────────────────────
    rows = []
    for label, df in results.items():
        if df.empty:
            rows.append({
                "subset": label,
                "input_genes": gene_counts[label],
                "status": "empty",
                "n_background": len(background),
                "pathways_tested": 0,
                "significant": 0,
                "top_pathway": "-",
                "top_padj": np.nan,
                "note": reasons[label],
            })
            continue
        sig = _hits(df)
        top = df.iloc[0]
        rows.append({
            "subset": label,
            "input_genes": gene_counts[label],
            "status": "ok" if len(sig) > 0 else "no_significant",
            "n_background": int(top["n_background"]),
            "pathways_tested": len(df),
            "significant": len(sig),
            "top_pathway": top["term"],
            "top_padj": float(f"{top['p_adj']:.6g}"),
            "note": "" if len(sig) > 0 else (
                f"No {label.lower()} pathways passed padj < {PADJ_CUTOFF} with >= {MIN_GENE_OVERLAP} overlapping genes."
            ),
        })

    summary = pd.DataFrame(rows)
    summary.to_csv(config.KEGG_SUMMARY_FILE, index=False)
    log.info(f"\n  Summary -> {config.KEGG_SUMMARY_FILE.name}")
    log.info("\n" + summary.to_string(index=False))

    _combined_overview(results)

    log.info("STEP 19 COMPLETE")
    return results


if __name__ == "__main__":
    run_kegg_enrichment()
