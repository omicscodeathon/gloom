"""
step19_kegg_enrichment.py - KEGG Pathway Enrichment Validation (Lung-Focused)
Validates novel candidate genes by testing enrichment in KEGG pathways,
then filters results to pathways directly relevant to lung cancer biology.

Three gene lists are tested:
  1. All high-confidence novel candidates
  2. Up-regulated DE-significant candidates
  3. Down-regulated DE-significant candidates

Each is filtered to a curated set of lung/cancer-relevant KEGG pathways
covering: direct lung cancer, core LUAD signaling, cancer hallmark processes,
and tumor microenvironment categories.

Outputs: filtered CSV tables + bar-plot and dot-plot figures.
Requires: pip install gseapy>=0.10.8  (uses Enrichr REST API)
"""

import logging
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config
config.create_output_dirs()

logging.basicConfig(
    level=getattr(logging, config.LOG_LEVEL),
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.FileHandler(config.LOG_FILE), logging.StreamHandler(sys.stdout)],
)
log = logging.getLogger(__name__)

# ── Parameters ────────────────────────────────────────────────────────────────
GENE_SET_LIBRARY = "KEGG_2021_Human"
PROB_THRESHOLD   = 0.80    # min predicted_prob to include
TOP_N_PLOT       = 20      # max pathways shown per figure
PADJ_CUTOFF      = 0.05    # significance threshold
MIN_GENE_OVERLAP = 3       # min genes overlapping a pathway to report

# ── Curated lung/cancer-relevant KEGG pathway keywords ───────────────────────
# These substring patterns are matched (case-insensitive) against KEGG pathway
# names. They cover four biological categories relevant to LUAD.
LUNG_RELEVANT_PATTERNS = [
    # ── 1. Direct lung cancer pathways ──────────────────────────────────────
    "non-small cell lung",
    "small cell lung",
    "lung cancer",

    # ── 2. Core LUAD oncogenic signaling ────────────────────────────────────
    # KRAS, EGFR, STK11, KEAP1, and their downstream cascades dominate LUAD
    "erbb signaling",           # EGFR family - most common LUAD driver
    "mapk signaling",           # KRAS/BRAF downstream
    "pi3k-akt signaling",       # mTOR/survival axis
    "mtor signaling",
    "ras signaling",
    "jak-stat signaling",       # immune + cytokine crosstalk
    "vegf signaling",           # angiogenesis in tumors
    "hippo signaling",          # YAP/TAZ - emerging LUAD role
    "wnt signaling",
    "tgf-beta signaling",       # EMT driver
    "notch signaling",
    "hedgehog signaling",
    "foxo signaling",           # apoptosis / cell cycle
    "ampk signaling",           # STK11/LKB1 downstream (major LUAD TSG)

    # ── 3. Cancer hallmark processes ─────────────────────────────────────────
    "pathways in cancer",       # master KEGG cancer map
    "p53 signaling",            # TP53 alterations common in LUAD
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
    "nf-kappa b signaling",     # inflammation / immune evasion
    "tnf signaling",
    "il-17 signaling",
    "pd-l1 expression",
    "natural killer cell",
    "t cell receptor signaling",
    "b cell receptor signaling",
    "toll-like receptor signaling",
    "complement and coagulation",
    "oxidative phosphorylation",   # metabolic reprogramming in tumors
    "glycolysis",
    "hypoxia",
]


def _matches_lung_pattern(pathway_name: str) -> bool:
    """Return True if pathway_name contains any curated lung/cancer keyword."""
    name_lower = pathway_name.lower()
    return any(pat in name_lower for pat in LUNG_RELEVANT_PATTERNS)


# ── Core enrichment ───────────────────────────────────────────────────────────

def _run_enrichr(gene_list: list, label: str) -> pd.DataFrame:
    try:
        import gseapy as gp
    except ImportError:
        log.error("gseapy is not installed. Run: pip install gseapy>=0.10.8")
        raise

    if len(gene_list) < MIN_GENE_OVERLAP:
        log.warning(f"  [{label}] Only {len(gene_list)} genes - skipping.")
        return pd.DataFrame()

    log.info(f"  [{label}] Querying Enrichr ({GENE_SET_LIBRARY}) with {len(gene_list)} genes …")
    try:
        enr = gp.enrichr(
            gene_list=gene_list,
            gene_sets=GENE_SET_LIBRARY,
            organism="human",
            outdir=None,
            verbose=False,
        )
    except Exception as exc:
        log.error(f"  [{label}] Enrichr query failed: {exc}")
        return pd.DataFrame()

    df = enr.results.copy()
    if df.empty:
        log.warning(f"  [{label}] No results returned from Enrichr.")
        return pd.DataFrame()

    # Normalise column names across gseapy versions
    df = df.rename(columns={
        "Term":             "pathway",
        "Overlap":          "overlap",
        "P-value":          "pvalue",
        "Adjusted P-value": "padj",
        "Odds Ratio":       "odds_ratio",
        "Combined Score":   "combined_score",
        "Genes":            "genes",
    })

    # Parse "5/134" overlap strings
    if "overlap" in df.columns and df["overlap"].dtype == object:
        split = df["overlap"].str.split("/", expand=True).astype(int)
        df["overlap_genes"]  = split[0]
        df["pathway_size"]   = split[1]

    df = df[df.get("overlap_genes", pd.Series(0, index=df.index)) >= MIN_GENE_OVERLAP]
    df["neg_log10_padj"] = -np.log10(df["padj"].clip(lower=1e-300))
    df["subset"] = label

    total_before = len(df)

    # ── Apply lung/cancer filter ───────────────────────────────────────────
    df = df[df["pathway"].apply(_matches_lung_pattern)].reset_index(drop=True)
    df = df.sort_values("padj").reset_index(drop=True)

    log.info(f"  [{label}] {total_before} pathways total → {len(df)} after lung/cancer filter")
    return df


# ── Visualisations ────────────────────────────────────────────────────────────

def _barplot(df: pd.DataFrame, label: str, out_path: Path) -> None:
    sig = df[df["padj"] < PADJ_CUTOFF].head(TOP_N_PLOT).sort_values("neg_log10_padj")
    if sig.empty:
        log.warning(f"  [{label}] No significant lung/cancer pathways - bar-plot skipped.")
        return

    cmap   = cm.get_cmap("RdYlGn")
    colors = cmap(sig["neg_log10_padj"] / sig["neg_log10_padj"].max())

    fig, ax = plt.subplots(figsize=(11, max(4, len(sig) * 0.40)))
    ax.barh(sig["pathway"], sig["neg_log10_padj"], color=colors, edgecolor="white", height=0.7)
    ax.axvline(-np.log10(PADJ_CUTOFF), color="grey", linestyle="--", lw=0.9,
               label=f"padj = {PADJ_CUTOFF}")
    ax.set_xlabel("-log₁₀(adjusted p-value)", fontsize=10)
    ax.set_title(
        f"KEGG Lung/Cancer Pathways - {label}\n"
        f"(top {len(sig)} significant, padj < {PADJ_CUTOFF})",
        fontsize=11,
    )
    ax.legend(fontsize=8)
    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    log.info(f"  [{label}] Bar-plot → {out_path.name}")


def _dotplot(df: pd.DataFrame, label: str, out_path: Path) -> None:
    sig = df[df["padj"] < PADJ_CUTOFF].head(TOP_N_PLOT).copy()
    if sig.empty or "overlap_genes" not in sig.columns:
        log.warning(f"  [{label}] Dot-plot skipped.")
        return

    sig["gene_ratio"] = sig["overlap_genes"] / sig["pathway_size"]
    sig = sig.sort_values("gene_ratio")

    fig, ax = plt.subplots(figsize=(11, max(4, len(sig) * 0.42)))
    sc = ax.scatter(
        sig["gene_ratio"],
        sig["pathway"],
        s=sig["overlap_genes"] * 14,
        c=sig["neg_log10_padj"],
        cmap="RdYlGn",
        edgecolors="grey",
        linewidths=0.4,
        alpha=0.85,
    )
    plt.colorbar(sc, ax=ax, label="-log₁₀(adjusted p-value)", shrink=0.6)
    ax.set_xlabel("Gene Ratio  (overlap / pathway size)", fontsize=10)
    ax.set_title(
        f"KEGG Lung/Cancer Pathways - {label}\n"
        f"(dot size = overlapping genes)",
        fontsize=11,
    )
    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    log.info(f"  [{label}] Dot-plot → {out_path.name}")


def _combined_overview(results: dict) -> None:
    subsets = {k: v for k, v in results.items() if not v.empty}
    if not subsets:
        return

    n   = len(subsets)
    fig, axes = plt.subplots(1, n, figsize=(9 * n, 8), squeeze=False)

    for ax, (label, df) in zip(axes[0], subsets.items()):
        sig = df[df["padj"] < PADJ_CUTOFF].head(15).sort_values("neg_log10_padj")
        if sig.empty:
            ax.text(0.5, 0.5, "No significant\npathways",
                    ha="center", va="center", transform=ax.transAxes, fontsize=11)
            ax.set_title(label)
            continue
        cmap   = cm.get_cmap("RdYlGn")
        colors = cmap(sig["neg_log10_padj"] / sig["neg_log10_padj"].max())
        ax.barh(sig["pathway"], sig["neg_log10_padj"],
                color=colors, edgecolor="white", height=0.7)
        ax.axvline(-np.log10(PADJ_CUTOFF), color="grey", linestyle="--", lw=0.9)
        ax.set_xlabel("-log₁₀(padj)", fontsize=9)
        ax.set_title(f"KEGG Lung - {label}", fontsize=10)
        ax.tick_params(axis="y", labelsize=8)

    plt.suptitle(
        f"KEGG Lung/Cancer Pathway Enrichment - Novel LUAD Candidates\n"
        f"(prob ≥ {PROB_THRESHOLD}, library: {GENE_SET_LIBRARY})",
        fontsize=12, y=1.01,
    )
    plt.tight_layout()
    out = config.FIGURES_DIR / "kegg_lung_overview.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    log.info(f"  Combined overview → {out.name}")


# ── Main ──────────────────────────────────────────────────────────────────────

def run_kegg_enrichment() -> dict:
    log.info("=" * 60)
    log.info("STEP 19 - KEGG LUNG/CANCER PATHWAY ENRICHMENT")
    log.info("=" * 60)

    nc_path = config.RESULTS_DIR / "novel_candidates.csv"
    if not nc_path.exists():
        raise FileNotFoundError(f"novel_candidates.csv not found: {nc_path}\nRun Step 14 first.")

    nc = pd.read_csv(nc_path, index_col=0)
    log.info(f"  Loaded {len(nc):,} novel candidates")

    high_conf  = nc[nc["predicted_prob"] >= PROB_THRESHOLD]
    all_genes  = high_conf.index.tolist()
    up_genes   = high_conf[
        (high_conf["is_de_significant"] == True) & (high_conf["direction"] == "up")
    ].index.tolist()
    down_genes = high_conf[
        (high_conf["is_de_significant"] == True) & (high_conf["direction"] == "down")
    ].index.tolist()

    log.info(f"  High-confidence (prob ≥ {PROB_THRESHOLD}): {len(high_conf):,} genes")
    log.info(f"  Gene lists - all: {len(all_genes)}  up: {len(up_genes)}  down: {len(down_genes)}")

    analyses = [
        (all_genes,  "All candidates", config.KEGG_ALL_FILE),
        (up_genes,   "Upregulated",    config.KEGG_UP_FILE),
        (down_genes, "Downregulated",  config.KEGG_DOWN_FILE),
    ]

    results = {}
    for gene_list, label, out_csv in analyses:
        df = _run_enrichr(gene_list, label)
        if not df.empty:
            df.to_csv(out_csv, index=False)
            log.info(f"  [{label}] Saved {len(df)} lung/cancer pathways → {out_csv.name}")
            slug = label.lower().replace(" ", "_")
            _barplot(df, label, config.FIGURES_DIR / f"kegg_barplot_{slug}.png")
            _dotplot(df, label, config.FIGURES_DIR / f"kegg_dotplot_{slug}.png")
        results[label] = df

    # ── Summary table ──────────────────────────────────────────────────────
    rows = []
    for label, df in results.items():
        if df.empty:
            rows.append({"subset": label, "lung_pathways_tested": 0,
                         "significant": 0, "top_pathway": "-", "top_padj": np.nan})
            continue
        sig = df[df["padj"] < PADJ_CUTOFF]
        top = df.iloc[0]
        rows.append({
            "subset":               label,
            "lung_pathways_tested": len(df),
            "significant":          len(sig),
            "top_pathway":          top["pathway"],
            "top_padj":             round(top["padj"], 6),
        })

    summary = pd.DataFrame(rows)
    summary.to_csv(config.KEGG_SUMMARY_FILE, index=False)
    log.info(f"\n  Summary → {config.KEGG_SUMMARY_FILE.name}")
    log.info("\n" + summary.to_string(index=False))

    _combined_overview(results)

    log.info("STEP 19 COMPLETE")
    return results


if __name__ == "__main__":
    run_kegg_enrichment()
