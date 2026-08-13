"""
step7b_differential_network_features.py — Differential Co-expression Features
==============================================================================
Computes tumor-vs-normal differential co-expression features by comparing
the tumor network (step6) and the normal network (step6b).

Biological rationale:
  Genes whose network connectivity *increases* in tumor vs. normal tissue are
  candidates for cancer-specific rewiring — a mechanistic signal distinct from
  simple differential expression. This complements the DE features in step4
  and adds genuine network biology that the tumor-only network in step7 misses.

Differential features computed per gene:
  delta_degree          = tumor_degree − normal_degree
  delta_clustering      = tumor_clustering_coeff − normal_clustering_coeff
  delta_mean_edge_weight = tumor_mean_edge_weight − normal_mean_edge_weight
  delta_betweenness     = tumor_betweenness − normal_betweenness
  tumor_specific_degree = edges in tumor AND NOT in normal (new connections)
  normal_specific_degree= edges in normal AND NOT in tumor (lost connections)
  rewiring_ratio        = |tumor_edges − normal_edges| / max(tumor_edges + normal_edges, 1)

Prerequisite: Run step6 AND step6b before this step.
Output: processed/differential_network_features.csv
        (automatically included by step8 if file exists)
"""
import logging, sys, time
from pathlib import Path
import networkx as nx
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config
config.create_output_dirs()
logging.basicConfig(
    level=getattr(logging, config.LOG_LEVEL),
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.FileHandler(config.LOG_FILE, encoding="utf-8"), logging.StreamHandler(sys.stdout)],
)
log = logging.getLogger(__name__)


def _load_graph(path, label):
    log.info(f"  Loading {label} network: {path}")
    if not Path(path).exists():
        raise FileNotFoundError(
            f"{label} network not found: {path}\n"
            f"Run {'step6' if 'normal' not in str(path) else 'step6b'} first."
        )
    G = nx.read_graphml(str(path))
    for u, v, data in G.edges(data=True):
        for attr in ("weight", "abs_weight"):
            if attr in data:
                data[attr] = float(data[attr])
    log.info(f"    {G.number_of_nodes():,} nodes, {G.number_of_edges():,} edges")
    return G


def _degree_dict(G):
    return {node: deg for node, deg in G.degree()}


# Thresholds above which clustering / betweenness are too slow to compute
# on a dense graph and are skipped (set to 0.0 for all genes).
MAX_EDGES_CLUSTERING   = 300_000
MAX_EDGES_BETWEENNESS  = 500_000


def _clustering_dict(G):
    n_edges = G.number_of_edges()
    if n_edges > MAX_EDGES_CLUSTERING:
        log.warning(f"  Skipping clustering coefficient — graph too dense "
                    f"({n_edges:,} edges > {MAX_EDGES_CLUSTERING:,} limit). Set to 0.0.")
        return {node: 0.0 for node in G.nodes()}
    return nx.clustering(G, weight="abs_weight")


def _betweenness_dict(G, k=500):
    n_edges = G.number_of_edges()
    if n_edges > MAX_EDGES_BETWEENNESS:
        log.warning(f"  Skipping betweenness centrality — graph too dense "
                    f"({n_edges:,} edges > {MAX_EDGES_BETWEENNESS:,} limit). Set to 0.0.")
        return {node: 0.0 for node in G.nodes()}
    n = G.number_of_nodes()
    if n > 3000:
        return nx.betweenness_centrality(G, k=k, normalized=True, seed=config.RANDOM_STATE)
    return nx.betweenness_centrality(G, normalized=True)


def _mean_edge_weight_dict(G):
    records = {}
    for node in G.nodes():
        weights = [abs(data.get("weight", 0.0)) for _, _, data in G.edges(node, data=True)]
        records[node] = float(np.mean(weights)) if weights else 0.0
    return records


def run_differential_network_features():
    log.info("=" * 60)
    log.info("STEP 7b — DIFFERENTIAL CO-EXPRESSION FEATURES")
    log.info("=" * 60)

    G_tumor  = _load_graph(config.NETWORK_GRAPH_FILE,        "Tumor")
    G_normal = _load_graph(config.NORMAL_NETWORK_GRAPH_FILE, "Normal")

    # Gene universe = union of both networks
    all_genes = sorted(set(G_tumor.nodes()) | set(G_normal.nodes()))
    log.info(f"  Gene universe (union): {len(all_genes):,}")

    t0 = time.time()

    log.info("  Computing degree features …")
    deg_t = _degree_dict(G_tumor)
    deg_n = _degree_dict(G_normal)

    log.info("  Computing clustering coefficients …")
    clust_t = _clustering_dict(G_tumor)
    clust_n = _clustering_dict(G_normal)

    log.info("  Computing betweenness centrality …")
    between_t = _betweenness_dict(G_tumor)
    between_n = _betweenness_dict(G_normal)

    log.info("  Computing mean edge weight …")
    mew_t = _mean_edge_weight_dict(G_tumor)
    mew_n = _mean_edge_weight_dict(G_normal)

    log.info("  Computing tumor-specific / normal-specific edges …")
    # Per-gene count of edges that appear only in tumor or only in normal
    tumor_specific_deg  = {}
    normal_specific_deg = {}
    rewiring_ratio      = {}
    for gene in all_genes:
        t_nbrs = set(G_tumor.neighbors(gene))  if G_tumor.has_node(gene)  else set()
        n_nbrs = set(G_normal.neighbors(gene)) if G_normal.has_node(gene) else set()
        ts = len(t_nbrs - n_nbrs)
        ns = len(n_nbrs - t_nbrs)
        total = len(t_nbrs | n_nbrs)
        tumor_specific_deg[gene]  = ts
        normal_specific_deg[gene] = ns
        rewiring_ratio[gene]      = (ts + ns) / total if total > 0 else 0.0

    elapsed = time.time() - t0
    log.info(f"  All features computed in {elapsed:.1f}s")

    # ── Assemble DataFrame ──────────────────────────────────────────────────────
    records = []
    for gene in all_genes:
        dt = float(deg_t.get(gene, 0))
        dn = float(deg_n.get(gene, 0))
        ct = float(clust_t.get(gene, 0.0))
        cn = float(clust_n.get(gene, 0.0))
        bt = float(between_t.get(gene, 0.0))
        bn = float(between_n.get(gene, 0.0))
        wt = float(mew_t.get(gene, 0.0))
        wn = float(mew_n.get(gene, 0.0))
        records.append({
            "gene":                    gene,
            "tumor_degree":            dt,
            "normal_degree":           dn,
            "delta_degree":            dt - dn,
            "tumor_clustering":        ct,
            "normal_clustering":       cn,
            "delta_clustering":        ct - cn,
            "tumor_mean_edge_weight":  wt,
            "normal_mean_edge_weight": wn,
            "delta_mean_edge_weight":  wt - wn,
            "tumor_betweenness":       bt,
            "normal_betweenness":      bn,
            "delta_betweenness":       bt - bn,
            "tumor_specific_degree":   float(tumor_specific_deg.get(gene, 0)),
            "normal_specific_degree":  float(normal_specific_deg.get(gene, 0)),
            "rewiring_ratio":          float(rewiring_ratio.get(gene, 0.0)),
        })

    features = pd.DataFrame(records).set_index("gene")
    features.replace([np.inf, -np.inf], 0.0, inplace=True)
    features.fillna(0.0, inplace=True)

    # Summary statistics
    gained = (features["delta_degree"] > 0).sum()
    lost   = (features["delta_degree"] < 0).sum()
    log.info(f"  Genes gaining connections in tumor: {gained:,}  losing: {lost:,}")
    log.info(f"  Mean rewiring_ratio: {features['rewiring_ratio'].mean():.4f}")

    features.to_csv(config.DIFFERENTIAL_NETWORK_FEATURES_FILE)
    log.info(f"  Differential features: {features.shape}")
    log.info(f"  Saved → {config.DIFFERENTIAL_NETWORK_FEATURES_FILE}")
    log.info("STEP 7b COMPLETE")
    return features


if __name__ == "__main__":
    df = run_differential_network_features()
    print(df.describe().round(4))
    print(f"\nShape: {df.shape}")
    print("\nTop 10 genes by delta_degree (most gained connections in tumor):")
    print(df["delta_degree"].nlargest(10).to_string())
