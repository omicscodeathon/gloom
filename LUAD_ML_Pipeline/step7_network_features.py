"""
step7_network_features.py
--------------------------
Extracts per-gene topological features from co-expression network:
degree, weighted_degree, avg_neighbor_degree, betweenness/closeness/
eigenvector centrality, clustering coefficient, edge weight stats,
component features.
Output: network_features.csv
"""
import logging, sys, time, warnings
from pathlib import Path
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config
config.create_output_dirs()
logging.basicConfig(level=getattr(logging, config.LOG_LEVEL),
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.FileHandler(config.LOG_FILE), logging.StreamHandler(sys.stdout)])
log = logging.getLogger(__name__)
BETWEENNESS_EXACT_LIMIT = 3000
BETWEENNESS_K_SAMPLES   = 500

def extract_degree_features(G):
    log.info("  Extracting degree features …")
    degree   = dict(G.degree())
    wdegree  = dict(G.degree(weight="abs_weight"))
    df = pd.DataFrame({"degree": degree, "weighted_degree": wdegree})
    df.index.name = "gene"
    return df

def extract_avg_neighbor_degree(G):
    log.info("  Extracting average neighbour degree …")
    avg_nd = nx.average_neighbor_degree(G)
    s = pd.Series(avg_nd, name="avg_neighbor_degree")
    s.index.name = "gene"; return s

def extract_betweenness_centrality(G):
    n = G.number_of_nodes()
    log.info(f"  Betweenness centrality (n={n:,}) …")
    if n > BETWEENNESS_EXACT_LIMIT:
        bc = nx.betweenness_centrality(G, k=BETWEENNESS_K_SAMPLES, normalized=True, seed=config.RANDOM_STATE)
    else:
        bc = nx.betweenness_centrality(G, normalized=True)
    s = pd.Series(bc, name="betweenness_centrality"); s.index.name = "gene"; return s

def extract_closeness_centrality(G):
    log.info("  Closeness centrality …")
    cc = nx.closeness_centrality(G)
    s  = pd.Series(cc, name="closeness_centrality"); s.index.name = "gene"; return s

def extract_eigenvector_centrality(G):
    log.info("  Eigenvector centrality …")
    try:
        ec = nx.eigenvector_centrality(G, max_iter=1000, tol=1e-6, weight="abs_weight")
    except nx.PowerIterationFailedConvergence:
        max_deg = max(dict(G.degree()).values()) or 1
        ec = {node: deg/max_deg for node, deg in G.degree()}
    s = pd.Series(ec, name="eigenvector_centrality"); s.index.name = "gene"; return s

def extract_clustering_coefficient(G):
    log.info("  Clustering coefficient …")
    clust = nx.clustering(G, weight="abs_weight")
    s     = pd.Series(clust, name="clustering_coefficient"); s.index.name = "gene"; return s

def extract_edge_weight_stats(G):
    log.info("  Edge weight statistics …")
    records = {}
    for node in G.nodes():
        weights = [abs(data.get("weight",0.0)) for _,_,data in G.edges(node, data=True)]
        if not weights:
            records[node] = {"mean_edge_weight":0.0,"max_edge_weight":0.0,"min_edge_weight":0.0,"std_edge_weight":0.0}
        else:
            w = np.array(weights, dtype=np.float32)
            records[node] = {"mean_edge_weight":float(np.mean(w)),"max_edge_weight":float(np.max(w)),
                             "min_edge_weight":float(np.min(w)),"std_edge_weight":float(np.std(w,ddof=0))}
    df = pd.DataFrame.from_dict(records, orient="index"); df.index.name = "gene"; return df

def extract_component_features(G):
    log.info("  Component features …")
    components  = list(nx.connected_components(G))
    largest_set = max(components, key=len)
    records = {}
    for comp in components:
        size = len(comp)
        for node in comp:
            records[node] = {"in_largest_component": int(node in largest_set), "component_size": size}
    df = pd.DataFrame.from_dict(records, orient="index"); df.index.name = "gene"; return df

def run_network_feature_extraction():
    log.info("="*60); log.info("STEP 7 — NETWORK FEATURE EXTRACTION"); log.info("="*60)
    G = nx.read_graphml(str(config.NETWORK_GRAPH_FILE))
    for u, v, data in G.edges(data=True):
        for attr in ("weight","abs_weight"):
            if attr in data: data[attr] = float(data[attr])
    log.info(f"  Graph: {G.number_of_nodes():,} nodes, {G.number_of_edges():,} edges")
    common_genes = pd.read_csv(config.PROCESSED_DIR/"common_genes.csv").squeeze("columns")
    t0 = time.time()
    f_deg     = extract_degree_features(G)
    f_avnd    = extract_avg_neighbor_degree(G)
    f_between = extract_betweenness_centrality(G)
    f_close   = extract_closeness_centrality(G)
    f_eigen   = extract_eigenvector_centrality(G)
    f_clust   = extract_clustering_coefficient(G)
    f_weights = extract_edge_weight_stats(G)
    f_comp    = extract_component_features(G)
    t_elapsed = time.time() - t0
    log.info(f"  All features extracted in {t_elapsed:.1f}s")
    features = pd.concat([f_deg, f_avnd, f_between, f_close, f_eigen, f_clust, f_weights, f_comp], axis=1)
    features.index.name = "gene"
    features = features.reindex(common_genes.values).fillna(0.0)
    features.replace([np.inf, -np.inf], 0.0, inplace=True)

    # Plots
    feat_cols = features.columns.tolist()
    n_cols = 4; n_rows = int(np.ceil(len(feat_cols)/n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols*3.5, n_rows*2.5))
    axes = np.array(axes).flatten()
    for i, col in enumerate(feat_cols):
        axes[i].hist(features[col].dropna(), bins=60, color="steelblue", alpha=0.8, edgecolor="none")
        axes[i].set_title(col, fontsize=8); axes[i].set_yticks([]); axes[i].tick_params(labelsize=6)
    for j in range(i+1, len(axes)): axes[j].set_visible(False)
    plt.tight_layout(); fig.savefig(config.FIGURES_DIR/"network_feature_distributions.png", dpi=120, bbox_inches="tight"); plt.close(fig)

    features.to_csv(config.NETWORK_FEATURES_FILE)
    log.info(f"  Network features: {features.shape}")
    log.info("STEP 7 COMPLETE")
    return features

if __name__ == "__main__":
    nf = run_network_feature_extraction()
    print(nf.head(5).round(6))
    print(f"Shape: {nf.shape}")
