"""
step6_coexpression_network.py
------------------------------
Builds a gene co-expression network from tumor expression data.
Uses chunked Pearson/Spearman correlation with a cutoff threshold.
Outputs: coexpression_network_edges.csv + coexpression_network.graphml
"""
import gc, logging, sys, time
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
    handlers=[logging.FileHandler(config.LOG_FILE, encoding="utf-8"), logging.StreamHandler(sys.stdout)])
log = logging.getLogger(__name__)

def _zscore_matrix(X):
    # NaN-safe: statistics ignore missing values and missing cells are set to the gene mean
    # (z = 0), so a few missing observations no longer turn a whole gene into NaN.
    mean = np.nanmean(X, axis=1, keepdims=True)
    std  = np.nanstd(X, axis=1, keepdims=True, ddof=1)
    std  = np.where((std == 0) | np.isnan(std), 1.0, std)
    Z = (X - mean) / std
    return np.where(np.isnan(Z), 0.0, Z)

def compute_pearson_correlation_chunked(expr_matrix, gene_names, cutoff, chunk_size=500):
    n_genes, n_samples = expr_matrix.shape
    log.info(f"  Matrix: {n_genes} genes x {n_samples} samples | cutoff |r|>={cutoff}")
    Z = _zscore_matrix(expr_matrix.astype(np.float32))
    n_chunks   = int(np.ceil(n_genes / chunk_size))
    idx_a_parts, idx_b_parts, corr_parts = [], [], []
    n_edges    = 0
    denom      = n_samples - 1
    t0 = time.time()
    for chunk_idx in range(n_chunks):
        start = chunk_idx * chunk_size
        end   = min(start + chunk_size, n_genes)
        corr_block = (Z[start:end] @ Z.T) / denom
        # Vectorised edge extraction: upper triangle only (column index > global row index)
        rows, cols = np.nonzero(np.abs(corr_block) >= cutoff)
        keep = cols > (rows + start)
        rows, cols = rows[keep], cols[keep]
        if len(rows):
            idx_a_parts.append(rows + start)
            idx_b_parts.append(cols)
            corr_parts.append(corr_block[rows, cols].astype(np.float64))
            n_edges += len(rows)
        if (chunk_idx + 1) % 10 == 0 or chunk_idx == n_chunks - 1:
            log.info(f"  Chunk {chunk_idx+1}/{n_chunks} edges so far: {n_edges:,} ({time.time()-t0:.1f}s)")
        del corr_block; gc.collect()
    log.info(f"  Total edges: {n_edges:,}")
    if not n_edges:
        return pd.DataFrame(columns=["gene_a","gene_b","correlation"])
    gene_arr = np.array(gene_names)
    idx_a = np.concatenate(idx_a_parts)
    idx_b = np.concatenate(idx_b_parts)
    return pd.DataFrame({
        "gene_a": gene_arr[idx_a],
        "gene_b": gene_arr[idx_b],
        "correlation": np.concatenate(corr_parts),
    })

def build_networkx_graph(edges_df):
    G = nx.Graph()
    all_genes = pd.unique(edges_df[["gene_a","gene_b"]].values.ravel())
    G.add_nodes_from(all_genes)
    for _, row in edges_df.iterrows():
        G.add_edge(row["gene_a"], row["gene_b"],
                   weight=float(row["correlation"]),
                   abs_weight=abs(float(row["correlation"])))
    log.info(f"  Graph: {G.number_of_nodes():,} nodes, {G.number_of_edges():,} edges")
    return G

def run_coexpression_network():
    log.info("="*60); log.info("STEP 6 — CO-EXPRESSION NETWORK CONSTRUCTION"); log.info("="*60)
    tumor_expr  = pd.read_csv(config.TUMOR_EXPR_HARMONIZED, index_col=0)
    gene_names  = tumor_expr.index.tolist()
    expr_matrix = tumor_expr.values.astype(np.float32)
    t0 = time.time()
    edges_df = compute_pearson_correlation_chunked(expr_matrix, gene_names,
                                                    cutoff=config.COEXPR_CORRELATION_CUTOFF)
    t_elapsed = time.time() - t0
    G = build_networkx_graph(edges_df)
    # Add isolated nodes
    for g in tumor_expr.index:
        if g not in G:
            G.add_node(g)

    # Plots
    degrees = [d for _, d in G.degree() if d > 0]
    if degrees:
        degree_counts = pd.Series(degrees).value_counts().sort_index()
        fig, axes = plt.subplots(1,2,figsize=(12,4))
        axes[0].bar(degree_counts.index, degree_counts.values, color="steelblue", alpha=0.8, width=0.8)
        axes[0].set_xlabel("Degree"); axes[0].set_ylabel("Count"); axes[0].set_title("Degree Distribution (linear)")
        axes[1].scatter(degree_counts.index, degree_counts.values, color="steelblue", s=12, alpha=0.8)
        axes[1].set_xscale("log"); axes[1].set_yscale("log")
        axes[1].set_xlabel("Degree (log)"); axes[1].set_ylabel("Count (log)"); axes[1].set_title("Degree Distribution (log-log)")
        plt.tight_layout(); fig.savefig(config.FIGURES_DIR/"network_degree_distribution.png", dpi=130, bbox_inches="tight"); plt.close(fig)

    edges_df.to_csv(config.NETWORK_EDGES_FILE, index=False)
    nx.write_graphml(G, str(config.NETWORK_GRAPH_FILE))
    log.info(f"  Runtime: {t_elapsed:.1f}s")
    log.info("STEP 6 COMPLETE")
    return {"edges_df": edges_df, "graph": G}

if __name__ == "__main__":
    r = run_coexpression_network()
    G = r["graph"]
    print(f"Nodes: {G.number_of_nodes():,}  Edges: {G.number_of_edges():,}")
    top_hubs = sorted(G.degree(), key=lambda x: x[1], reverse=True)[:10]
    for gene, deg in top_hubs:
        print(f"  {gene:<15} degree={deg}")
