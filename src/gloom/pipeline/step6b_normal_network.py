"""
step6b_normal_network.py — Normal Tissue Co-expression Network
==============================================================
Builds a gene co-expression network from NORMAL (GTEx lung) expression.
Mirrors the construction in step6 but applied to NORMAL_EXPR_HARMONIZED.

Used by step7b to compute differential co-expression features:
  delta_degree = tumor_degree − normal_degree   (positive → gained connections in tumor)
  etc.

Outputs:
  network/normal_coexpression_network_edges.csv
  network/normal_coexpression_network.graphml
"""
import logging, sys
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

# Reuse the correlation + graph functions from step6 to avoid duplication
from step6_coexpression_network import compute_pearson_correlation_chunked, build_networkx_graph


def run_normal_coexpression_network():
    log.info("=" * 60)
    log.info("STEP 6b — NORMAL TISSUE CO-EXPRESSION NETWORK")
    log.info("=" * 60)

    normal_expr = pd.read_csv(config.NORMAL_EXPR_HARMONIZED, index_col=0)
    gene_names  = normal_expr.index.tolist()
    expr_matrix = normal_expr.values.astype(np.float32)

    log.info(f"  Normal expression: {expr_matrix.shape[0]:,} genes × {expr_matrix.shape[1]:,} samples")
    log.info(f"  Correlation cutoff: |r| >= {config.COEXPR_CORRELATION_CUTOFF}")

    edges_df = compute_pearson_correlation_chunked(
        expr_matrix, gene_names, cutoff=config.COEXPR_CORRELATION_CUTOFF
    )

    G = build_networkx_graph(edges_df)

    # Add isolated nodes (genes with no edges in normal network)
    for g in normal_expr.index:
        if g not in G:
            G.add_node(g)

    edges_df.to_csv(config.NORMAL_NETWORK_EDGES_FILE, index=False)
    nx.write_graphml(G, str(config.NORMAL_NETWORK_GRAPH_FILE))

    log.info(f"  Normal network: {G.number_of_nodes():,} nodes, {G.number_of_edges():,} edges")
    log.info(f"  Edges saved  → {config.NORMAL_NETWORK_EDGES_FILE}")
    log.info(f"  Graph saved  → {config.NORMAL_NETWORK_GRAPH_FILE}")
    log.info("STEP 6b COMPLETE")
    return {"edges_df": edges_df, "graph": G}


if __name__ == "__main__":
    r = run_normal_coexpression_network()
    G = r["graph"]
    print(f"Normal network: {G.number_of_nodes():,} nodes  {G.number_of_edges():,} edges")
    top_hubs = sorted(G.degree(), key=lambda x: x[1], reverse=True)[:10]
    for gene, deg in top_hubs:
        print(f"  {gene:<15} degree={deg}")
