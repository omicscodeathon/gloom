"""
step15_network_annotation.py - Network Annotation
Annotates network nodes with ML results, DE stats, network features.
Assigns colours and sizes for visualization.
Outputs: annotated_network.graphml, annotated_nodes.csv, annotated_edges.csv
"""
import logging, sys
from pathlib import Path
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
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

COLOR_CGC="#E8524A"; COLOR_NOVEL="#3A7DBF"; COLOR_DE_ONLY="#F0882A"; COLOR_NS="#CCCCCC"

def assign_node_color(row):
    if row.get("is_cgc_gene",False): return COLOR_CGC
    if row.get("novel_candidate",False): return COLOR_NOVEL
    if row.get("is_de_significant",False): return COLOR_DE_ONLY
    return COLOR_NS

def assign_node_size(prob, min_size=5.0, max_size=50.0):
    return float(min_size + prob*(max_size-min_size))

def run_network_annotation():
    log.info("="*60); log.info("STEP 15 — NETWORK ANNOTATION"); log.info("="*60)
    G = nx.read_graphml(str(config.NETWORK_GRAPH_FILE))
    for u,v,data in G.edges(data=True):
        for attr in ("weight","abs_weight"):
            if attr in data: data[attr] = float(data[attr])
    log.info(f"  Graph: {G.number_of_nodes():,} nodes, {G.number_of_edges():,} edges")
    ranking      = pd.read_csv(config.GENE_RANKINGS_FILE, index_col=0)
    net_features = pd.read_csv(config.NETWORK_FEATURES_FILE, index_col=0)
    annotation   = pd.read_csv(config.PROCESSED_DIR/"gene_annotation_table.csv", index_col=0)
    graph_nodes  = list(G.nodes())
    n_genes = len(graph_nodes)

    # Build node table
    node_df = pd.DataFrame(index=graph_nodes); node_df.index.name = "gene"
    for col in ["predicted_prob","predicted_label","rank","percentile","novel_candidate"]:
        if col in ranking.columns: node_df[col] = ranking[col].reindex(graph_nodes)
    node_df["predicted_prob"]  = node_df.get("predicted_prob",pd.Series(0.0,index=graph_nodes)).fillna(0.0)
    node_df["rank"]            = node_df.get("rank",pd.Series(n_genes,index=graph_nodes)).fillna(n_genes).astype(int)
    node_df["novel_candidate"] = node_df.get("novel_candidate",pd.Series(False,index=graph_nodes)).fillna(False).astype(bool)
    for col in ["is_cgc_gene","log2fc","pvalue_adj","neg_log10_padj","direction","is_de_significant","mean_tumor","mean_normal"]:
        if col in annotation.columns: node_df[col] = annotation[col].reindex(graph_nodes)
    node_df["is_cgc_gene"]       = node_df.get("is_cgc_gene",pd.Series(False,index=graph_nodes)).fillna(False).astype(bool)
    node_df["is_de_significant"] = node_df.get("is_de_significant",pd.Series(False,index=graph_nodes)).fillna(False).astype(bool)
    node_df["log2fc"]            = node_df.get("log2fc",pd.Series(0.0,index=graph_nodes)).fillna(0.0)
    node_df["direction"]         = node_df.get("direction",pd.Series("ns",index=graph_nodes)).fillna("ns")
    for col in ["degree","weighted_degree","betweenness_centrality","closeness_centrality","eigenvector_centrality","clustering_coefficient","in_largest_component","mean_edge_weight"]:
        if col in net_features.columns: node_df[col] = net_features[col].reindex(graph_nodes).fillna(0.0)
    node_df["node_color"] = node_df.apply(assign_node_color, axis=1)
    node_df["node_size"]  = node_df["predicted_prob"].apply(assign_node_size)

    # Annotate nodes
    def _cast(val):
        if isinstance(val,(bool,np.bool_)): return bool(val)
        if isinstance(val,(np.integer,)): return int(val)
        if isinstance(val,(np.floating,)): return 0.0 if (np.isnan(float(val)) or np.isinf(float(val))) else float(val)
        if isinstance(val,float) and (np.isnan(val) or np.isinf(val)): return 0.0
        return val
    for gene in G.nodes():
        if gene in node_df.index:
            for col, val in node_df.loc[gene].items():
                G.nodes[gene][col] = _cast(val)

    # Edge annotation
    n_pos=0; n_neg=0
    for u,v,data in G.edges(data=True):
        w = float(data.get("weight",0.0))
        G[u][v]["edge_type"] = "positive" if w>=0 else "negative"
        if w>=0: n_pos+=1
        else:    n_neg+=1
    log.info(f"  Positive edges: {n_pos:,}  Negative: {n_neg:,}")

    # Edge table
    rows = []
    for u,v,data in G.edges(data=True):
        ud=G.nodes[u]; vd=G.nodes[v]
        uc = ud.get("is_cgc_gene",False) or ud.get("novel_candidate",False)
        vc = vd.get("is_cgc_gene",False) or vd.get("novel_candidate",False)
        rows.append({"gene_a":u,"gene_b":v,"weight":float(data.get("weight",0.0)),
                     "abs_weight":float(data.get("abs_weight",0.0)),"edge_type":data.get("edge_type","positive"),
                     "gene_a_rank":int(ud.get("rank",0)),"gene_b_rank":int(vd.get("rank",0)),
                     "gene_a_prob":float(ud.get("predicted_prob",0.0)),"gene_b_prob":float(vd.get("predicted_prob",0.0)),
                     "both_candidates":bool(uc and vc)})
    edge_df = pd.DataFrame(rows)

    # Degree vs rank plot
    if "degree" in node_df.columns and "rank" in node_df.columns:
        fig, ax = plt.subplots(figsize=(9,6))
        categories = [
            ("CGC",   node_df["is_cgc_gene"],                                                   COLOR_CGC,     60, 0.85),
            ("Novel", node_df.get("novel_candidate",pd.Series(False,index=node_df.index)),       COLOR_NOVEL,   40, 0.80),
            ("Other", ~node_df.get("is_cgc_gene",pd.Series(False,index=node_df.index)) & ~node_df.get("novel_candidate",pd.Series(False,index=node_df.index)), COLOR_NS, 5, 0.20),
        ]
        for cat_name,(mask,color,size,alpha) in [(c[0],(c[1],c[2],c[3],c[4])) for c in categories]:
            sub = node_df[mask] if isinstance(mask,pd.Series) else node_df
            if len(sub)>0:
                ax.scatter(sub["degree"],sub["rank"],c=color,s=size,alpha=alpha,label=cat_name,linewidths=0)
        ax.set_xlabel("Network Degree"); ax.set_ylabel("ML Rank"); ax.invert_yaxis()
        ax.set_title("Network Degree vs ML Rank"); ax.legend(fontsize=8,markerscale=2)
        plt.tight_layout(); fig.savefig(config.FIGURES_DIR/"annotated_degree_vs_rank.png",dpi=130,bbox_inches="tight"); plt.close(fig)

    nx.write_graphml(G, str(config.ANNOTATED_NETWORK_FILE))
    node_df.to_csv(config.ANNOTATED_NODES_FILE)
    edge_df.to_csv(config.ANNOTATED_EDGES_FILE, index=False)
    log.info(f"  CGC nodes: {node_df['is_cgc_gene'].sum()}  Novel: {node_df.get('novel_candidate',pd.Series()).sum()}")
    log.info("STEP 15 COMPLETE")
    return {"graph":G,"node_df":node_df,"edge_df":edge_df}

if __name__ == "__main__":
    r = run_network_annotation()
    print(r["node_df"].nsmallest(10,"rank")[["rank","predicted_prob","is_cgc_gene","novel_candidate","log2fc","direction"]].round(4))
