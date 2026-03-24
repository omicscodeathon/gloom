"""
step16_network_export.py - Network Export
Exports annotated network in multiple formats:
GraphML, GML, edge TSV, node TSV, Cytoscape JSON, sub-networks.
Also computes and saves network statistics report.
"""
import json, logging, sys
from pathlib import Path
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from networkx.readwrite import json_graph

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config
config.create_output_dirs()
logging.basicConfig(level=getattr(logging, config.LOG_LEVEL),
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.FileHandler(config.LOG_FILE), logging.StreamHandler(sys.stdout)])
log = logging.getLogger(__name__)
EXPORT_DIR = config.NETWORK_DIR/"exports"

def sanitise_graph(G):
    H = G.copy()
    def _cast(val):
        if val is None: return ""
        if isinstance(val,(np.integer,)): return int(val)
        if isinstance(val,(np.floating,)):
            v=float(val); return 0.0 if (np.isnan(v) or np.isinf(v)) else v
        if isinstance(val,(np.bool_,)): return bool(val)
        if isinstance(val,float): return 0.0 if (np.isnan(val) or np.isinf(val)) else val
        return val
    for node,attrs in H.nodes(data=True):
        for k in list(attrs.keys()): H.nodes[node][k] = _cast(attrs[k])
    for u,v,attrs in H.edges(data=True):
        for k in list(attrs.keys()): H[u][v][k] = _cast(attrs[k])
    return H

def run_network_export():
    log.info("="*60); log.info("STEP 16 — NETWORK EXPORT"); log.info("="*60)
    EXPORT_DIR.mkdir(parents=True, exist_ok=True)
    G = nx.read_graphml(str(config.ANNOTATED_NETWORK_FILE))
    for u,v,data in G.edges(data=True):
        for attr in ("weight","abs_weight"):
            if attr in data: data[attr] = float(data[attr])
    log.info(f"  Graph: {G.number_of_nodes():,} nodes, {G.number_of_edges():,} edges")
    node_df = pd.read_csv(config.ANNOTATED_NODES_FILE, index_col=0)
    G_clean = sanitise_graph(G)

    # Full exports
    nx.write_graphml(G_clean, str(EXPORT_DIR/"network_full.graphml"))
    log.info(f"  GraphML saved: {EXPORT_DIR}/network_full.graphml")
    # GML
    H = nx.relabel_nodes(G_clean, {n: str(n).replace(" ","_") for n in G_clean.nodes()})
    for node,attrs in H.nodes(data=True):
        for k,v in list(attrs.items()):
            if isinstance(v,bool): H.nodes[node][k] = int(v)
    for u,v,attrs in H.edges(data=True):
        for k,val in list(attrs.items()):
            if isinstance(val,bool): H[u][v][k] = int(val)
    nx.write_gml(H, str(EXPORT_DIR/"network_full.gml"))
    # Edge TSV
    rows = [{"gene_a":u,"gene_b":v,"weight":round(float(data.get("weight",0.0)),6),
             "abs_weight":round(float(data.get("abs_weight",0.0)),6),"edge_type":data.get("edge_type","positive")}
            for u,v,data in G_clean.edges(data=True)]
    pd.DataFrame(rows).to_csv(EXPORT_DIR/"network_full_edgelist.tsv", sep="\t", index=False)
    # Node TSV
    node_df.to_csv(EXPORT_DIR/"network_nodes.tsv", sep="\t")
    # Cytoscape JSON
    with open(EXPORT_DIR/"network_cytoscape.json","w") as fh:
        json.dump(json_graph.node_link_data(G_clean), fh, indent=2)
    log.info("  Full network exports complete.")

    # Sub-networks
    subnetwork_sizes = {}
    subdefs = [
        ("top100",     None,               None, 100),
        ("cgc",        "is_cgc_gene",      True, None),
        ("novel",      "novel_candidate",  True, None),
        ("candidates", "predicted_label",  1,    None),
    ]
    for label, mask_col, mask_val, top_n in subdefs:
        if mask_col and mask_col in node_df.columns:
            selected = node_df[node_df[mask_col]==mask_val].index.tolist()
        elif top_n and "rank" in node_df.columns:
            selected = node_df[node_df.index.isin(G_clean.nodes())].nsmallest(top_n,"rank").index.tolist()
        else:
            continue
        selected = [n for n in selected if n in G_clean.nodes()]
        SG = G_clean.subgraph(selected).copy()
        subnetwork_sizes[label] = {"nodes":SG.number_of_nodes(),"edges":SG.number_of_edges()}
        if SG.number_of_nodes()==0: continue
        nx.write_graphml(SG, str(EXPORT_DIR/f"subnetwork_{label}.graphml"))
        rows2 = [{"gene_a":u,"gene_b":v,"weight":round(float(data.get("weight",0.0)),6),"abs_weight":round(float(data.get("abs_weight",0.0)),6),"edge_type":data.get("edge_type","positive")} for u,v,data in SG.edges(data=True)]
        pd.DataFrame(rows2).to_csv(EXPORT_DIR/f"subnetwork_{label}_edgelist.tsv", sep="\t", index=False)
        log.info(f"  [{label}] {SG.number_of_nodes()} nodes, {SG.number_of_edges()} edges")

    # Statistics
    degrees = [d for _,d in G_clean.degree()]
    components = list(nx.connected_components(G_clean))
    largest_cc = max(components, key=len) if components else set()
    try: assortativity = nx.degree_assortativity_coefficient(G_clean)
    except: assortativity = float("nan")
    try: avg_clust = nx.average_clustering(G_clean)
    except: avg_clust = float("nan")
    stats = {"n_nodes":G_clean.number_of_nodes(),"n_edges":G_clean.number_of_edges(),
             "density":round(nx.density(G_clean),8),"n_connected_components":len(components),
             "largest_component_size":len(largest_cc),"avg_clustering_coeff":round(avg_clust,6),
             "avg_degree":round(float(np.mean(degrees)),4),"max_degree":int(np.max(degrees)),
             "degree_assortativity":round(assortativity,6)}
    if "is_cgc_gene" in node_df.columns:   stats["n_cgc_nodes"]   = int(node_df["is_cgc_gene"].sum())
    if "novel_candidate" in node_df.columns: stats["n_novel_nodes"] = int(node_df["novel_candidate"].sum())
    pd.DataFrame(list(stats.items()),columns=["metric","value"]).to_csv(EXPORT_DIR/"network_statistics_report.csv",index=False)
    degree_series = pd.Series(dict(G_clean.degree()),name="degree")
    top_hubs = degree_series.nlargest(20).reset_index(); top_hubs.columns=["gene","degree"]
    lines = ["="*60,"LUAD CO-EXPRESSION NETWORK — STATISTICS REPORT","="*60,""]
    for k,v in stats.items(): lines.append(f"  {k:<35}: {v}")
    lines += ["","TOP 20 HUB GENES","="*60]
    for _,row in top_hubs.iterrows(): lines.append(f"  {row['gene']:<20} degree={int(row['degree'])}")
    (EXPORT_DIR/"network_statistics_report.txt").write_text("\n".join(lines))

    log.info(f"\n  Sub-networks:")
    for label,sizes in subnetwork_sizes.items():
        log.info(f"    {label:<20} nodes={sizes['nodes']:>5,}  edges={sizes['edges']:>7,}")
    log.info("STEP 16 COMPLETE")
    return {"stats":stats,"top_hubs":top_hubs,"subnetwork_sizes":subnetwork_sizes}

if __name__ == "__main__":
    r = run_network_export()
    print("\n--- Network Statistics ---")
    for k,v in r["stats"].items(): print(f"  {k:<35}: {v}")
    print("\n--- Top Hub Genes ---")
    print(r["top_hubs"].to_string(index=False))
