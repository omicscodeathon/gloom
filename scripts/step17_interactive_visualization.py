"""
step17_interactive_visualization.py - Interactive Visualization
Produces interactive HTML visualizations using Plotly:
1. Volcano plot  2. Gene ranking  3. Network  4. Feature importance
5. Combined tabbed dashboard (interactive_dashboard.html)
Requires: pip install plotly
"""
import logging, sys
from pathlib import Path
import networkx as nx
import numpy as np
import pandas as pd

try:
    import plotly.graph_objects as go
    import plotly.express as px
    import plotly.io as pio
    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config
config.create_output_dirs()
logging.basicConfig(level=getattr(logging, config.LOG_LEVEL),
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.FileHandler(config.LOG_FILE), logging.StreamHandler(sys.stdout)])
log = logging.getLogger(__name__)

COLOR_CGC="#E8524A"; COLOR_NOVEL="#3A7DBF"; COLOR_DE_ONLY="#F0882A"; COLOR_NS="#AAAAAA"
COLOR_MAP = {"CGC gene":COLOR_CGC,"Novel candidate":COLOR_NOVEL,"DE significant":COLOR_DE_ONLY,"Not significant":COLOR_NS}

def assign_category(row):
    if row.get("is_cgc_gene",False): return "CGC gene"
    if row.get("novel_candidate",False): return "Novel candidate"
    if row.get("is_de_significant",False): return "DE significant"
    return "Not significant"

def make_interactive_volcano(de_df, ranking, out_path):
    merged = de_df[["log2fc","pvalue_adj","neg_log10_padj","direction","significant"]].copy()
    merged = merged.join(ranking[["predicted_prob","rank","is_cgc_gene","novel_candidate","is_de_significant"]], how="left")
    merged["category"] = merged.apply(assign_category, axis=1)
    merged["gene"]     = merged.index
    merged["neg_log10_padj"] = merged["neg_log10_padj"].replace([np.inf,-np.inf],np.nan).fillna(merged["neg_log10_padj"].quantile(0.999))
    merged["hover"] = merged.apply(lambda r: f"<b>{r['gene']}</b><br>Log2FC: {r['log2fc']:.3f}<br>Adj.P: {r['pvalue_adj']:.2e}<br>Dir: {r['direction']}<br>Prob: {r.get('predicted_prob',0):.4f}<br>Category: {r['category']}", axis=1)
    fig = go.Figure()
    for cat in ["CGC gene","Novel candidate","DE significant","Not significant"]:
        sub = merged[merged["category"]==cat]
        if not len(sub): continue
        fig.add_trace(go.Scatter(x=sub["log2fc"],y=sub["neg_log10_padj"],mode="markers",name=f"{cat} (n={len(sub):,})",
                                  marker=dict(color=COLOR_MAP[cat],size=9 if cat!="Not significant" else 4,opacity=0.80,line=dict(width=0)),
                                  text=sub["hover"],hovertemplate="%{text}<extra></extra>"))
    fig.add_vline(x=config.DE_LOG2FC_THRESHOLD,  line_dash="dash",line_color="black",line_width=1,opacity=0.6)
    fig.add_vline(x=-config.DE_LOG2FC_THRESHOLD, line_dash="dash",line_color="black",line_width=1,opacity=0.6)
    fig.add_hline(y=-np.log10(config.DE_PVALUE_THRESHOLD),line_dash="dot",line_color="black",line_width=1,opacity=0.6)
    fig.update_layout(title="<b>Volcano Plot - LUAD Tumor vs GTEx Normal</b>",xaxis_title="Log2 Fold-Change",yaxis_title="-log10(adj P)",
                      hovermode="closest",plot_bgcolor="white",paper_bgcolor="white",width=900,height=650)
    fig.write_html(str(out_path), include_plotlyjs="cdn")
    log.info(f"  Volcano saved: {out_path}")
    return fig

def make_interactive_ranking(ranking, out_path):
    df = ranking.copy(); df["gene"] = df.index
    df["category"]    = df.apply(assign_category, axis=1)
    df["abs_log2fc"]  = df["log2fc"].abs().fillna(0)
    df["marker_size"] = (df["abs_log2fc"].clip(0,8)*2+4).round(1)
    df["hover"] = df.apply(lambda r: f"<b>{r['gene']}</b><br>Rank: {int(r['rank'])}<br>Prob: {r['predicted_prob']:.4f}<br>Log2FC: {r['log2fc']:.3f}<br>Dir: {r['direction']}<br>Category: {r['category']}", axis=1)
    fig = go.Figure()
    for cat in ["CGC gene","Novel candidate","DE significant","Not significant"]:
        sub = df[df["category"]==cat]
        if not len(sub): continue
        fig.add_trace(go.Scatter(x=sub["rank"],y=sub["predicted_prob"],mode="markers",name=f"{cat} (n={len(sub):,})",
                                  marker=dict(color=COLOR_MAP[cat],size=sub["marker_size"],opacity=0.75,line=dict(width=0)),
                                  text=sub["hover"],hovertemplate="%{text}<extra></extra>"))
    fig.update_layout(title="<b>Gene Ranking - Predicted Probability by Rank</b>",xaxis_title="Rank Position",yaxis_title="Predicted Probability",
                      hovermode="closest",plot_bgcolor="white",paper_bgcolor="white",width=1000,height=600)
    fig.write_html(str(out_path), include_plotlyjs="cdn")
    log.info(f"  Ranking plot saved: {out_path}")
    return fig

def make_interactive_network(node_df, edge_df, out_path, top_n=150):
    top_nodes = node_df[node_df["rank"].notna()].nsmallest(top_n,"rank").index.tolist()
    if not top_nodes: return go.Figure()
    top_set = set(top_nodes)
    sub_edges = edge_df[edge_df["gene_a"].isin(top_set) & edge_df["gene_b"].isin(top_set)]
    SG = nx.Graph(); SG.add_nodes_from(top_nodes)
    for _,row in sub_edges.iterrows(): SG.add_edge(row["gene_a"],row["gene_b"],weight=float(row["abs_weight"]))
    import networkx as nx2
    pos = nx.spring_layout(SG, k=1.8, iterations=80, seed=config.RANDOM_STATE, weight="weight")
    fig = go.Figure()
    for _,row in sub_edges.iterrows():
        ga,gb = row["gene_a"],row["gene_b"]
        if ga not in pos or gb not in pos: continue
        x0,y0=pos[ga]; x1,y1=pos[gb]
        fig.add_trace(go.Scatter(x=[x0,x1,None],y=[y0,y1,None],mode="lines",
                                  line=dict(width=0.8,color="#4A90D9" if row["edge_type"]=="positive" else "#E05A5A"),
                                  opacity=float(np.clip(row["abs_weight"]*0.8,0.05,0.6)),hoverinfo="none",showlegend=False))
    df_plot = node_df.loc[[n for n in top_nodes if n in node_df.index]].copy()
    df_plot["x"] = df_plot.index.map(lambda g: pos.get(g,(0,0))[0])
    df_plot["y"] = df_plot.index.map(lambda g: pos.get(g,(0,0))[1])
    df_plot["category"] = df_plot.apply(assign_category, axis=1)
    df_plot["size"]     = (df_plot["predicted_prob"]*25+6).clip(6,31)
    df_plot["hover"]    = df_plot.apply(lambda r: f"<b>{r.name}</b><br>Rank: {int(r.get('rank',0))}<br>Prob: {r['predicted_prob']:.4f}<br>Dir: {r.get('direction','ns')}<br>Degree: {int(r.get('degree',0))}<br>Cat: {r['category']}", axis=1)
    for cat in ["Not significant","DE significant","Novel candidate","CGC gene"]:
        sub = df_plot[df_plot["category"]==cat]
        if not len(sub): continue
        fig.add_trace(go.Scatter(x=sub["x"],y=sub["y"],mode="markers",name=f"{cat} (n={len(sub)})",
                                  marker=dict(color=COLOR_MAP[cat],size=sub["size"],opacity=0.88,line=dict(width=0.3,color="white")),
                                  hovertemplate=sub["hover"]+"<extra></extra>"))
    fig.update_layout(title=f"<b>Co-expression Sub-network - Top {top_n} Ranked Genes</b>",showlegend=True,hovermode="closest",
                      plot_bgcolor="white",paper_bgcolor="white",xaxis=dict(showgrid=False,zeroline=False,showticklabels=False),
                      yaxis=dict(showgrid=False,zeroline=False,showticklabels=False),width=1000,height=800)
    fig.write_html(str(out_path), include_plotlyjs="cdn")
    log.info(f"  Network plot saved: {out_path}")
    return fig

def make_interactive_feature_importance(imp_table, out_path):
    method_cols = [c for c in imp_table.columns if c not in ("rank","mean_importance","std_importance","feature_group")]
    df = imp_table.reset_index().rename(columns={"index":"feature"}) if "feature" not in imp_table.columns else imp_table.reset_index()
    if "feature" not in df.columns: df = df.rename(columns={"index":"feature"})
    df = df.sort_values("mean_importance", ascending=True)
    groups = df["feature_group"].unique().tolist()
    palette = px.colors.qualitative.Set2
    group_color = {g: palette[i%len(palette)] for i,g in enumerate(groups)}
    df["color"] = df["feature_group"].map(group_color)
    fig = go.Figure()
    for group in groups:
        sub = df[df["feature_group"]==group]
        if not len(sub): continue
        feat_col = "feature" if "feature" in sub.columns else sub.columns[0]
        fig.add_trace(go.Bar(y=sub[feat_col],x=sub["mean_importance"],name=group,orientation="h",
                              marker_color=group_color[group],opacity=0.85,
                              error_x=dict(type="data",array=sub["std_importance"].values,visible=True,color="grey",thickness=1.2,width=3)))
    fig.update_layout(title="<b>Feature Importance - Mean Across Methods</b>",xaxis_title="Normalised Importance",yaxis_title="Feature",
                      barmode="overlay",plot_bgcolor="white",paper_bgcolor="white",height=max(500,len(df)*22),width=950,margin=dict(l=220,r=40,t=80,b=60))
    fig.write_html(str(out_path), include_plotlyjs="cdn")
    log.info(f"  Feature importance saved: {out_path}")
    return fig

def make_combined_dashboard(fig_volcano, fig_ranking, fig_network, fig_importance, out_path):
    def _div(fig, div_id):
        return pio.to_html(fig, full_html=False, include_plotlyjs=False, div_id=div_id, config={"responsive":True})
    vol_div = _div(fig_volcano,"volcano_div"); rank_div = _div(fig_ranking,"ranking_div")
    net_div = _div(fig_network,"network_div"); imp_div  = _div(fig_importance,"importance_div")
    html = f"""<!DOCTYPE html><html lang="en"><head><meta charset="UTF-8">
<title>LUAD ML Pipeline Dashboard</title>
<script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
<style>*{{box-sizing:border-box;margin:0;padding:0}}body{{font-family:Arial,sans-serif;background:#F8F9FA}}
header{{background:linear-gradient(135deg,#1A3A5C,#2D6A9F);color:white;padding:18px 30px}}
header h1{{font-size:1.6em;font-weight:700}}header p{{font-size:.95em;opacity:.85;margin-top:4px}}
.tab-bar{{display:flex;background:#1A3A5C;padding:0 20px}}.tab{{padding:12px 22px;cursor:pointer;color:#AAC4E0;font-size:.95em;border-bottom:3px solid transparent;transition:all .2s;user-select:none}}
.tab:hover{{color:white}}.tab.active{{color:white;border-bottom:3px solid #5BB4F0;background:rgba(255,255,255,.08)}}
.tab-content{{display:none;padding:20px 30px}}.tab-content.active{{display:block}}
.plot-card{{background:white;border-radius:8px;box-shadow:0 2px 8px rgba(0,0,0,.08);padding:16px;margin-bottom:20px}}
.plot-desc{{font-size:.88em;color:#666;margin-bottom:10px;line-height:1.5}}
footer{{text-align:center;padding:16px;font-size:.82em;color:#888;border-top:1px solid #DDD;margin-top:20px}}</style></head><body>
<header><h1>LUAD ML Pipeline - Interactive Results Dashboard</h1>
<p>Machine Learning Pipeline for Lung Adenocarcinoma Candidate Gene Prioritization</p></header>
<div class="tab-bar">
<div class="tab active" onclick="showTab('volcano',event)">Volcano Plot</div>
<div class="tab" onclick="showTab('ranking',event)">Gene Ranking</div>
<div class="tab" onclick="showTab('network',event)">Network</div>
<div class="tab" onclick="showTab('importance',event)">Feature Importance</div></div>
<div id="tab-volcano" class="tab-content active"><div class="plot-card"><div class="plot-desc"><b>Volcano Plot</b> - Differential expression LUAD tumor vs GTEx normal. Hover for gene details.</div>{vol_div}</div></div>
<div id="tab-ranking" class="tab-content"><div class="plot-card"><div class="plot-desc"><b>Gene Ranking</b> - All genes ranked by ML predicted probability. Hover for full annotation.</div>{rank_div}</div></div>
<div id="tab-network" class="tab-content"><div class="plot-card"><div class="plot-desc"><b>Co-expression Network</b> - Top-150 ranked genes. Node size = predicted probability.</div>{net_div}</div></div>
<div id="tab-importance" class="tab-content"><div class="plot-card"><div class="plot-desc"><b>Feature Importance</b> - Normalised importance across methods. Error bars = std.</div>{imp_div}</div></div>
<footer>LUAD ML Pipeline · step17_interactive_visualization.py</footer>
<script>function showTab(name,event){{document.querySelectorAll(".tab-content").forEach(el=>el.classList.remove("active"));document.querySelectorAll(".tab").forEach(el=>el.classList.remove("active"));document.getElementById("tab-"+name).classList.add("active");if(event&&event.target)event.target.classList.add("active");}}</script>
</body></html>"""
    out_path.write_text(html, encoding="utf-8")
    log.info(f"  Combined dashboard saved: {out_path}")

def run_interactive_visualization():
    log.info("="*60); log.info("STEP 17 - INTERACTIVE VISUALIZATION"); log.info("="*60)
    if not PLOTLY_AVAILABLE:
        log.error("Plotly not installed. Run: pip install plotly"); return {}
    ranking   = pd.read_csv(config.GENE_RANKINGS_FILE,    index_col=0)
    de_df     = pd.read_csv(config.DE_RESULTS_FILE,        index_col=0)
    node_df   = pd.read_csv(config.ANNOTATED_NODES_FILE,   index_col=0)
    edge_df   = pd.read_csv(config.ANNOTATED_EDGES_FILE)
    imp_table = pd.read_csv(config.FEATURE_IMPORTANCE_FILE, index_col=0) if config.FEATURE_IMPORTANCE_FILE.exists() else None
    for col in ["is_cgc_gene","novel_candidate","is_de_significant","predicted_prob","rank"]:
        if col not in de_df.columns and col in ranking.columns:
            de_df[col] = ranking[col].reindex(de_df.index)
    log.info("[1] Volcano …")
    fig_v = make_interactive_volcano(de_df, ranking, config.FIGURES_DIR/"interactive_volcano.html")
    log.info("[2] Ranking …")
    fig_r = make_interactive_ranking(ranking, config.FIGURES_DIR/"interactive_ranking.html")
    log.info("[3] Network …")
    fig_n = make_interactive_network(node_df, edge_df, config.FIGURES_DIR/"interactive_network.html")
    log.info("[4] Feature importance …")
    if imp_table is not None and len(imp_table)>0:
        fig_i = make_interactive_feature_importance(imp_table, config.FIGURES_DIR/"interactive_feature_importance.html")
    else:
        fig_i = go.Figure(); fig_i.update_layout(title="Feature importance not available",height=400)
    log.info("[5] Dashboard …")
    make_combined_dashboard(fig_v, fig_r, fig_n, fig_i, config.FIGURES_DIR/"interactive_dashboard.html")
    log.info("STEP 17 COMPLETE")
    return {"fig_volcano":fig_v,"fig_ranking":fig_r,"fig_network":fig_n,"fig_importance":fig_i}

if __name__ == "__main__":
    r = run_interactive_visualization()
    if r:
        for fname in ["interactive_volcano.html","interactive_ranking.html","interactive_network.html","interactive_feature_importance.html","interactive_dashboard.html"]:
            fpath = config.FIGURES_DIR/fname
            if fpath.exists():
                print(f"  {fname:<45} {fpath.stat().st_size//1024:>6,} KB")
