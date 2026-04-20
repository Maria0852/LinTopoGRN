#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx

# ----------------------
# Config
# ----------------------
input_csv = 'results/main/scLineageGRN/predicted_edges.csv'
outdir = 'results/main/scLineageGRN/'
top_n = 1000  # number of top edges to consider for TF ranking
sub_top = 30  # number of edges per TF to draw
os.makedirs(outdir, exist_ok=True)

# ----------------------
# Load data
# ----------------------
print(f"[load] reading {input_csv} ...")
df = pd.read_csv(input_csv)
if 'prob_T' not in df.columns:
    raise ValueError("predicted_edges.csv must contain column 'prob_T'")
print(f"[info] total predicted edges: {len(df):,}")

# ----------------------
# 1. Histogram of prob_T
# ----------------------
plt.figure(figsize=(6,4))
plt.hist(df['prob_T'], bins=100, color='steelblue', edgecolor='white')
plt.xlabel('Calibrated Probability (prob_T)')
plt.ylabel('Count')
plt.title('Distribution of predicted edge probabilities')
plt.tight_layout()
hist_path = os.path.join(outdir, 'prob_T_hist.png')
plt.savefig(hist_path, dpi=200)
plt.close()
print(f"[plot] saved histogram -> {hist_path}")

# ----------------------
# 2. Top-N edges export
# ----------------------
df_top = df.sort_values('prob_T', ascending=False).head(top_n)
top_path = os.path.join(outdir, f'top{top_n}_edges.csv')
df_top.to_csv(top_path, index=False)
print(f"[export] saved top-{top_n} edges -> {top_path}")

# ----------------------
# 3. Identify most active TFs
# ----------------------
tf_counts = df_top['gene_i'].value_counts()
top_tfs = tf_counts.head(3).index.tolist()
print("\n[top TFs detected]")
for i, tf in enumerate(top_tfs, 1):
    print(f" {i}. {tf} (edges={tf_counts[tf]})")

# ----------------------
# 4. Draw subnetwork for each top TF
# ----------------------
for tf in top_tfs:
    sub_df = df[df['gene_i'] == tf].sort_values('prob_T', ascending=False).head(sub_top)
    if len(sub_df) == 0:
        print(f"[warn] no edges for {tf}")
        continue

    G = nx.DiGraph()
    for _, row in sub_df.iterrows():
        G.add_edge(row['gene_i'], row['gene_j'], weight=row['prob_T'])

    plt.figure(figsize=(8,6))
    pos = nx.spring_layout(G, seed=42, k=0.8)
    edge_weights = [G[u][v]['weight'] for u,v in G.edges()]
    nodesize = [500 if n==tf else 300 for n in G.nodes()]
    nodecolor = ['tomato' if n==tf else 'skyblue' for n in G.nodes()]

    nx.draw_networkx_nodes(G, pos, node_size=nodesize, node_color=nodecolor, alpha=0.9)
    nx.draw_networkx_labels(G, pos, font_size=8)
    nx.draw_networkx_edges(G, pos, arrows=True, edge_color=edge_weights, edge_cmap=plt.cm.Blues, width=2)
    plt.title(f'{tf} predicted regulatory subnetwork (Top {sub_top})')
    plt.axis('off')
    sub_path = os.path.join(outdir, f'{tf}_subnetwork.png')
    plt.tight_layout()
    plt.savefig(sub_path, dpi=300)
    plt.close()
    print(f"[plot] saved {tf} subnetwork -> {sub_path}")

# ----------------------
# 5. Summary statistics
# ----------------------
mean_prob = df['prob_T'].mean()
median_prob = df['prob_T'].median()
max_prob = df['prob_T'].max()
thr_count = (df['prob_T'] > 0.9).sum()

print("\n===== SUMMARY =====")
print(f"Mean prob_T:    {mean_prob:.4f}")
print(f"Median prob_T:  {median_prob:.4f}")
print(f"Max prob_T:     {max_prob:.4f}")
print(f"Edges with prob_T > 0.9: {thr_count:,}")
print(f"Top-{top_n} average prob_T: {df_top['prob_T'].mean():.4f}")
print(f"Outputs saved to: {outdir}")
