#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import pandas as pd
import requests
from tqdm import tqdm
import gseapy as gp

# ----------------------
# Config
# ----------------------
input_csv = 'results/main/scLineageGRN/top1000_edges.csv'
outdir = 'results/main/scLineageGRN/validation_report/'
os.makedirs(outdir, exist_ok=True)

# ----------------------
# 1. Load predicted edges
# ----------------------
print(f"[load] reading predicted edges from {input_csv} ...")
df_pred = pd.read_csv(input_csv)
required_cols = {'gene_i', 'gene_j'}
if not required_cols.issubset(df_pred.columns):
    raise ValueError(f"Input must contain columns {required_cols}")

pred_edges = set(zip(df_pred['gene_i'], df_pred['gene_j']))
print(f"[info] total predicted edges: {len(pred_edges):,}")

# ----------------------
# 2. Download TF-target reference sets (TRRUST, RegNetwork)
# ----------------------
print("[download] fetching TF-target reference data ...")
try:
    url_trrust = 'https://www.grnpedia.org/trrust/download/mouse_trrust_rawdata.txt'
    df_trrust = pd.read_csv(url_trrust, sep='\t', header=None, names=['gene_i','gene_j','mode','ref'])
    df_trrust = df_trrust[['gene_i','gene_j']].drop_duplicates()
    print(f"  TRRUST loaded: {len(df_trrust):,} edges")
except Exception as e:
    print(f"  [warn] failed to load TRRUST: {e}")
    df_trrust = pd.DataFrame(columns=['gene_i','gene_j'])

try:
    url_regnet = 'https://regnetworkweb.org/download/human_mus_musculus.txt'
    df_regnet = pd.read_csv(url_regnet, sep='\t', header=None, names=['gene_i','gene_j','type','source'])
    df_regnet = df_regnet[df_regnet['type'].str.contains('TF', na=False)][['gene_i','gene_j']].drop_duplicates()
    print(f"  RegNetwork loaded: {len(df_regnet):,} edges")
except Exception as e:
    print(f"  [warn] failed to load RegNetwork: {e}")
    df_regnet = pd.DataFrame(columns=['gene_i','gene_j'])

# Merge all reference edges
df_ref = pd.concat([df_trrust, df_regnet], ignore_index=True).drop_duplicates()
ref_edges = set(zip(df_ref['gene_i'], df_ref['gene_j']))
print(f"[info] total reference edges combined: {len(ref_edges):,}")

# ----------------------
# 3. Cross-validation
# ----------------------
print("[match] comparing predicted vs reference edges ...")
matched = [edge for edge in tqdm(pred_edges) if edge in ref_edges]
match_rate = len(matched) / len(pred_edges) if pred_edges else 0.0

print(f"[result] matched {len(matched):,} edges ({match_rate*100:.2f}%)")

df_matched = pd.DataFrame(matched, columns=['gene_i','gene_j'])
match_path = os.path.join(outdir, 'matched_edges.csv')
df_matched.to_csv(match_path, index=False)
print(f"[export] matched edges saved -> {match_path}")

# ----------------------
# 4. Global enrichment (all predicted targets)
# ----------------------
print("[enrich] running GO/KEGG enrichment for all predicted targets ...")
pred_genes = list(set(df_pred['gene_j']))

def run_enrich(gene_list, gene_sets, out_csv):
    try:
        enr = gp.enrichr(gene_list=list(gene_list),
                         gene_sets=gene_sets,
                         organism='Mouse',
                         outdir=None, cutoff=0.05)
        df_enr = enr.results[['Gene_set','Term','Overlap','Adjusted P-value','Combined Score']]
        df_enr.to_csv(out_csv, index=False)
        print(f"  [export] enrichment saved -> {out_csv}")
        return df_enr
    except Exception as e:
        print(f"  [warn] enrichment failed: {e}")
        return pd.DataFrame()

global_enrich_path = os.path.join(outdir, 'enrichment_results.csv')
df_enr_global = run_enrich(pred_genes, ['GO_Biological_Process_2023','KEGG_2021_Mouse'], global_enrich_path)

# ----------------------
# 5. Per-TF enrichment for top-3 TFs
# ----------------------
print("[per-TF] identifying top-3 TFs by outgoing edges among top predictions ...")
df_rank = df_pred.sort_values('prob_T', ascending=False).head(1000)
tf_counts = df_rank['gene_i'].value_counts()

top_tfs = tf_counts.head(3).index.tolist()
print("  Top TFs:")
for i, tf in enumerate(top_tfs, 1):
    print(f"   {i}. {tf} (edges={tf_counts[tf]})")

per_tf_summaries = []
for tf in top_tfs:
    tg = df_pred.loc[df_pred['gene_i']==tf, 'gene_j'].unique().tolist()
    # save target list
    tgt_path = os.path.join(outdir, f'{tf}_targets.txt')
    with open(tgt_path, 'w') as f:
        for g in tg: f.write(g+'\n')
    print(f"  [export] {tf} targets -> {tgt_path} (n={len(tg)})")

    # run GO BP enrichment
    tf_enrich_path = os.path.join(outdir, f'{tf}_enrichment.csv')
    df_enr_tf = run_enrich(tg, ['GO_Biological_Process_2023'], tf_enrich_path)

    # summary top terms
    top_terms = []
    if not df_enr_tf.empty:
        for _, row in df_enr_tf.head(5).iterrows():
            top_terms.append(f"{row['Term']} (p_adj={row['Adjusted P-value']:.2e})")
    per_tf_summaries.append((tf, len(tg), top_terms))

# ----------------------
# 6. Summary report
# ----------------------
summary_path = os.path.join(outdir, 'validation_summary.txt')
with open(summary_path, 'w') as f:
    f.write("=== GRN Validation Summary ===\n")
    f.write(f"Predicted edges: {len(pred_edges):,}\n")
    f.write(f"Reference edges: {len(ref_edges):,}\n")
    f.write(f"Matched edges: {len(matched):,} ({match_rate*100:.2f}%)\n\n")
    if not df_enr_global.empty:
        f.write("[Global Enrichment: top 10]\n")
        for _, row in df_enr_global.head(10).iterrows():
            f.write(f"  {row['Gene_set']} :: {row['Term']} (p_adj={row['Adjusted P-value']:.2e})\n")
        f.write("\n")
    f.write("[Top TFs] by outgoing edges among top-1000 predictions\n")
    for tf, n_tg, terms in per_tf_summaries:
        f.write(f"  {tf}: {n_tg} targets\n")
        for t in terms:
            f.write(f"    - {t}\n")
print(f"[summary] report written -> {summary_path}")

print("\nValidation complete. All outputs in:", outdir)
