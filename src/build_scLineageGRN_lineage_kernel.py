"""Build the lineage kernel used by scLineageGRN from a Newick tree."""

import argparse, os, numpy as np, pandas as pd
from scLineageGRN_model import build_lineage_K_from_newick

def main():
    parser = argparse.ArgumentParser(description="Build the lineage kernel used by scLineageGRN from a Newick tree.")
    parser.add_argument('--newick', required=True, help='谱系树文件 (.newick)')
    parser.add_argument('--expr', required=True, help='表达矩阵 CSV (genes x cells)')
    parser.add_argument('--out', default='lineage_K.npy', help='输出文件路径 (.npy 或 .npz)')
    parser.add_argument('--tau', type=float, default=4.0, help='谱系距离衰减系数 τ')
    parser.add_argument('--k_top', type=int, default=32, help='每个细胞保留的最近叶子数')
    args = parser.parse_args()

    print(f"[load] reading expression: {args.expr}")
    expr = pd.read_csv(args.expr, index_col=0)
    cell_names = expr.columns.astype(str).tolist()
    print(f"[expr] cells={len(cell_names)}, genes={expr.shape[0]}")

    print(f"[build] constructing lineage kernel from {args.newick} (tau={args.tau}, k_top={args.k_top})")
    K = build_lineage_K_from_newick(
        newick_path=args.newick,
        cell_names=cell_names,
        tau=args.tau,
        k_top=args.k_top
    )

    print(f"[save] saving lineage_K, shape={K.shape} → {args.out}")
    if args.out.endswith('.npz'):
        np.savez(args.out, K=K.astype(np.float32))
    else:
        np.save(args.out, K.astype(np.float32))

    print("[done] OK.")

if __name__ == "__main__":
    main()
