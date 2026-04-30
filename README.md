# LinTopoGRN

This repository contains the code, configuration, and organized result outputs for the manuscript:

**Lineage-aware inference of developmental gene regulatory networks from lineage-resolved single-cell data with LinTopoGRN**

`LinTopoGRN` is a lineage-aware framework for inferring developmental gene regulatory networks from lineage-resolved single-cell expression data. The code in this repository covers lineage-kernel construction, model training, downstream network analysis, and prediction validation.

## Environment Setup

Create the Conda environment with the pinned dependencies used for the main workflow:

```bash
conda env create -f environment.yml
conda activate LinTopoGRN
```

## Repository Layout

- `src/`
  Main publication-facing source code for `LinTopoGRN`.
- `data/`
  Input datasets used for training and evaluation.
- `data/reference_networks/`
  Reference regulatory networks used as priors or external validation resources.
- `configs/`
  Historical configuration files retained for reference.
- `results/main/LinTopoGRN/`
  Main model outputs used for manuscript figures, tables, and interpretation.
- `results/main/LinTopoGRN/validation_report/`
  Validation and enrichment outputs derived from top-ranked predictions.
- `results/supplementary/c_elegans/`
  Supplementary outputs for the C. elegans analysis.
- `results/baseline/`
  Baseline outputs from earlier comparison runs.
- `results/archive/LinTopoGRN_runs/`
  Archived historical runs and cached intermediate outputs.
- `docs/`
  Supporting documentation, including a compact results index.


## Required Inputs

The main training script expects the following inputs:

1. An expression matrix in CSV format.
   Rows should correspond to genes and columns to cells or samples.
2. A reference regulatory network in CSV format.
   The script accepts columns such as `gene_i/gene_j`, `tf/target`, `src/dst`, or `gene1/gene2`.
3. A lineage kernel in `.npy` format.
   This should be a cell-by-cell matrix aligned to the columns of the expression matrix.

The current default inputs in the main script are configured for the mouse embryo dataset included in `data/`.

## Quick Start

Build the lineage kernel:

```bash
python src/build_LinTopoGRN_lineage_kernel.py \
  --newick data/Chan_mouse_embryo/embryo1/embryo1_all.newick \
  --expr data/Chan_mouse_embryo/embryo1/preembryo1_sample.csv \
  --out data/Chan_mouse_embryo/embryo1/lineage_K.npy \
```

Run `LinTopoGRN`:

```bash
python src/run_LinTopoGRN.py
```



For a concise map of output files, see [docs/RESULTS_INDEX.md](/home/mengrui/gnn_grn%20copy/docs/RESULTS_INDEX.md:1).
