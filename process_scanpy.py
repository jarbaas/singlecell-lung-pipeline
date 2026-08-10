import logging
import argparse
import anndata as ad
import scanpy as sc
import pandas as pd
from pathlib import Path
import celltypist

# Enable logging for tracking in the cloud
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

# argparse Setup
parser = argparse.ArgumentParser(description="Downstream scRNA-seq processing pipeline")
parser.add_argument("-i", "--input", required=True, help="Path to matrix folder or .h5ad file")
parser.add_argument("-mod", "--model", required=True, help="Path to CellTypist model PKL")
parser.add_argument("-m", "--metadata", required=False, help="Optional path to metadata CSV")
parser.add_argument("-b", "--batch-key", default="Sample_Name", help="Metadata column defining the batch")
parser.add_argument("-t", "--threads", type=int, default=1, help="Number of CPU threads")
parser.add_argument("-s", "--seed", type=int, default=42, help="Random seed for reproducibility")
args = parser.parse_args()

# Global Scanpy settings
sc.settings.verbosity = 0
sc.settings.n_jobs = args.threads
sc.settings.seed = args.seed
sc.settings.set_figure_params(dpi=100, fontsize=10, dpi_save=400, facecolor='white', figsize=(6,6), format='png')

logging.info(f"Initialized pipeline with {args.threads} threads and seed {args.seed}")

# Data Ingestion (supports matrix or .h5ad inputs)
input_path = Path(args.input)
if input_path.is_file() and input_path.suffix == '.h5ad':
    logging.info(f"Reading .h5ad file: {input_path}")
    adata = sc.read_h5ad(input_path)
elif input_path.is_dir():
    logging.info(f"Reading 10x matrix directory: {input_path}")
    adata = sc.read_10x_mtx(input_path, var_names='gene_symbols', make_unique=True, cache=False)
else:
    raise ValueError(f"Invalid input path: {input_path}")

# Metadata Integration
if args.metadata:
    logging.info(f"Merging external metadata: {args.metadata}")
    metadata = pd.read_csv(args.metadata, index_col='Unnamed: 0')
    adata.obs = adata.obs.join(metadata)

# Data Validation
if adata.n_obs == 0 or adata.n_vars == 0:
    raise ValueError("Matrix is empty. Terminating pipeline.")

if args.batch_key not in adata.obs.columns:
    raise KeyError(f"'{args.batch_key}' missing from adata.obs. Required for batch processing.")

logging.info(f"Validation passed. Matrix: {adata.n_obs} cells x {adata.n_vars} genes.")

# Perform QC
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, log1p=True)

# Doublet detection
sc.pp.scrublet(adata, batch_key=args.batch_key)

# Save QC plots
sc.pl.violin(adata, ["n_genes_by_counts", "total_counts", "pct_counts_mt"], multi_panel=True, save="_pre_qc.png")

logging.info("QC completed and plots saved.")

#Filtering
sc.pp.filter_cells(adata, min_genes=100)
sc.pp.filter_genes(adata, min_cells=3)
adata = adata[adata.obs['pct_counts_mt'] < 15].copy()
adata = adata[~adata.obs['predicted_doublet']].copy()
logging.info(f"Filtering complete. {adata.n_obs} cells retained.")

# Saving count data
adata.layers["counts"] = adata.X.copy()

# Normalize
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

logging.info("Normalization completed")

# Select features
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key=args.batch_key)
sc.pl.highly_variable_genes(adata, save="_hvg.png")
logging.info(f"Selected top 2000 variable genes and plot saved.")

# Dimensionality Reduction
sc.tl.pca(adata)
sc.pp.neighbors(adata)
logging.info(f"Dimensionality reduction completed and plot saved.")

sc.tl.umap(adata)
sc.pl.umap(adata, color=args.batch_key, size=2, save="_sample_umap.png")
logging.info(f"UMAP completed and plot saved.")

# Automatic Clustering
logging.info(f"Performing automatic clustering...")
predictions = celltypist.annotate(adata, model=args.model, majority_voting=True)
adata = predictions.to_adata()
sc.pl.umap(adata, color='majority_voting', save="_labeled_umap.png")

# Extract markers for each annotated cell type and save to CSV
logging.info(f"Clustering success. Writing markers for cell types to CSV...")
sc.tl.rank_genes_groups(adata, groupby='majority_voting', method='wilcoxon')
sc.get.rank_genes_groups_df(adata, group=None).to_csv("marker_genes.csv", index=False)

adata.write_h5ad("analyzed_object.h5ad")
logging.info(f"Pipeline completed. Analyzed h5ad object saved.")






