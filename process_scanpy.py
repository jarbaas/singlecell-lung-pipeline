import logging
import argparse
import random
import numpy as np
import anndata as ad
import scanpy as sc
import pandas as pd
from pathlib import Path
import celltypist
from celltypist import models
import gc

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
parser.add_argument('--max_mt', type=float, default=15.0, help="Maximum allowable mitochondrial percentage, 15.0 by default")
parser.add_argument("-t", "--threads", type=int, default=1, help="Number of CPU threads")
parser.add_argument("-s", "--seed", type=int, default=42, help="Random seed for reproducibility")
args = parser.parse_args()

# Global Scanpy settings
sc.settings.verbosity = 0
sc.settings.n_jobs = args.threads
np.random.seed(args.seed)
random.seed(args.seed)
sc.set_figure_params(dpi=100, fontsize=10, dpi_save=400, facecolor='white', figsize=6, format='png')

logging.info(f"Initialized pipeline with {args.threads} threads and seed {args.seed}")

# Data Ingestion (matrix or .h5ad inputs)
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
    metadata = pd.read_csv(args.metadata, index_col=0)
    adata.obs = adata.obs.join(metadata)
    if args.batch_key not in adata.obs.columns:
        raise KeyError(f"'{args.batch_key}' not found after merging metadata. Check column name in {args.metadata}.")
    if adata.obs[args.batch_key].isna().all():
        raise ValueError(f"Metadata merge failed. All values for '{args.batch_key}' are NaN. Verify your matrix barcodes and metadata index match exactly.")
        
# Data Validation
if adata.n_obs == 0 or adata.n_vars == 0:
    raise ValueError("Matrix is empty. Terminating pipeline.")

logging.info(f"Validation passed. Matrix: {adata.n_obs} cells x {adata.n_vars} genes.")

# Perform QC
adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, log1p=True)

# Drop cells missing from metadata file
adata = adata[adata.obs[args.batch_key].notna()].copy()
logging.info(f"Filtered matrix to cells with valid metadata. New shape: {adata.n_obs} cells x {adata.n_vars} genes.")

# Save QC plots
logging.info("Generating QC plots...")
sc.pl.violin(adata, ["n_genes_by_counts", "total_counts", "pct_counts_mt"], log=True, multi_panel=True, jitter=False, save="_pre_qc.png", rasterized=True)
logging.info("QC completed and plots saved.")

# Basic Filtering
sc.pp.filter_cells(adata, min_genes=100)
sc.pp.filter_genes(adata, min_cells=3)
adata = adata[adata.obs['pct_counts_mt'] < args.max_mt].copy()
gc.collect()

# Doublet detection 
logging.info("Running doublet detection...")
sc.pp.scrublet(adata, batch_key=args.batch_key, random_state=args.seed)

# Filter doublets
adata = adata[~adata.obs['predicted_doublet']].copy()
gc.collect()
logging.info(f"Filtering complete. {adata.n_obs} cells retained.")

# Saving count data
adata.layers["counts"] = adata.X.copy()

# Normalize
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
logging.info("Normalization completed")

# Select features
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key=args.batch_key)
sc.pl.highly_variable_genes(adata, save="_hvg.png")
logging.info("Selected top 2000 variable genes and plot saved.")

# PCA
logging.info("Running PCA...")
sc.tl.pca(adata, svd_solver='arpack', random_state=args.seed)
    
# Batch Correction and Neighborhood Graph
if args.batch_key in adata.obs.columns and adata.obs[args.batch_key].nunique() > 1:
    logging.info(f"Running native Harmony integration on {args.batch_key}...")
    sc.external.pp.harmony_integrate(adata, args.batch_key, random_state=args.seed)
        
    logging.info("Building neighborhood graph (Harmony-adjusted)...")
    sc.pp.neighbors(adata, use_rep='X_pca_harmony', random_state=args.seed)
else:
    logging.info(f"Single batch detected for {args.batch_key}. Skipping integration.")
    
    logging.info("Building neighborhood graph (PCA)...")
    sc.pp.neighbors(adata, random_state=args.seed)

gc.collect()
logging.info("Dimensionality reduction completed.")

# Calculate and plot UMAP
logging.info("Calculating UMAP...")
sc.tl.umap(adata, random_state=args.seed)
sc.pl.umap(adata, color=args.batch_key, size=2, legend_loc='none', save="_sample_umap.png")
logging.info("UMAP completed and plot saved.")

# Automatic Clustering
logging.info("Pre-computing Leiden over-clustering for majority voting...")
sc.tl.leiden(adata, resolution=5.0, flavor="igraph", n_iterations=2, directed=False, key_added='leiden_over_clustering', random_state=args.seed) 

logging.info("Performing automatic clustering via CellTypist...")
model = models.Model.load(model=args.model)
predictions = celltypist.annotate(
    adata, 
    model=model, 
    majority_voting=True, 
    over_clustering='leiden_over_clustering'
)
adata = predictions.to_adata()

logging.info("Saving UMAP visualizations...")
sc.pl.umap(adata, color='majority_voting', legend_loc='on data', legend_fontsize=8, legend_fontoutline=2, save="_labeled_umap.png")

# Extract markers for each annotated cell type and save to CSV
logging.info("Clustering success. Computing marker genes via Welch's t-test...")
sc.tl.rank_genes_groups(adata, groupby='majority_voting', method='wilcoxon', use_raw=False)
sc.get.rank_genes_groups_df(adata, group=None).to_csv("marker_genes.csv", index=False)

# Save the final object
logging.info("Saving analyzed data to h5ad...")
adata.write_h5ad("analyzed_object.h5ad")
logging.info("Pipeline completed. Analyzed h5ad object saved.")