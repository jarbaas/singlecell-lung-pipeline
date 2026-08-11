# Automatic Nextflow scRNA-seq Pipeline for IPF Lung Tissue Analysis

An automated, cloud-native downstream single-cell RNA-seq pipeline designed to process and annotate large datasets, with lung tissue samples in mind. Built for deployment on AWS, this pipeline orchestrates a Python-based Scanpy workflow featuring QC, dimensionality reduction, UMAP and automated cell type annotation within an isolated Docker environment.

## Architecture

* **Workflow Manager:** Nextflow (DSL2)

* **Compute & Storage:** AWS EC2 and Amazon S3

* **Containerization:** Docker (custom file with python:3.14-slim base)

* **Analytical Stack:** Scanpy, Scrublet, Harmonypy, CellTypist  

## Pipeline Workflow

* **Data Ingestion:** Dynamically pulls raw .mtx matrices or pre-compiled .h5ad objects and metadata directly from S3.

* **Quality Control:** Filters low-quality cells, thresholds mitochondrial reads, and removes predicted doublets via Scrublet.

* **Dimensionality Reduction & Batch Correction:** Performs PCA, integrates patient batches using Harmony, and generates UMAP projections.

* **Optimized Clustering:** Utilizes the igraph flavor of the Leiden algorithm at an optimized resolution to reduce processing time.

* **Automated Annotation:** Maps cell identities using CellTypist and performs majority-voting refinement, accepting use of a custom pre-trained input model.

* **Differential Expression:** Extracts cluster-specific marker genes using Welch's t-test.

* **Deliverables:** Pushes the finalized .h5ad object, marker gene CSV, and UMAP visualizations back to S3.

## Quick Start

To execute this pipeline on an AWS EC2 instance:

1) **Container Registry Setup:** Build the Docker image from the provided Dockerfile and push it to your Amazon Elastic Container Registry (ECR).

2) **Execution:** Run the pipeline natively. Nextflow will automatically pull the repository and configuration from GitHub. Pass your specific ECR image URI and S3 bucket paths via the command line:

```
nextflow run jarbaas/singlecell-lung-pipeline \
  --container_image '<your-ecr-image-uri>' \
  --input 's3://<your-bucket-name>/raw_data/' \
  --metadata 's3://<your-bucket-name>/metadata.csv' \
  --model 's3://<your-bucket-name>/model_noca.pkl' \
  --outdir 's3://<your-bucket-name>/results/'
```

## Production Optimizations
* **Flexible Data Ingestion:** The pipeline dynamically accepts either raw 10x `.mtx` directories or pre-compiled `.h5ad` objects,  increasing flexibility in processing upstream data formats.
* **Dynamic Metadata Handling:** Incorporates a customizable `batch_key` argument to manage varying external metadata setups, ensuring robust Harmony batch correction regardless of the dataset's specific column naming conventions.
* **Pre-Clustering for Annotation:** Performs Leiden pre-clustering to establish community labels at a defined resolution. CellTypist leverages these pre-computed labels for its majority-voting mechanism, bypassing redundant calculations and significantly optimizing runtime efficiency.

## Limitations
* An important consideration with this pipeline is the automatic processing of the genes.tsv (or features.tsv) file. The script will search the second column of the file, however, testing revealed certain files uploaded to GEO contain only one column. To fix, simply duplicate the column.
