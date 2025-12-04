# Quick Start Guide for Charlotte

Welcome! This guide will help you get started with the ChIP-seq and single-cell analysis container on RunAI.

## 1. Connect to RunAI

First, submit an interactive job on RunAI:

```bash
runai submit charlotte-analysis \
  --image ghcr.io/salk-harnessing-plants-initiative/chipseq-singlecell-pipeline:latest \
  --gpu 1 \
  --cpu 16 \
  --memory 64G \
  --pvc your-data-pvc:/data \
  --pvc reference-pvc:/references \
  --interactive \
  -- tunnel
```

## 2. Connect VS Code

Once the job starts, you'll see output like:

```
To grant access to the server, please log into https://github.com/login/device and use code XXXX-XXXX
```

1. Go to the URL shown
2. Enter the code
3. Open VS Code on your laptop
4. Press `Ctrl+Shift+P` (or `Cmd+Shift+P` on Mac)
5. Type "Remote-Tunnels: Connect to Tunnel"
6. Select your tunnel

You're now connected!

## 3. Run ChIP-seq Analysis

### Set up your analysis

```bash
# Go to your data directory
cd /data

# Check your FASTQ files are there
ls *.fastq.gz
```

### Configure the pipeline

Edit the environment variables in the script or export them:

```bash
# Set your reference genome path
export REFERENCE=/references/TAIR10/bowtie2_index

# Set number of CPUs (adjust based on your allocation)
export CPU_CORES=16

# Paired-end or single-end
export READ_TYPE=PE

# Your sample names
export IP_SAMPLE=IP_col0
export INPUT_SAMPLE=input_col0
```

### Run the pipeline

```bash
/scripts/chipseq/ChIProfiler_M00_Basic.sh
```

This will:
1. Trim adapters (trim_galore)
2. Align reads (bowtie2)
3. Remove duplicates (samtools)
4. Create tag directories (HOMER)
5. Generate browser tracks
6. Call peaks

### Check outputs

- `*_sorted_rmdup_uniq.bam` - Final aligned BAM files
- `*_TagDir/` - HOMER tag directories
- `*.tdf` - IGV tracks
- `*testPeaks*.bed` - Called peaks
- `1_res_read_mapping_stat.txt` - Mapping statistics

## 4. Run Single-Cell Analysis

### Python (Scanpy) - For large datasets

```bash
# Start Python
python

# Or use Jupyter
jupyter lab --ip=0.0.0.0 --port=8888
```

```python
import scanpy as sc

# Load your data
adata = sc.read_10x_h5('/data/filtered_feature_bc_matrix.h5')

# Basic QC
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# Normalize
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Find variable genes
sc.pp.highly_variable_genes(adata, n_top_genes=2000)

# PCA
sc.tl.pca(adata)

# Neighbors and clustering
sc.pp.neighbors(adata, n_pcs=30)
sc.tl.leiden(adata)

# UMAP
sc.tl.umap(adata)

# Save
adata.write('/outputs/my_analysis.h5ad')
```

### R (Seurat)

```bash
# Start R
R
```

```r
library(Seurat)

# Load data
data <- Read10X('/data/filtered_feature_bc_matrix/')
seurat <- CreateSeuratObject(data, min.cells = 3, min.features = 200)

# QC
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")
seurat <- subset(seurat, subset = nFeature_RNA > 200 & percent.mt < 20)

# Normalize and find variable features
seurat <- NormalizeData(seurat)
seurat <- FindVariableFeatures(seurat)

# Scale and PCA
seurat <- ScaleData(seurat)
seurat <- RunPCA(seurat)

# Cluster
seurat <- FindNeighbors(seurat, dims = 1:30)
seurat <- FindClusters(seurat, resolution = 0.5)
seurat <- RunUMAP(seurat, dims = 1:30)

# Save
saveRDS(seurat, '/outputs/seurat_analysis.rds')
```

## 5. Tips for Large Datasets (3M cells)

For very large datasets, use these memory-efficient approaches:

### Scanpy with backed mode

```python
# Read in backed mode (keeps data on disk)
adata = sc.read_h5ad('/data/large_data.h5ad', backed='r')
```

### Batch processing

```python
# Process in chunks
import anndata as ad

# Concatenate multiple samples efficiently
adatas = [sc.read_h5ad(f) for f in sample_files]
adata = ad.concat(adatas, join='outer')
```

### Request more resources

```bash
# For 3M cells, request at least 128GB RAM
runai submit charlotte-bigdata \
  --image ghcr.io/salk-harnessing-plants-initiative/chipseq-singlecell-pipeline:latest \
  --gpu 1 \
  --cpu 32 \
  --memory 128G \
  ...
```

## 6. File Transfer

Since you have server access mounted to the cluster:

1. Upload files to the mounted server path
2. They'll appear in `/data` inside the container
3. Save outputs to `/outputs`
4. Download from the mounted server path

## Need Help?

- Check the full documentation: [README.md](README.md)
- Container help: Run `help` inside the container
- ChIP-seq questions: Ask Guanghui
- Container issues: Ask Elizabeth
