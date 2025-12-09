# ChIP-seq & Single-Cell Analysis Pipeline

Docker container for ChIP-seq and single-cell RNA-seq analysis, optimized for deployment on RunAI with VS Code Remote Tunnel support.

## Features

### ChIP-seq Analysis
- **trim_galore** - Adapter trimming
- **bowtie2** - Read alignment
- **samtools** - BAM manipulation
- **HOMER** - Peak calling (makeTagDirectory, findPeaks, makeUCSCfile)
- **igvtools** - TDF generation for IGV
- **fastqc** - Quality control
- **cutadapt** - Adapter removal
- **deeptools** - BAM/bigWig analysis
- **multiqc** - Aggregate QC reports
- **GNU parallel** - Job parallelization

### Single-Cell Analysis (Python)
- **scanpy** - Single-cell analysis framework
- **anndata** - Data structures for single-cell
- **leidenalg** - Community detection
- **bbknn** - Batch correction
- **harmonypy** - Batch correction

### Single-Cell Analysis (R)
- **Seurat** - Single-cell analysis
- **SingleCellExperiment** - Bioconductor framework
- **scater/scran** - QC and normalization
- **SingleR** - Cell type annotation

### 10x Genomics
- **cellranger** - Process 10x data (requires separate installation due to license)

## Quick Start

### Pull the image

```bash
docker pull ghcr.io/salk-harnessing-plants-initiative/chipseq-singlecell-pipeline:latest
```

### Run interactively

```bash
docker run -it \
  -v /path/to/data:/data \
  -v /path/to/outputs:/outputs \
  -v /path/to/references:/references \
  ghcr.io/salk-harnessing-plants-initiative/chipseq-singlecell-pipeline:latest
```

### Start VS Code Remote Tunnel

```bash
docker run -it \
  -v /path/to/data:/data \
  ghcr.io/salk-harnessing-plants-initiative/chipseq-singlecell-pipeline:latest \
  tunnel
```

Follow the authentication prompts to connect your VS Code.

## ChIP-seq Pipeline

The container includes Guanghui's ChIP-seq analysis scripts:

```bash
# Set environment variables
export REFERENCE=/references/bowtie2_index
export CPU_CORES=12
export READ_TYPE=PE
export IP_SAMPLE=IP_col0
export INPUT_SAMPLE=input_col0

# Run the pipeline
cd /data
/scripts/chipseq/ChIProfiler_M00_Basic.sh
```

### Environment Variables

| Variable | Default | Description |
|----------|---------|-------------|
| `REFERENCE` | `/references/bowtie2_index` | Path to bowtie2 index |
| `CPU_CORES` | `12` | Number of CPU cores for parallel |
| `READ_TYPE` | `PE` | `PE` (paired-end) or `SE` (single-end) |
| `IP_SAMPLE` | `IP_col0` | IP sample name prefix |
| `INPUT_SAMPLE` | `input_col0` | Input/control sample name prefix |

## Single-Cell Analysis

### Python (Scanpy)

Use `uv run` to execute Python commands:

```bash
# Run Python interactively
uv run python

# Run a script
uv run python /data/analysis.py

# Run Python tools
uv run multiqc .
uv run cutadapt --help
```

Example scanpy workflow:

```python
import scanpy as sc
import anndata as ad

# Load data
adata = sc.read_h5ad('/data/my_data.h5ad')

# Standard workflow
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata)
sc.tl.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata)
```

### R (Seurat)

```r
library(Seurat)

# Load data
seurat_obj <- Read10X('/data/filtered_feature_bc_matrix/')
seurat_obj <- CreateSeuratObject(seurat_obj)

# Standard workflow
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)
seurat_obj <- FindNeighbors(seurat_obj)
seurat_obj <- FindClusters(seurat_obj)
seurat_obj <- RunUMAP(seurat_obj)
```

## RunAI Deployment

Uses **RunAI CLI v2** syntax. See [docs/RUNAI_QUICK_REFERENCE.md](docs/RUNAI_QUICK_REFERENCE.md) for full reference.

### Submit VS Code Tunnel Workspace

```bash
runai workspace submit vscode-tunnel \
  --project talmo-lab \
  --image ghcr.io/salk-harnessing-plants-initiative/chipseq-singlecell-pipeline:latest \
  --cpu-core-request 8 \
  --cpu-memory-request 32G \
  --host-path path=/hpi/path/to/data,mount=/data,mount-propagation=HostToContainer \
  --command -- tunnel
```

### Submit ChIP-seq Analysis

```bash
runai workspace submit chipseq-analysis \
  --project talmo-lab \
  --image ghcr.io/salk-harnessing-plants-initiative/chipseq-singlecell-pipeline:latest \
  --cpu-core-request 16 \
  --cpu-memory-request 64G \
  --host-path path=/hpi/path/to/data,mount=/data,mount-propagation=HostToContainer \
  --host-path path=/hpi/path/to/outputs,mount=/outputs,mount-propagation=HostToContainer,readwrite \
  --environment REFERENCE=/references/bowtie2_index \
  --environment CPU_CORES=16 \
  --command -- /scripts/chipseq/ChIProfiler_M00_Basic.sh
```

### Monitor and Manage Jobs

```bash
runai workspace list                              # List all workspaces
runai workspace logs <NAME> -p talmo-lab --follow # View logs
runai workspace delete <NAME> -p talmo-lab        # Delete workspace
```

## Mount Points

| Path | Purpose |
|------|---------|
| `/workspace` | Working directory |
| `/data` | Input data (FASTQ, h5ad, etc.) |
| `/outputs` | Analysis outputs |
| `/references` | Reference genomes, indices |
| `/logs` | Log files |

## cellranger Installation

Due to 10x Genomics licensing, cellranger must be installed separately:

1. Download from [10x Genomics](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest)
2. Mount to `/opt/cellranger` or add to PATH

## Development

### Build locally

```bash
docker build -t chipseq-singlecell:dev .
```

### Run tests

```bash
docker run --rm chipseq-singlecell:dev help
```

## License

MIT License - See [LICENSE](LICENSE) for details.

## Acknowledgments

- ChIP-seq scripts by Guanghui Xu (Law Lab)
- Container development by Salk Harnessing Plants Initiative
