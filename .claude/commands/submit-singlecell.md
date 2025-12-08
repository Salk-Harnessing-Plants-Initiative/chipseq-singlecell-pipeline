# Submit Single-Cell Analysis Job

Submit a single-cell RNA-seq analysis job to the RunAI cluster.

## Python/Scanpy Analysis

```bash
runai workspace submit singlecell-scanpy-<NAME> \
  --project talmo-lab \
  --image ghcr.io/salk-harnessing-plants-initiative/chipseq-singlecell-pipeline:latest \
  --cpu-core-request 32 \
  --cpu-memory-request 128G \
  --gpu-devices-request 1 \
  --host-path path=<DATA_PATH>,mount=/data,mount-propagation=HostToContainer \
  --host-path path=<OUTPUT_PATH>,mount=/outputs,mount-propagation=HostToContainer,readwrite \
  --host-path path=<WORKSPACE_PATH>,mount=/workspace,mount-propagation=HostToContainer,readwrite \
  --command -- python /workspace/<SCRIPT_NAME>.py
```

## R/Seurat Analysis

```bash
runai workspace submit singlecell-seurat-<NAME> \
  --project talmo-lab \
  --image ghcr.io/salk-harnessing-plants-initiative/chipseq-singlecell-pipeline:latest \
  --cpu-core-request 32 \
  --cpu-memory-request 128G \
  --host-path path=<DATA_PATH>,mount=/data,mount-propagation=HostToContainer \
  --host-path path=<OUTPUT_PATH>,mount=/outputs,mount-propagation=HostToContainer,readwrite \
  --host-path path=<WORKSPACE_PATH>,mount=/workspace,mount-propagation=HostToContainer,readwrite \
  --command -- Rscript /workspace/<SCRIPT_NAME>.R
```

## Parameters to Replace

| Placeholder | Description | Example |
|-------------|-------------|---------|
| `<NAME>` | Job name suffix | `arabidopsis-root` |
| `<DATA_PATH>` | Path to h5ad/10x files | `/hpi/hpi_dev/users/eberrigan/singlecell/data` |
| `<OUTPUT_PATH>` | Path for outputs | `/hpi/hpi_dev/users/eberrigan/singlecell/outputs` |
| `<WORKSPACE_PATH>` | Path to analysis scripts | `/hpi/hpi_dev/users/eberrigan/singlecell/scripts` |
| `<SCRIPT_NAME>` | Your analysis script | `analyze_root_cells` |

## Resource Recommendations

| Dataset Size | CPU | Memory | GPU |
|-------------|-----|--------|-----|
| < 100K cells | 16 | 64G | Optional |
| 100K - 500K cells | 32 | 128G | Recommended |
| 500K - 1M cells | 32 | 192G | Recommended |
| > 1M cells | 32 | 256G | Required |

## Example: Scanpy Analysis (500K cells)

```bash
runai workspace submit singlecell-scanpy-root-atlas \
  --project talmo-lab \
  --image ghcr.io/salk-harnessing-plants-initiative/chipseq-singlecell-pipeline:latest \
  --cpu-core-request 32 \
  --cpu-memory-request 192G \
  --gpu-devices-request 1 \
  --host-path path=/hpi/hpi_dev/users/eberrigan/singlecell/data,mount=/data,mount-propagation=HostToContainer \
  --host-path path=/hpi/hpi_dev/users/eberrigan/singlecell/outputs,mount=/outputs,mount-propagation=HostToContainer,readwrite \
  --host-path path=/hpi/hpi_dev/users/eberrigan/singlecell/scripts,mount=/workspace,mount-propagation=HostToContainer,readwrite \
  --command -- python /workspace/analyze_root_atlas.py
```

## Example: Seurat Analysis

```bash
runai workspace submit singlecell-seurat-leaf \
  --project talmo-lab \
  --image ghcr.io/salk-harnessing-plants-initiative/chipseq-singlecell-pipeline:latest \
  --cpu-core-request 32 \
  --cpu-memory-request 128G \
  --host-path path=/hpi/hpi_dev/users/eberrigan/singlecell/data,mount=/data,mount-propagation=HostToContainer \
  --host-path path=/hpi/hpi_dev/users/eberrigan/singlecell/outputs,mount=/outputs,mount-propagation=HostToContainer,readwrite \
  --host-path path=/hpi/hpi_dev/users/eberrigan/singlecell/scripts,mount=/workspace,mount-propagation=HostToContainer,readwrite \
  --command -- Rscript /workspace/analyze_leaf.R
```

## Available Python Packages

- **Core**: scanpy, anndata, numpy, scipy, pandas
- **Deep learning**: scvi-tools (GPU-accelerated)
- **Integration**: bbknn, harmonypy
- **Clustering**: leidenalg
- **Visualization**: matplotlib, seaborn, plotly

## Available R Packages

- **Core**: Seurat, SeuratObject, Matrix
- **Bioconductor**: SingleCellExperiment, scater, scran, DropletUtils
- **Annotation**: SingleR, celldex
- **Visualization**: ggplot2, patchwork, cowplot

## Monitor Job

```bash
# Check status
runai workspace list | grep singlecell

# View logs
runai workspace logs singlecell-scanpy-root-atlas -p talmo-lab --follow

# Detailed status
runai workspace describe singlecell-scanpy-root-atlas -p talmo-lab
```

## Sample Analysis Script (Python)

```python
#!/usr/bin/env python
# /workspace/analyze_root_atlas.py

import scanpy as sc
import anndata as ad

# Load data
adata = sc.read_h5ad('/data/root_atlas.h5ad')

# Standard workflow
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata)

# Save results
adata.write('/outputs/root_atlas_analyzed.h5ad')
sc.pl.umap(adata, color='leiden', save='/outputs/umap_clusters.png')
```

## Sample Analysis Script (R)

```r
#!/usr/bin/env Rscript
# /workspace/analyze_leaf.R

library(Seurat)

# Load data
data <- Read10X('/data/leaf_10x/')
seurat_obj <- CreateSeuratObject(counts = data, min.cells = 3, min.features = 200)

# Standard workflow
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)
seurat_obj <- FindNeighbors(seurat_obj)
seurat_obj <- FindClusters(seurat_obj)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)

# Save results
saveRDS(seurat_obj, '/outputs/leaf_analyzed.rds')
pdf('/outputs/umap_clusters.pdf')
DimPlot(seurat_obj, reduction = 'umap')
dev.off()
```

## Troubleshooting

### Out of memory (OOMKilled)
- Increase `--cpu-memory-request`
- For large datasets, use sparse matrices
- Consider subsetting data for initial exploration

### GPU not available
- Verify GPU request: `--gpu-devices-request 1`
- Check cluster GPU availability
- scvi-tools will fall back to CPU if no GPU

### Script not found
- Verify `WORKSPACE_PATH` contains your script
- Check script permissions: `chmod +x script.py`

## Related Commands

- `/submit-tunnel` - Interactive VS Code session
- `/monitor-jobs` - Monitor job status
- `/cleanup-jobs` - Clean up completed jobs
