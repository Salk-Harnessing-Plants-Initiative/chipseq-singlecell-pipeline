# Manual Execution Using RunAI

This guide shows how to run ChIP-seq and single-cell analysis jobs using RunAI commands on the Salk HPI cluster.

## Prerequisites

- RunAI CLI authenticated: `runai login`
- Access to `talmo-lab` project (Kubernetes namespace: `runai-talmo-lab`)
- Data files accessible on cluster storage

## Quick Start - VS Code Tunnel

The fastest way to get started is with an interactive VS Code tunnel:

```bash
runai workspace submit vscode-myname \
  --project talmo-lab \
  --image ghcr.io/salk-harnessing-plants-initiative/chipseq-singlecell-pipeline:latest \
  --cpu-core-request 8 \
  --cpu-memory-request 32G \
  --host-path path=/hpi/path/to/your/data,mount=/data,mount-propagation=HostToContainer \
  --host-path path=/hpi/path/to/your/outputs,mount=/outputs,mount-propagation=HostToContainer,readwrite \
  --command -- tunnel
```

**Connect to the tunnel:**
1. Run `runai workspace logs vscode-myname -p talmo-lab --follow`
2. Look for the GitHub device authentication URL
3. Complete authentication in browser
4. Open VS Code → Remote Explorer → Tunnels → Connect

---

## ChIP-seq Analysis

### Step 1: Prepare Data

Ensure your data is organized on the cluster:
```
/hpi/path/to/data/
├── IP_col0_R1.fq.gz        # IP sample, read 1
├── IP_col0_R2.fq.gz        # IP sample, read 2
├── input_col0_R1.fq.gz     # Input control, read 1
└── input_col0_R2.fq.gz     # Input control, read 2

/hpi/path/to/references/
└── TAIR10/                  # Bowtie2 index
    ├── TAIR10.1.bt2
    ├── TAIR10.2.bt2
    ├── TAIR10.3.bt2
    ├── TAIR10.4.bt2
    ├── TAIR10.rev.1.bt2
    └── TAIR10.rev.2.bt2
```

### Step 2: Submit ChIP-seq Job

```bash
runai workspace submit chipseq-col0-h3k4me3 \
  --project talmo-lab \
  --image ghcr.io/salk-harnessing-plants-initiative/chipseq-singlecell-pipeline:latest \
  --cpu-core-request 16 \
  --cpu-memory-request 64G \
  --host-path path=/hpi/hpi_dev/users/myname/chipseq/data,mount=/data,mount-propagation=HostToContainer \
  --host-path path=/hpi/hpi_dev/users/myname/chipseq/outputs,mount=/outputs,mount-propagation=HostToContainer,readwrite \
  --host-path path=/hpi/hpi_dev/references/arabidopsis,mount=/references,mount-propagation=HostToContainer \
  --environment REFERENCE=/references/TAIR10 \
  --environment CPU_CORES=16 \
  --environment READ_TYPE=PE \
  --environment IP_SAMPLE=IP_col0_H3K4me3 \
  --environment INPUT_SAMPLE=input_col0 \
  --command -- /scripts/chipseq/ChIProfiler_M00_Basic.sh
```

### Step 3: Monitor Progress

```bash
# Check status
runai workspace list | grep chipseq

# View logs in real-time
runai workspace logs chipseq-col0-h3k4me3 -p talmo-lab --follow

# Detailed status
runai workspace describe chipseq-col0-h3k4me3 -p talmo-lab
```

### Step 4: Check Results

After completion, your outputs will contain:
```
/outputs/
├── IP_col0_H3K4me3_sorted_rmdup_uniq.bam      # Deduplicated alignments
├── IP_col0_H3K4me3_sorted_rmdup_uniq.bam.bai  # BAM index
├── IP_col0_H3K4me3_sorted_rmdup_uniq_TagDir/  # HOMER tag directory
├── IP_col0_H3K4me3.bedGraph.gz                # Normalized coverage
├── IP_col0_H3K4me3.tdf                        # IGV visualization
├── input_col0_sorted_rmdup_uniq.bam           # Input control
├── input_col0_sorted_rmdup_uniq_TagDir/       # Input tag directory
├── 1_res_read_mapping_stat.txt                # Mapping statistics
└── *_testPeaks*.bed                           # Called peaks
```

### Step 5: Cleanup

```bash
runai workspace delete chipseq-col0-h3k4me3 -p talmo-lab
```

---

## Single-Cell Analysis

### Python/Scanpy Workflow

#### Step 1: Create Analysis Script

Save this as `/hpi/path/to/scripts/analyze_cells.py`:

```python
#!/usr/bin/env python
import scanpy as sc
import anndata as ad

# Load data
adata = sc.read_h5ad('/data/my_dataset.h5ad')
print(f"Loaded {adata.n_obs} cells, {adata.n_vars} genes")

# Standard preprocessing
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Feature selection and dimensionality reduction
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata)

# Save results
adata.write('/outputs/analyzed.h5ad')
sc.pl.umap(adata, color='leiden', save='/outputs/umap_clusters.png')
print("Analysis complete!")
```

#### Step 2: Submit Job

```bash
runai workspace submit singlecell-root-atlas \
  --project talmo-lab \
  --image ghcr.io/salk-harnessing-plants-initiative/chipseq-singlecell-pipeline:latest \
  --cpu-core-request 32 \
  --cpu-memory-request 128G \
  --gpu-devices-request 1 \
  --host-path path=/hpi/hpi_dev/users/myname/singlecell/data,mount=/data,mount-propagation=HostToContainer \
  --host-path path=/hpi/hpi_dev/users/myname/singlecell/outputs,mount=/outputs,mount-propagation=HostToContainer,readwrite \
  --host-path path=/hpi/hpi_dev/users/myname/singlecell/scripts,mount=/workspace,mount-propagation=HostToContainer \
  --command -- python /workspace/analyze_cells.py
```

### R/Seurat Workflow

#### Step 1: Create Analysis Script

Save this as `/hpi/path/to/scripts/analyze_cells.R`:

```r
#!/usr/bin/env Rscript
library(Seurat)
library(ggplot2)

# Load 10x data
data <- Read10X('/data/10x_output/filtered_feature_bc_matrix/')
seurat_obj <- CreateSeuratObject(counts = data, min.cells = 3, min.features = 200)
print(seurat_obj)

# Standard workflow
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)
seurat_obj <- FindClusters(seurat_obj)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)

# Save results
saveRDS(seurat_obj, '/outputs/analyzed.rds')
pdf('/outputs/umap_clusters.pdf', width = 10, height = 8)
DimPlot(seurat_obj, reduction = 'umap', label = TRUE)
dev.off()

print("Analysis complete!")
```

#### Step 2: Submit Job

```bash
runai workspace submit singlecell-seurat-leaf \
  --project talmo-lab \
  --image ghcr.io/salk-harnessing-plants-initiative/chipseq-singlecell-pipeline:latest \
  --cpu-core-request 32 \
  --cpu-memory-request 128G \
  --host-path path=/hpi/hpi_dev/users/myname/singlecell/data,mount=/data,mount-propagation=HostToContainer \
  --host-path path=/hpi/hpi_dev/users/myname/singlecell/outputs,mount=/outputs,mount-propagation=HostToContainer,readwrite \
  --host-path path=/hpi/hpi_dev/users/myname/singlecell/scripts,mount=/workspace,mount-propagation=HostToContainer \
  --command -- Rscript /workspace/analyze_cells.R
```

---

## Resource Recommendations

### ChIP-seq Jobs

| Dataset Size | CPU | Memory | Typical Runtime |
|-------------|-----|--------|-----------------|
| < 50M reads | 8 | 32G | 1-2 hours |
| 50-100M reads | 12 | 48G | 2-4 hours |
| > 100M reads | 16 | 64G | 4-8 hours |

### Single-Cell Jobs

| Cells | CPU | Memory | GPU | Typical Runtime |
|-------|-----|--------|-----|-----------------|
| < 50K | 16 | 64G | Optional | 30 min - 1 hour |
| 50-200K | 32 | 128G | Recommended | 1-2 hours |
| 200K-1M | 32 | 192G | Recommended | 2-4 hours |
| > 1M | 32 | 256G | Required | 4-8 hours |

---

## Monitoring Jobs

### Real-time Dashboard

```bash
./scripts/monitor-runai-jobs.sh --watch
```

### Manual Monitoring

```bash
# List all jobs
runai workspace list

# Filter by type
runai workspace list | grep chipseq
runai workspace list | grep singlecell

# Count by status
runai workspace list | awk '{print $2}' | sort | uniq -c

# Detailed status
runai workspace describe <JOB_NAME> -p talmo-lab

# Follow logs
runai workspace logs <JOB_NAME> -p talmo-lab --follow
```

---

## Cleanup

### Using Helper Script (Recommended)

```bash
# Preview cleanup
./scripts/cleanup-runai.sh --completed --dry-run

# Delete completed jobs
./scripts/cleanup-runai.sh --completed --yes

# Delete failed jobs
./scripts/cleanup-runai.sh --failed --yes

# Delete all finished jobs
./scripts/cleanup-runai.sh --all --yes
```

### Manual Cleanup

```bash
# Delete single job
runai workspace delete <JOB_NAME> -p talmo-lab

# Delete all completed jobs
runai workspace list | grep "Succeeded" | awk '{print $1}' | \
  xargs -I {} runai workspace delete {} -p talmo-lab
```

---

## Troubleshooting

### Job Fails Immediately

```bash
# Check logs
runai workspace logs <JOB_NAME> -p talmo-lab

# Common issues:
# - "Image pull error" → Check image name/tag
# - "Permission denied" → Check host path permissions
# - "Command not found" → Check command syntax
```

### Job Stuck in Pending

```bash
# Check detailed status
runai workspace describe <JOB_NAME> -p talmo-lab

# Common issues:
# - Insufficient cluster resources → Wait or reduce resource request
# - Invalid host path → Verify path exists on cluster
```

### Out of Memory (OOMKilled)

Increase memory allocation:
```bash
--cpu-memory-request 128G   # Instead of 64G
--cpu-memory-request 256G   # For very large datasets
```

### VS Code Tunnel Not Connecting

1. Verify workspace is running:
   ```bash
   runai workspace list | grep vscode
   ```

2. Check logs for authentication URL:
   ```bash
   runai workspace logs vscode-myname -p talmo-lab --follow
   ```

3. Complete GitHub device authentication in browser

4. Refresh VS Code Remote Explorer

### Script Errors

1. Test script locally first:
   ```bash
   docker run -it --rm \
     -v /local/data:/data \
     ghcr.io/.../chipseq-singlecell-pipeline:latest \
     bash
   # Then run your script interactively
   ```

2. Check file permissions on cluster

3. Verify input paths are correct

---

## Additional Resources

- [Quick Reference](RUNAI_QUICK_REFERENCE.md) - Command cheat sheet
- [RunAI Documentation](https://docs.run.ai/)
- [Salk HPI RunAI Guide](https://researchit.salk.edu/runai/)
- [Container Image](https://github.com/salk-harnessing-plants-initiative/chipseq-singlecell-pipeline/pkgs/container/chipseq-singlecell-pipeline)
