# Test Docker Container

Test the ChIP-seq & Single-Cell Pipeline Docker container.

## Quick Test

```bash
# Test help command
docker run --rm chipseq-singlecell-pipeline:local help
```

## Test All Tools

```bash
# ChIP-seq tools
docker run --rm chipseq-singlecell-pipeline:local bowtie2 --version
docker run --rm chipseq-singlecell-pipeline:local samtools --version
docker run --rm chipseq-singlecell-pipeline:local trim_galore --version
docker run --rm chipseq-singlecell-pipeline:local fastqc --version
docker run --rm chipseq-singlecell-pipeline:local findPeaks 2>&1 | head -1
docker run --rm chipseq-singlecell-pipeline:local igvtools version

# Python packages
docker run --rm chipseq-singlecell-pipeline:local python -c "import scanpy; print(f'scanpy: {scanpy.__version__}')"
docker run --rm chipseq-singlecell-pipeline:local python -c "import anndata; print(f'anndata: {anndata.__version__}')"
docker run --rm chipseq-singlecell-pipeline:local python -c "import scvi; print(f'scvi-tools: {scvi.__version__}')"

# R packages
docker run --rm chipseq-singlecell-pipeline:local R -e "library(Seurat); packageVersion('Seurat')"
docker run --rm chipseq-singlecell-pipeline:local R -e "library(SingleCellExperiment); packageVersion('SingleCellExperiment')"

# VS Code CLI
docker run --rm chipseq-singlecell-pipeline:local code --version
```

## Interactive Test

```bash
# Start interactive shell
docker run -it --rm chipseq-singlecell-pipeline:local bash

# Inside container, test tools manually
bowtie2 --version
python -c "import scanpy as sc; print(sc.__version__)"
R -e "library(Seurat)"
```

## Test with Mounted Data

```bash
# Mount test data directory
docker run -it --rm \
  -v /path/to/test/data:/data:ro \
  -v /path/to/test/outputs:/outputs \
  chipseq-singlecell-pipeline:local bash

# Inside container, verify mounts
ls /data
ls /outputs
```

## Test Entrypoint Commands

```bash
# Test help
docker run --rm chipseq-singlecell-pipeline:local help

# Test bash (default)
docker run -it --rm chipseq-singlecell-pipeline:local

# Test tunnel (will prompt for auth)
docker run -it --rm chipseq-singlecell-pipeline:local tunnel

# Test jupyter
docker run --rm -p 8888:8888 chipseq-singlecell-pipeline:local jupyter
```

## Test ChIP-seq Pipeline

```bash
# Verify pipeline script exists and is executable
docker run --rm chipseq-singlecell-pipeline:local ls -la /scripts/chipseq/

# Check script syntax
docker run --rm chipseq-singlecell-pipeline:local bash -n /scripts/chipseq/ChIProfiler_M00_Basic.sh
```

## Test Python Environment

```bash
docker run --rm chipseq-singlecell-pipeline:local python << 'EOF'
import scanpy as sc
import anndata as ad
import numpy as np

# Create test AnnData object
adata = ad.AnnData(np.random.randn(100, 50))
print(f"Created AnnData: {adata.shape}")

# Test basic scanpy operations
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.pca(adata)
print("Scanpy operations successful")
EOF
```

## Test R Environment

```bash
docker run --rm chipseq-singlecell-pipeline:local R --vanilla << 'EOF'
library(Seurat)
library(Matrix)

# Create test Seurat object
counts <- Matrix(rpois(1000, 5), nrow = 100, ncol = 10, sparse = TRUE)
rownames(counts) <- paste0("gene", 1:100)
colnames(counts) <- paste0("cell", 1:10)
seurat_obj <- CreateSeuratObject(counts = counts)
print(seurat_obj)
print("Seurat operations successful")
EOF
```

## Comprehensive Test Script

```bash
#!/bin/bash
# test-container.sh

IMAGE="${1:-chipseq-singlecell-pipeline:local}"
echo "Testing image: $IMAGE"

echo ""
echo "=== Testing ChIP-seq Tools ==="
docker run --rm $IMAGE bowtie2 --version | head -1
docker run --rm $IMAGE samtools --version | head -1
docker run --rm $IMAGE trim_galore --version 2>&1 | head -1
docker run --rm $IMAGE fastqc --version

echo ""
echo "=== Testing Python Packages ==="
docker run --rm $IMAGE python -c "
import scanpy, anndata, scvi, pandas, numpy
print(f'scanpy: {scanpy.__version__}')
print(f'anndata: {anndata.__version__}')
print(f'scvi-tools: {scvi.__version__}')
print(f'pandas: {pandas.__version__}')
print(f'numpy: {numpy.__version__}')
"

echo ""
echo "=== Testing R Packages ==="
docker run --rm $IMAGE R --vanilla -e "
packages <- c('Seurat', 'SingleCellExperiment', 'ggplot2')
for (pkg in packages) {
  cat(sprintf('%s: %s\n', pkg, as.character(packageVersion(pkg))))
}
"

echo ""
echo "=== Testing Entrypoint ==="
docker run --rm $IMAGE help | head -10

echo ""
echo "=== All tests passed! ==="
```

## Expected Output

A successful test should show:
- All ChIP-seq tools reporting versions
- Python packages importing without errors
- R packages loading without errors
- Entrypoint responding to `help` command

## Troubleshooting

### "Command not found"
- Tool may not be in PATH
- Check Dockerfile for installation step
- Verify build completed without errors

### Python import errors
- Package may not be installed
- Check `uv pip list` in container
- Verify python-env build stage

### R package errors
- Package may not be installed
- Check `installed.packages()` in R
- Verify r-env build stage

### Permission denied
- Check file permissions in container
- Scripts should be executable (`chmod +x`)

## Related Commands

- `/docker-build` - Build the image
- `/validate-bash` - Validate scripts
