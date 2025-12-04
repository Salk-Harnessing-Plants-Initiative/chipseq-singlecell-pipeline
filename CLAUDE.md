# Claude Code Instructions

## Project Overview

This repository contains a Docker container for ChIP-seq and single-cell RNA-seq analysis, designed for deployment on RunAI with VS Code Remote Tunnel support.

## Key Components

- **Dockerfile**: Multi-stage build with Ubuntu 22.04 base
- **scripts/chipseq/**: Guanghui's ChIP-seq analysis pipeline
- **scripts/entrypoint.sh**: Container entrypoint with command routing
- **.github/workflows/**: CI/CD for Docker builds and script validation
- **.devcontainer/**: VS Code devcontainer configuration

## Technology Stack

- **Python**: Managed by `uv`, includes scanpy/anndata/scvi-tools
- **R**: 4.3+ with Seurat and Bioconductor packages
- **ChIP-seq**: trim_galore, bowtie2, samtools, HOMER, igvtools
- **VS Code**: Remote Tunnel via `code tunnel` command

## CI/CD

- Images pushed to `ghcr.io/salk-harnessing-plants-initiative/chipseq-singlecell-pipeline`
- Workflow validates bash scripts with ShellCheck
- Tests verify all tools are installed correctly

## Common Tasks

### Adding new Python packages
```dockerfile
RUN uv pip install package-name
```

### Adding new R packages
```dockerfile
RUN R -e "install.packages('package-name', repos='https://cloud.r-project.org')"
```

### Updating ChIP-seq scripts
Scripts in `scripts/chipseq/` should maintain compatibility with:
- Environment variables for configuration
- `set -euo pipefail` for safety
- Proper quoting for paths with spaces
