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

## RunAI Integration

The `.claude/skills/runai.md` skill provides RunAI CLI v2 command reference for job submission.

### Available Slash Commands

| Command | Purpose |
|---------|---------|
| `/submit-chipseq` | Submit ChIP-seq analysis job |
| `/submit-singlecell` | Submit single-cell analysis job |
| `/submit-tunnel` | Submit VS Code tunnel workspace |
| `/monitor-jobs` | Monitor RunAI job status |
| `/cleanup-jobs` | Clean up completed/failed jobs |
| `/docker-build` | Build Docker image locally |
| `/docker-test` | Test Docker container |
| `/validate-bash` | Validate shell scripts with ShellCheck |

### OpenSpec Commands

| Command | Purpose |
|---------|---------|
| `/openspec:proposal` | Scaffold new change proposal |
| `/openspec:apply` | Implement approved change |
| `/openspec:archive` | Archive deployed change |

See `docs/RUNAI_QUICK_REFERENCE.md` for command examples.
