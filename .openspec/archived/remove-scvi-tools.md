# OpenSpec Proposal: Remove scvi-tools Dependency

**Change ID**: `remove-scvi-tools`
**Created**: 2025-12-09
**Status**: Implemented
**Approved**: 2025-12-09
**Implemented**: 2025-12-09

## Summary

Remove scvi-tools from the Python dependencies to eliminate the PyTorch requirement, significantly reducing Docker image build time and size.

## Motivation / Problem Statement

The current Docker build is failing on GitHub Actions due to disk space constraints:

1. **scvi-tools requires PyTorch** (~800MB CPU-only, ~5GB with CUDA)
2. **GitHub Actions runners have limited disk space** (~14GB free)
3. **The Docker image is already large** due to R packages, HOMER, and other bioinformatics tools
4. **PyTorch CUDA dependencies cause "No space left on device" errors** during build

### Build Failure Evidence

```
× Failed to download `nvidia-cublas-cu12==12.8.4.1`
╰─▶ failed to write to file: No space left on device (os error 28)
help: `nvidia-cublas-cu12` was included because `scvi-tools` depends on `torch`
```

### Usage Analysis

scvi-tools is used for advanced deep learning single-cell analysis:
- Probabilistic modeling (scVI, scANVI)
- VAE-based batch correction
- Differential expression with uncertainty

**However**, for most single-cell workflows in this pipeline:
- **scanpy + anndata** provide core functionality
- **bbknn + harmonypy** provide batch correction
- **leidenalg** provides clustering

scvi-tools can be installed separately in a dedicated GPU container if needed.

## Proposed Solution

### 1. Remove scvi-tools from pyproject.toml

Remove the `scvi-tools` dependency and all PyTorch-related configuration.

### 2. Fix Container Tool Access

Ensure all tools are properly accessible:
- Fix entrypoint.sh venv path (currently `/opt/venv`, should be `/opt/python-env/.venv`)
- Fix VS Code CLI download URL (returns 404)

### 3. Update Documentation

- Add note to README.md about scvi-tools availability
- Update entrypoint help text to remove scvi-tools reference
- Document how to install scvi-tools separately if needed

## Tool Access Summary

| Tool Type | Location | Access Method |
|-----------|----------|---------------|
| ChIP-seq CLI | `/opt/bowtie2`, `/opt/samtools/bin`, `/opt/homer/bin` | `$PATH` |
| Python packages | `/opt/python-env/.venv` | `uv run` or activated venv |
| Python CLI tools | `/opt/python-env/.venv/bin` | `uv run <tool>` (cutadapt, multiqc, etc.) |
| R packages | System R library | `R` or `Rscript` |
| TrimGalore | `/usr/local/bin/trim_galore` | Direct CLI |
| VS Code | `/usr/local/bin/code` | Direct CLI |

### uv Best Practices

1. **Use `uv run` for Python commands** - Automatically uses the project's virtual environment
2. **Use `uv sync` for installation** - Reads pyproject.toml, creates lockfile
3. **Keep pyproject.toml as source of truth** - All dependencies defined there
4. **No manual venv activation needed** - `uv run python` handles it

Example usage in container:
```bash
# Run Python script
uv run python script.py

# Run Python tool
uv run cutadapt --help
uv run multiqc .

# Interactive Python
uv run python

# Add new package at runtime (if needed)
uv add package-name
```

## Files Affected

### Modified Files

| File | Changes |
|------|---------|
| `pyproject.toml` | Remove scvi-tools, remove PyTorch index config |
| `Dockerfile` | Fix VS Code CLI download URL |
| `scripts/entrypoint.sh` | Fix venv path, update help text |
| `README.md` | Remove scvi-tools from features list |

## Testing Plan

1. **Docker Build Validation**
   - [ ] Verify GitHub Actions build completes without disk space errors
   - [ ] Verify image is pushed to GHCR
   - [ ] Verify image can be pulled on RunAI

2. **Python Environment Validation**
   - [ ] Verify scanpy imports correctly
   - [ ] Verify anndata imports correctly
   - [ ] Verify batch correction tools work (bbknn, harmonypy)

3. **Documentation Validation**
   - [ ] README accurately reflects available tools
   - [ ] Instructions for optional scvi-tools installation are clear

## Rollback Plan

If scvi-tools is required:
1. Re-add `"scvi-tools"` to pyproject.toml dependencies
2. Use a larger runner or split the build into multiple stages
3. Consider a separate GPU-enabled container for deep learning workflows

## Dependencies

None - this is a removal that simplifies the build.

## Alternatives Considered

1. **CPU-only PyTorch** - Attempted but still causes disk space issues
2. **Split build stages** - Complex and may still hit limits
3. **Larger runners** - Costs money, not ideal for open source
4. **Separate GPU container** - Good future option, out of scope for this change

## Approval Checklist

- [ ] Technical approach reviewed
- [ ] Impact on existing workflows assessed
- [ ] Documentation updates planned
- [ ] Testing plan sufficient
