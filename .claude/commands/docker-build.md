# Build Docker Image

Build the ChIP-seq & Single-Cell Pipeline Docker image locally.

## Quick Build

```bash
docker build -t chipseq-singlecell-pipeline:local .
```

## Build with Specific Tag

```bash
docker build -t chipseq-singlecell-pipeline:dev .
docker build -t chipseq-singlecell-pipeline:$(git rev-parse --short HEAD) .
```

## Build for Registry

```bash
# Build with full registry path
docker build -t ghcr.io/salk-harnessing-plants-initiative/chipseq-singlecell-pipeline:latest .

# Build with version tag
docker build -t ghcr.io/salk-harnessing-plants-initiative/chipseq-singlecell-pipeline:v1.0.0 .
```

## Build Options

| Option | Description | Example |
|--------|-------------|---------|
| `-t` | Tag the image | `-t myimage:v1` |
| `--no-cache` | Build without cache | `--no-cache` |
| `--platform` | Target platform | `--platform linux/amd64` |
| `--progress` | Build output style | `--progress=plain` |

## Full Build Command

```bash
docker build \
  --platform linux/amd64 \
  --progress=plain \
  -t ghcr.io/salk-harnessing-plants-initiative/chipseq-singlecell-pipeline:latest \
  .
```

## Build Stages

The Dockerfile uses multi-stage builds:

1. **base** - Ubuntu 22.04 with system dependencies
2. **chipseq-tools** - Bowtie2, Samtools, HOMER, IGVtools, FastQC, TrimGalore
3. **python-env** - Python 3.11 with uv, scanpy, scvi-tools
4. **r-env** - R 4.3+ with Seurat, Bioconductor
5. **cellranger-env** - Placeholder for 10x cellranger
6. **runtime** - Final image with VS Code CLI

## Build Specific Stage

```bash
# Build only up to python-env stage
docker build --target python-env -t chipseq-python:test .

# Build only up to chipseq-tools stage
docker build --target chipseq-tools -t chipseq-tools:test .
```

## Verify Build

After building, verify tools are installed:

```bash
# Run test container
docker run --rm chipseq-singlecell-pipeline:local help

# Check specific tools
docker run --rm chipseq-singlecell-pipeline:local bowtie2 --version
docker run --rm chipseq-singlecell-pipeline:local samtools --version
docker run --rm chipseq-singlecell-pipeline:local python -c "import scanpy; print(scanpy.__version__)"
docker run --rm chipseq-singlecell-pipeline:local R -e "library(Seurat); packageVersion('Seurat')"
```

## Push to Registry

```bash
# Login to GHCR
echo $GITHUB_TOKEN | docker login ghcr.io -u USERNAME --password-stdin

# Push image
docker push ghcr.io/salk-harnessing-plants-initiative/chipseq-singlecell-pipeline:latest
```

## CI/CD Build

The GitHub Actions workflow automatically builds on:
- Push to `main` branch
- Version tags (`v*.*.*`)
- Pull requests

Images are pushed to:
```
ghcr.io/salk-harnessing-plants-initiative/chipseq-singlecell-pipeline
```

## Troubleshooting

### Build fails at package installation
- Check network connectivity
- Try `--no-cache` to rebuild from scratch
- Check if package repositories are available

### Out of disk space
```bash
# Clean up Docker resources
docker system prune -a
docker builder prune
```

### Build is slow
- Use BuildKit: `DOCKER_BUILDKIT=1 docker build ...`
- Ensure layer caching is working
- Consider using `--mount=type=cache` for package managers

## Related Commands

- `/docker-test` - Test the built image
- `/validate-bash` - Validate scripts before building
