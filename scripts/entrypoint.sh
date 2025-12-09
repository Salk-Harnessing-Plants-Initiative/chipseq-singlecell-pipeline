#!/bin/bash
# ==============================================================================
# Entrypoint script for ChIP-seq & Single-Cell Pipeline container
# ==============================================================================
set -euo pipefail

# Change to the uv project directory for uv run to work
cd /opt/python-env

# Handle commands
case "${1:-bash}" in
    tunnel|code-tunnel)
        echo "Starting VS Code Remote Tunnel..."
        echo "Follow the authentication prompts to connect."
        exec code tunnel --accept-server-license-terms
        ;;
    jupyter|lab)
        echo "Starting JupyterLab..."
        exec uv run jupyter lab --ip=0.0.0.0 --port=8888 --no-browser --allow-root
        ;;
    python)
        # Run Python with uv
        shift
        exec uv run python "$@"
        ;;
    help|--help|-h)
        cat << 'EOF'
ChIP-seq & Single-Cell Analysis Pipeline Container

USAGE:
    docker run [options] <image> [command]

COMMANDS:
    bash            Start interactive bash shell (default)
    tunnel          Start VS Code Remote Tunnel
    jupyter         Start JupyterLab server
    python          Run Python with uv (e.g., python script.py)
    help            Show this help message

PYTHON (using uv):
    # Run Python script
    uv run python script.py

    # Run Python tool
    uv run cutadapt --help
    uv run multiqc .

    # Interactive Python
    uv run python

    # Add package at runtime
    uv add package-name

EXAMPLES:
    # Interactive shell
    docker run -it -v /data:/data <image>

    # VS Code Remote Tunnel
    docker run -it <image> tunnel

    # JupyterLab
    docker run -p 8888:8888 <image> jupyter

    # Run Python script
    docker run -v /data:/data <image> python /data/analysis.py

AVAILABLE TOOLS:
    ChIP-seq:       trim_galore, bowtie2, samtools, findPeaks, makeTagDirectory,
                    makeUCSCfile, igvtools, fastqc, cutadapt, multiqc, deeptools
    Single-cell:    scanpy, anndata, bbknn, harmonypy, Seurat
    General:        parallel, R, python, uv

CHIP-SEQ SCRIPTS:
    /scripts/chipseq/ChIProfiler_M00_Basic.sh
    /scripts/chipseq/ChIProfiler_M01_Stat_PE.sh
    /scripts/chipseq/ChIProfiler_M03_CallPeaks.sh

MOUNT POINTS:
    /workspace      Working directory
    /data           Input data
    /outputs        Analysis outputs
    /references     Reference genomes
    /logs           Log files

EOF
        ;;
    *)
        exec "$@"
        ;;
esac
