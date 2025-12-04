#!/bin/bash
# ==============================================================================
# Entrypoint script for ChIP-seq & Single-Cell Pipeline container
# ==============================================================================
set -euo pipefail

# Activate virtual environment
source /opt/venv/bin/activate

# Handle commands
case "${1:-bash}" in
    tunnel|code-tunnel)
        echo "Starting VS Code Remote Tunnel..."
        echo "Follow the authentication prompts to connect."
        exec code tunnel --accept-server-license-terms
        ;;
    jupyter|lab)
        echo "Starting JupyterLab..."
        exec jupyter lab --ip=0.0.0.0 --port=8888 --no-browser --allow-root
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
    help            Show this help message

EXAMPLES:
    # Interactive shell
    docker run -it -v /data:/data <image>

    # VS Code Remote Tunnel
    docker run -it <image> tunnel

    # JupyterLab
    docker run -p 8888:8888 <image> jupyter

AVAILABLE TOOLS:
    ChIP-seq:       trim_galore, bowtie2, samtools, findPeaks, makeTagDirectory,
                    makeUCSCfile, igvtools, fastqc, cutadapt, multiqc, deeptools
    Single-cell:    scanpy, anndata, scvi-tools, Seurat
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
