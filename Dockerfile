# ==============================================================================
# ChIP-seq & Single-Cell Analysis Pipeline Container
# ==============================================================================
# Multi-stage build for bioinformatics environment
# Includes: ChIP-seq tools, Scanpy, Seurat, cellranger, VS Code Remote Tunnel
# Author: Elizabeth (Salk HPI)
# ==============================================================================

# ------------------------------------------------------------------------------
# Stage 1: Base system with build tools
# ------------------------------------------------------------------------------
FROM ubuntu:22.04 AS base

ENV DEBIAN_FRONTEND=noninteractive \
    TZ=America/Los_Angeles \
    LANG=C.UTF-8 \
    LC_ALL=C.UTF-8

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    # Build essentials
    build-essential \
    gfortran \
    pkg-config \
    cmake \
    # Compression and archives
    zip \
    unzip \
    gzip \
    bzip2 \
    libbz2-dev \
    liblzma-dev \
    zlib1g-dev \
    # Network and download
    curl \
    wget \
    ca-certificates \
    # Version control
    git \
    # Python build deps
    libssl-dev \
    libffi-dev \
    libncurses5-dev \
    libncursesw5-dev \
    libreadline-dev \
    libsqlite3-dev \
    # XML and web
    libxml2-dev \
    libcurl4-openssl-dev \
    # Graphics
    libpng-dev \
    libjpeg-dev \
    libtiff-dev \
    libfreetype6-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    # HDF5 for single-cell
    libhdf5-dev \
    # Linear algebra
    libopenblas-dev \
    liblapack-dev \
    # Java (for igvtools)
    default-jre \
    # Utilities
    parallel \
    pigz \
    vim \
    less \
    htop \
    procps \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# ------------------------------------------------------------------------------
# Stage 2: ChIP-seq bioinformatics tools
# ------------------------------------------------------------------------------
FROM base AS chipseq-tools

WORKDIR /opt

# Install Bowtie2
ARG BOWTIE2_VERSION=2.5.2
RUN curl -fsSL "https://github.com/BenLangmead/bowtie2/releases/download/v${BOWTIE2_VERSION}/bowtie2-${BOWTIE2_VERSION}-linux-x86_64.zip" -o bowtie2.zip \
    && unzip bowtie2.zip \
    && mv bowtie2-${BOWTIE2_VERSION}-linux-x86_64 /opt/bowtie2 \
    && rm bowtie2.zip
ENV PATH="/opt/bowtie2:${PATH}"

# Install Samtools
ARG SAMTOOLS_VERSION=1.19
RUN curl -fsSL "https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2" -o samtools.tar.bz2 \
    && tar -xjf samtools.tar.bz2 \
    && cd samtools-${SAMTOOLS_VERSION} \
    && ./configure --prefix=/opt/samtools \
    && make -j$(nproc) \
    && make install \
    && cd .. \
    && rm -rf samtools-${SAMTOOLS_VERSION} samtools.tar.bz2
ENV PATH="/opt/samtools/bin:${PATH}"

# Install HOMER
RUN mkdir -p /opt/homer \
    && cd /opt/homer \
    && curl -fsSL http://homer.ucsd.edu/homer/configureHomer.pl -o configureHomer.pl \
    && perl configureHomer.pl -install homer \
    && perl configureHomer.pl -install hg38 \
    && perl configureHomer.pl -install mm10
ENV PATH="/opt/homer/bin:${PATH}"

# Install IGVtools
ARG IGVTOOLS_VERSION=2.16.2
RUN curl -fsSL "https://data.broadinstitute.org/igv/projects/downloads/2.16/IGV_${IGVTOOLS_VERSION}.zip" -o igv.zip \
    && unzip igv.zip \
    && mv IGV_${IGVTOOLS_VERSION} /opt/igv \
    && rm igv.zip \
    && ln -s /opt/igv/igvtools /usr/local/bin/igvtools
ENV PATH="/opt/igv:${PATH}"

# Install FastQC
ARG FASTQC_VERSION=0.12.1
RUN curl -fsSL "https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v${FASTQC_VERSION}.zip" -o fastqc.zip \
    && unzip fastqc.zip \
    && chmod +x FastQC/fastqc \
    && mv FastQC /opt/fastqc \
    && rm fastqc.zip \
    && ln -s /opt/fastqc/fastqc /usr/local/bin/fastqc
ENV PATH="/opt/fastqc:${PATH}"

# Install Cutadapt and TrimGalore (Python-based, install later with uv)

# ------------------------------------------------------------------------------
# Stage 3: Python environment with uv
# ------------------------------------------------------------------------------
FROM chipseq-tools AS python-env

# Install uv
RUN curl -LsSf https://astral.sh/uv/install.sh | env CARGO_HOME=/opt/uv UV_INSTALL_DIR=/opt/uv sh
ENV PATH="/opt/uv:/opt/uv/bin:/root/.local/bin:${PATH}"

# Set up Python environment with uv
WORKDIR /opt/python-env
COPY pyproject.toml .

# Install Python packages using uv sync (best practice)
# Dependencies defined in pyproject.toml
RUN uv sync --python 3.11 \
    && rm -rf /root/.cache/uv

ENV VIRTUAL_ENV="/opt/python-env/.venv" \
    PATH="/opt/python-env/.venv/bin:${PATH}"

# Install TrimGalore
RUN curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/refs/tags/0.6.10.tar.gz -o trimgalore.tar.gz \
    && tar -xzf trimgalore.tar.gz \
    && mv TrimGalore-0.6.10/trim_galore /usr/local/bin/ \
    && chmod +x /usr/local/bin/trim_galore \
    && rm -rf TrimGalore-0.6.10 trimgalore.tar.gz

# ------------------------------------------------------------------------------
# Stage 4: R environment with Seurat
# ------------------------------------------------------------------------------
FROM python-env AS r-env

# Install R 4.3
RUN apt-get update && apt-get install -y --no-install-recommends \
    software-properties-common \
    dirmngr \
    gnupg \
    && wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc \
    && add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/" \
    && apt-get update \
    && apt-get install -y --no-install-recommends \
    r-base \
    r-base-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install R packages
RUN R -e "install.packages(c( \
    'devtools', \
    'BiocManager', \
    'Seurat', \
    'SeuratObject', \
    'hdf5r', \
    'ggplot2', \
    'dplyr', \
    'tidyr', \
    'patchwork', \
    'cowplot', \
    'Matrix', \
    'data.table', \
    'future', \
    'reticulate' \
    ), repos='https://cloud.r-project.org', Ncpus=$(nproc))"

# Install Bioconductor packages for single-cell
RUN R -e "BiocManager::install(c( \
    'SingleCellExperiment', \
    'scater', \
    'scran', \
    'DropletUtils', \
    'celldex', \
    'SingleR', \
    'DESeq2', \
    'edgeR', \
    'limma' \
    ), Ncpus=$(nproc))"

# ------------------------------------------------------------------------------
# Stage 5: cellranger
# ------------------------------------------------------------------------------
FROM r-env AS cellranger-env

# Note: cellranger requires accepting 10x Genomics license
# Download link changes - users may need to update this or mount cellranger
ARG CELLRANGER_VERSION=8.0.1
# cellranger must be downloaded manually from 10x Genomics due to license
# This creates a placeholder - users should mount or install separately
RUN mkdir -p /opt/cellranger \
    && echo "cellranger ${CELLRANGER_VERSION} placeholder - mount or install separately" > /opt/cellranger/README.txt
ENV PATH="/opt/cellranger:${PATH}"

# Alternative: If cellranger tarball is available locally during build:
# COPY cellranger-${CELLRANGER_VERSION}.tar.gz /tmp/
# RUN tar -xzf /tmp/cellranger-${CELLRANGER_VERSION}.tar.gz -C /opt/ \
#     && rm /tmp/cellranger-${CELLRANGER_VERSION}.tar.gz
# ENV PATH="/opt/cellranger-${CELLRANGER_VERSION}:${PATH}"

# ------------------------------------------------------------------------------
# Stage 6: Runtime with VS Code support
# ------------------------------------------------------------------------------
FROM cellranger-env AS runtime

# Install VS Code CLI for Remote Tunnel
# Use the direct download URL which is more stable
RUN curl -fsSL "https://update.code.visualstudio.com/latest/cli-linux-x64/stable" -o /tmp/vscode_cli.tar.gz \
    && tar -xzf /tmp/vscode_cli.tar.gz -C /usr/local/bin \
    && rm /tmp/vscode_cli.tar.gz

# Create workspace directories
RUN mkdir -p /workspace /data /outputs /logs /references

# Copy scripts
COPY scripts/ /scripts/
RUN chmod +x /scripts/*.sh /scripts/**/*.sh 2>/dev/null || true

# Set working directory
WORKDIR /workspace

# Copy entrypoint
COPY scripts/entrypoint.sh /entrypoint.sh
RUN chmod +x /entrypoint.sh

ENTRYPOINT ["/entrypoint.sh"]
CMD ["bash"]

# ==============================================================================
# Metadata
# ==============================================================================
LABEL org.opencontainers.image.title="ChIP-seq & Single-Cell Pipeline" \
      org.opencontainers.image.description="Bioinformatics container for ChIP-seq and single-cell RNA-seq analysis" \
      org.opencontainers.image.authors="Salk Harnessing Plants Initiative" \
      org.opencontainers.image.source="https://github.com/Salk-Harnessing-Plants-Initiative/chipseq-singlecell-pipeline" \
      org.opencontainers.image.version="1.0.0" \
      python.version="3.11" \
      r.version="4.3+" \
      uv.managed="true"
