# RunAI CLI v2 Skill for ChIP-seq & Single-Cell Pipeline

This skill provides comprehensive guidance for running ChIP-seq and single-cell analysis jobs on the Salk HPI RunAI cluster.

## Project Configuration

- **Project**: `talmo-lab`
- **Namespace**: `runai-talmo-lab`
- **Image**: `ghcr.io/salk-harnessing-plants-initiative/chipseq-singlecell-pipeline:latest`

## RunAI CLI v2 Command Reference

### Core Commands

| Command | Description |
|---------|-------------|
| `runai workspace submit <name>` | Submit a new workspace/job |
| `runai workspace list` | List all workspaces |
| `runai workspace describe <name>` | Show workspace details |
| `runai workspace logs <name>` | View workspace logs |
| `runai workspace delete <name>` | Delete a workspace |
| `runai workspace exec <name>` | Execute command in workspace |

### Common Parameters

| Parameter | Description | Example |
|-----------|-------------|---------|
| `--project` | RunAI project name | `--project talmo-lab` |
| `--image` | Container image | `--image ghcr.io/.../chipseq-singlecell-pipeline:latest` |
| `--cpu-core-request` | CPU cores | `--cpu-core-request 16` |
| `--cpu-memory-request` | Memory | `--cpu-memory-request 64G` |
| `--gpu-devices-request` | GPUs | `--gpu-devices-request 1` |
| `--host-path` | Mount host directory | `--host-path path=/src,mount=/dst,mount-propagation=HostToContainer` |
| `--environment` | Environment variable | `--environment VAR=value` |
| `--command` | Container command | `--command -- /scripts/entrypoint.sh` |

---

## Job Templates

### 1. ChIP-seq Analysis Job

Submit a complete ChIP-seq analysis pipeline:

```bash
runai workspace submit chipseq-<SAMPLE_NAME> \
  --project talmo-lab \
  --image ghcr.io/salk-harnessing-plants-initiative/chipseq-singlecell-pipeline:latest \
  --cpu-core-request 16 \
  --cpu-memory-request 64G \
  --host-path path=/hpi/path/to/data,mount=/data,mount-propagation=HostToContainer \
  --host-path path=/hpi/path/to/outputs,mount=/outputs,mount-propagation=HostToContainer,readwrite \
  --host-path path=/hpi/path/to/references,mount=/references,mount-propagation=HostToContainer \
  --environment REFERENCE=/references/bowtie2_index \
  --environment CPU_CORES=16 \
  --environment READ_TYPE=PE \
  --environment IP_SAMPLE=<IP_SAMPLE_NAME> \
  --environment INPUT_SAMPLE=<INPUT_SAMPLE_NAME> \
  --command -- /scripts/chipseq/ChIProfiler_M00_Basic.sh
```

**Environment Variables for ChIP-seq:**
| Variable | Description | Default |
|----------|-------------|---------|
| `REFERENCE` | Path to bowtie2 index | `/references/bowtie2_index` |
| `CPU_CORES` | Number of parallel threads | `12` |
| `READ_TYPE` | `PE` (paired-end) or `SE` (single-end) | `PE` |
| `IP_SAMPLE` | IP sample file prefix | Required |
| `INPUT_SAMPLE` | Input/control sample prefix | Required |
| `TAG_SUFFIX` | Tag directory suffix | `_sorted_rmdup_uniq_TagDir` |

### 2. Single-Cell Analysis (Python/Scanpy)

Submit a Scanpy-based single-cell analysis:

```bash
runai workspace submit singlecell-scanpy-<NAME> \
  --project talmo-lab \
  --image ghcr.io/salk-harnessing-plants-initiative/chipseq-singlecell-pipeline:latest \
  --cpu-core-request 32 \
  --cpu-memory-request 128G \
  --gpu-devices-request 1 \
  --host-path path=/hpi/path/to/data,mount=/data,mount-propagation=HostToContainer \
  --host-path path=/hpi/path/to/outputs,mount=/outputs,mount-propagation=HostToContainer,readwrite \
  --command -- python /workspace/analysis.py
```

**For large datasets (3M+ cells):**
- Use `--cpu-memory-request 256G`
- Use `--gpu-devices-request 1` for scvi-tools GPU acceleration

### 3. Single-Cell Analysis (R/Seurat)

Submit a Seurat-based single-cell analysis:

```bash
runai workspace submit singlecell-seurat-<NAME> \
  --project talmo-lab \
  --image ghcr.io/salk-harnessing-plants-initiative/chipseq-singlecell-pipeline:latest \
  --cpu-core-request 32 \
  --cpu-memory-request 128G \
  --host-path path=/hpi/path/to/data,mount=/data,mount-propagation=HostToContainer \
  --host-path path=/hpi/path/to/outputs,mount=/outputs,mount-propagation=HostToContainer,readwrite \
  --command -- Rscript /workspace/analysis.R
```

### 4. VS Code Remote Tunnel

Submit an interactive workspace with VS Code tunnel:

```bash
runai workspace submit vscode-<USERNAME> \
  --project talmo-lab \
  --image ghcr.io/salk-harnessing-plants-initiative/chipseq-singlecell-pipeline:latest \
  --cpu-core-request 8 \
  --cpu-memory-request 32G \
  --host-path path=/hpi/path/to/data,mount=/data,mount-propagation=HostToContainer \
  --host-path path=/hpi/path/to/outputs,mount=/outputs,mount-propagation=HostToContainer,readwrite \
  --host-path path=/hpi/path/to/references,mount=/references,mount-propagation=HostToContainer \
  --command -- tunnel
```

**After submission:**
1. Check logs: `runai workspace logs vscode-<USERNAME> -p talmo-lab --follow`
2. Look for the GitHub authentication URL
3. Complete authentication in browser
4. Connect via VS Code: Remote Explorer → Tunnels

### 5. JupyterLab Session

Submit JupyterLab for interactive analysis:

```bash
runai workspace submit jupyter-<USERNAME> \
  --project talmo-lab \
  --image ghcr.io/salk-harnessing-plants-initiative/chipseq-singlecell-pipeline:latest \
  --cpu-core-request 16 \
  --cpu-memory-request 64G \
  --host-path path=/hpi/path/to/data,mount=/data,mount-propagation=HostToContainer \
  --host-path path=/hpi/path/to/outputs,mount=/outputs,mount-propagation=HostToContainer,readwrite \
  --command -- jupyter
```

---

## Job Management

### List Jobs

```bash
# List all workspaces
runai workspace list

# Filter to specific jobs
runai workspace list | grep chipseq
runai workspace list | grep singlecell
runai workspace list | grep vscode
```

### Check Job Status

```bash
# Detailed status
runai workspace describe <WORKSPACE_NAME> -p talmo-lab

# Quick status check
runai workspace list | grep <WORKSPACE_NAME>
```

### View Logs

```bash
# Follow logs in real-time
runai workspace logs <WORKSPACE_NAME> -p talmo-lab --follow

# Get last N lines
runai workspace logs <WORKSPACE_NAME> -p talmo-lab --tail 100
```

### Execute Commands in Running Workspace

```bash
# Interactive shell
runai workspace exec <WORKSPACE_NAME> -it -- bash

# Run specific command
runai workspace exec <WORKSPACE_NAME> -- ls /outputs
```

### Delete Workspace

```bash
# Delete single workspace
runai workspace delete <WORKSPACE_NAME> -p talmo-lab

# Delete multiple (filter first!)
runai workspace list | grep "chipseq.*Succeeded" | awk '{print $1}' | \
  xargs -I {} runai workspace delete {} -p talmo-lab
```

---

## Troubleshooting

### Job Stuck in "Pending"

```bash
# Check detailed status
runai workspace describe <WORKSPACE_NAME> -p talmo-lab

# Common causes:
# - Insufficient cluster resources
# - Invalid image name
# - Host path doesn't exist
```

### Job Fails Immediately

```bash
# Check logs for errors
runai workspace logs <WORKSPACE_NAME> -p talmo-lab

# Common causes:
# - Image pull error (check image name/tag)
# - Permission denied on host paths
# - Invalid command
```

### Out of Memory (OOMKilled)

Increase memory allocation:
```bash
--cpu-memory-request 128G  # Instead of 64G
```

### VS Code Tunnel Not Connecting

1. Check tunnel is running: `runai workspace logs vscode-<NAME> -p talmo-lab`
2. Look for "Open this link in your browser" message
3. Complete GitHub device authentication
4. Verify workspace status: `runai workspace describe vscode-<NAME> -p talmo-lab`

---

## CLI v2 Migration Reference

If you see old v1 commands, convert them:

| Old (v1) | New (v2) |
|----------|----------|
| `runai submit` | `runai workspace submit` |
| `runai list jobs` | `runai workspace list` |
| `runai describe job` | `runai workspace describe` |
| `runai logs` | `runai workspace logs` |
| `runai delete job` | `runai workspace delete` |
| `--cpu 12` | `--cpu-core-request 12` |
| `--memory 32G` | `--cpu-memory-request 32G` |
| `--gpu 1` | `--gpu-devices-request 1` |
| `--host-path /src:/dst:ro` | `--host-path path=/src,mount=/dst,mount-propagation=HostToContainer` |
| `--project runai-talmo-lab` | `--project talmo-lab` |

---

## Mount Points

The container expects these mount points:

| Mount | Purpose | Mode |
|-------|---------|------|
| `/data` | Input data (FASTQ, h5ad, etc.) | Read-only |
| `/outputs` | Analysis outputs | Read-write |
| `/references` | Reference genomes, indices | Read-only |
| `/workspace` | Working directory, scripts | Read-write |

---

## Resources

- [RunAI Documentation](https://docs.run.ai/)
- [Salk HPI RunAI Guide](https://researchit.salk.edu/runai/)
- [Container Image on GHCR](https://github.com/salk-harnessing-plants-initiative/chipseq-singlecell-pipeline/pkgs/container/chipseq-singlecell-pipeline)
