# RunAI Quick Reference for ChIP-seq & Single-Cell Pipeline

Quick command reference for running analysis jobs on the Salk HPI RunAI cluster.

## Project Info

- **Project name**: `talmo-lab`
- **Kubernetes namespace**: `runai-talmo-lab`
- **Image**: `ghcr.io/salk-harnessing-plants-initiative/chipseq-singlecell-pipeline:latest`

## Quick Commands

### Submit ChIP-seq Analysis Job

```bash
runai workspace submit chipseq-<SAMPLE> \
  --project talmo-lab \
  --image ghcr.io/salk-harnessing-plants-initiative/chipseq-singlecell-pipeline:latest \
  --cpu-core-request 16 \
  --cpu-memory-request 64G \
  --host-path path=<DATA_PATH>,mount=/data,mount-propagation=HostToContainer \
  --host-path path=<OUTPUT_PATH>,mount=/outputs,mount-propagation=HostToContainer,readwrite \
  --host-path path=<REFERENCE_PATH>,mount=/references,mount-propagation=HostToContainer \
  --environment REFERENCE=/references/<INDEX_NAME> \
  --environment CPU_CORES=16 \
  --environment READ_TYPE=PE \
  --environment IP_SAMPLE=<IP_SAMPLE_NAME> \
  --environment INPUT_SAMPLE=<INPUT_SAMPLE_NAME> \
  --command -- /scripts/chipseq/ChIProfiler_M00_Basic.sh
```

### Submit Single-Cell Analysis (Python/Scanpy)

```bash
runai workspace submit singlecell-<NAME> \
  --project talmo-lab \
  --image ghcr.io/salk-harnessing-plants-initiative/chipseq-singlecell-pipeline:latest \
  --cpu-core-request 32 \
  --cpu-memory-request 128G \
  --gpu-devices-request 1 \
  --host-path path=<DATA_PATH>,mount=/data,mount-propagation=HostToContainer \
  --host-path path=<OUTPUT_PATH>,mount=/outputs,mount-propagation=HostToContainer,readwrite \
  --command -- python /workspace/<SCRIPT>.py
```

### Submit Single-Cell Analysis (R/Seurat)

```bash
runai workspace submit singlecell-<NAME> \
  --project talmo-lab \
  --image ghcr.io/salk-harnessing-plants-initiative/chipseq-singlecell-pipeline:latest \
  --cpu-core-request 32 \
  --cpu-memory-request 128G \
  --host-path path=<DATA_PATH>,mount=/data,mount-propagation=HostToContainer \
  --host-path path=<OUTPUT_PATH>,mount=/outputs,mount-propagation=HostToContainer,readwrite \
  --command -- Rscript /workspace/<SCRIPT>.R
```

### Submit VS Code Tunnel

```bash
runai workspace submit vscode-<USERNAME> \
  --project talmo-lab \
  --image ghcr.io/salk-harnessing-plants-initiative/chipseq-singlecell-pipeline:latest \
  --cpu-core-request 8 \
  --cpu-memory-request 32G \
  --host-path path=<DATA_PATH>,mount=/data,mount-propagation=HostToContainer \
  --host-path path=<OUTPUT_PATH>,mount=/outputs,mount-propagation=HostToContainer,readwrite \
  --command -- tunnel
```

### List All Workloads

```bash
runai workspace list
```

Filter by type:
```bash
runai workspace list | grep chipseq
runai workspace list | grep singlecell
runai workspace list | grep vscode
```

### Check Specific Workspace

```bash
runai workspace describe <WORKSPACE_NAME> -p talmo-lab
```

### View Logs

```bash
runai workspace logs <WORKSPACE_NAME> -p talmo-lab --follow
```

### Delete Workspace

```bash
runai workspace delete <WORKSPACE_NAME> -p talmo-lab
```

## Common Workflows

### Test ChIP-seq Pipeline

```bash
# Submit test job
runai workspace submit chipseq-test \
  --project talmo-lab \
  --image ghcr.io/salk-harnessing-plants-initiative/chipseq-singlecell-pipeline:latest \
  --cpu-core-request 16 \
  --cpu-memory-request 64G \
  --host-path path=/hpi/path/to/data,mount=/data,mount-propagation=HostToContainer \
  --host-path path=/hpi/path/to/outputs,mount=/outputs,mount-propagation=HostToContainer,readwrite \
  --host-path path=/hpi/path/to/references,mount=/references,mount-propagation=HostToContainer \
  --environment REFERENCE=/references/TAIR10 \
  --environment CPU_CORES=16 \
  --environment READ_TYPE=PE \
  --environment IP_SAMPLE=IP_col0 \
  --environment INPUT_SAMPLE=input_col0 \
  --command -- /scripts/chipseq/ChIProfiler_M00_Basic.sh

# Monitor
runai workspace list | grep chipseq-test

# View logs
runai workspace logs chipseq-test -p talmo-lab --follow

# Clean up
runai workspace delete chipseq-test -p talmo-lab
```

### Interactive Analysis Session

```bash
# Submit VS Code tunnel
runai workspace submit vscode-myname \
  --project talmo-lab \
  --image ghcr.io/salk-harnessing-plants-initiative/chipseq-singlecell-pipeline:latest \
  --cpu-core-request 16 \
  --cpu-memory-request 64G \
  --host-path path=/hpi/path/to/data,mount=/data,mount-propagation=HostToContainer \
  --host-path path=/hpi/path/to/outputs,mount=/outputs,mount-propagation=HostToContainer,readwrite \
  --command -- tunnel

# Get authentication URL from logs
runai workspace logs vscode-myname -p talmo-lab --follow

# Connect via VS Code Remote Explorer → Tunnels
# Clean up when done
runai workspace delete vscode-myname -p talmo-lab
```

### Monitor Jobs

```bash
# Use monitoring script
./scripts/monitor-runai-jobs.sh

# Or watch command
watch -n 10 'runai workspace list'

# Count by status
runai workspace list | awk '{print $2}' | sort | uniq -c
```

### Cleanup Completed Jobs

```bash
# Preview what will be deleted
./scripts/cleanup-runai.sh --completed --dry-run

# Delete completed jobs
./scripts/cleanup-runai.sh --completed --yes

# Delete failed jobs
./scripts/cleanup-runai.sh --failed --yes
```

## Troubleshooting

### Job stuck in "Pending"

```bash
# Get detailed status
runai workspace describe <WORKSPACE_NAME> -p talmo-lab

# Common causes:
# - Insufficient cluster resources
# - Invalid image name
# - Host path doesn't exist
```

### Job fails immediately

```bash
# Check logs for errors
runai workspace logs <WORKSPACE_NAME> -p talmo-lab

# Common causes:
# - Image pull error (wrong tag)
# - Permission denied on paths
# - Script error
```

### Out of memory (OOMKilled)

Increase memory:
```bash
--cpu-memory-request 128G  # Instead of 64G
```

### VS Code tunnel not connecting

1. Check logs: `runai workspace logs vscode-<NAME> -p talmo-lab`
2. Look for GitHub authentication URL
3. Complete device authentication
4. Refresh VS Code Remote Explorer

## Helper Scripts

```bash
# Monitor jobs with dashboard
./scripts/monitor-runai-jobs.sh --watch

# Safe cleanup with confirmation
./scripts/cleanup-runai.sh --all
```

## Key Differences from Old RunAI CLI

| Old Command | New Command |
|-------------|-------------|
| `runai submit` | `runai workspace submit` |
| `runai list jobs` | `runai workspace list` |
| `runai describe job` | `runai workspace describe` |
| `runai logs` | `runai workspace logs` |
| `runai delete job` | `runai workspace delete` |
| `--cpu 12` | `--cpu-core-request 12` |
| `--memory 32G` | `--cpu-memory-request 32G` |
| `--gpu 1` | `--gpu-devices-request 1` |
| `--host-path /path:/mount:ro` | `--host-path path=/path,mount=/mount,mount-propagation=HostToContainer` |
| `--project runai-talmo-lab` | `--project talmo-lab` |

## Additional Resources

- Full guide: [MANUAL_RUNAI_EXECUTION.md](MANUAL_RUNAI_EXECUTION.md)
- RunAI docs: https://docs.run.ai/
- Salk RunAI guide: https://researchit.salk.edu/runai/
