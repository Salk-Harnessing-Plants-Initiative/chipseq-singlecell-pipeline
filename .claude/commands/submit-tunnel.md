# Submit VS Code Tunnel Workspace

Submit an interactive VS Code Remote Tunnel workspace to the RunAI cluster.

## Quick Command

```bash
runai workspace submit vscode-<USERNAME> \
  --project talmo-lab \
  --image ghcr.io/salk-harnessing-plants-initiative/chipseq-singlecell-pipeline:latest \
  --cpu-core-request 8 \
  --cpu-memory-request 32G \
  --host-path path=<DATA_PATH>,mount=/data,mount-propagation=HostToContainer \
  --host-path path=<OUTPUT_PATH>,mount=/outputs,mount-propagation=HostToContainer,readwrite \
  --host-path path=<WORKSPACE_PATH>,mount=/workspace,mount-propagation=HostToContainer,readwrite \
  --command -- tunnel
```

## Parameters to Replace

| Placeholder | Description | Example |
|-------------|-------------|---------|
| `<USERNAME>` | Your username | `charlotte`, `eberrigan` |
| `<DATA_PATH>` | Path to your data | `/hpi/hpi_dev/users/eberrigan/data` |
| `<OUTPUT_PATH>` | Path for outputs | `/hpi/hpi_dev/users/eberrigan/outputs` |
| `<WORKSPACE_PATH>` | Path to scripts/notebooks | `/hpi/hpi_dev/users/eberrigan/workspace` |

## Example: Interactive Analysis Session

```bash
runai workspace submit vscode-charlotte \
  --project talmo-lab \
  --image ghcr.io/salk-harnessing-plants-initiative/chipseq-singlecell-pipeline:latest \
  --cpu-core-request 16 \
  --cpu-memory-request 64G \
  --host-path path=/hpi/hpi_dev/users/charlotte/data,mount=/data,mount-propagation=HostToContainer \
  --host-path path=/hpi/hpi_dev/users/charlotte/outputs,mount=/outputs,mount-propagation=HostToContainer,readwrite \
  --host-path path=/hpi/hpi_dev/users/charlotte/workspace,mount=/workspace,mount-propagation=HostToContainer,readwrite \
  --command -- tunnel
```

## With GPU (for scvi-tools)

```bash
runai workspace submit vscode-charlotte-gpu \
  --project talmo-lab \
  --image ghcr.io/salk-harnessing-plants-initiative/chipseq-singlecell-pipeline:latest \
  --cpu-core-request 16 \
  --cpu-memory-request 128G \
  --gpu-devices-request 1 \
  --host-path path=/hpi/hpi_dev/users/charlotte/data,mount=/data,mount-propagation=HostToContainer \
  --host-path path=/hpi/hpi_dev/users/charlotte/outputs,mount=/outputs,mount-propagation=HostToContainer,readwrite \
  --host-path path=/hpi/hpi_dev/users/charlotte/workspace,mount=/workspace,mount-propagation=HostToContainer,readwrite \
  --command -- tunnel
```

## Connection Steps

### 1. Submit the workspace
```bash
runai workspace submit vscode-charlotte \
  --project talmo-lab \
  ...
```

### 2. Monitor logs for authentication URL
```bash
runai workspace logs vscode-charlotte -p talmo-lab --follow
```

Look for output like:
```
To grant access to the server, please log into https://github.com/login/device
and use code XXXX-XXXX
```

### 3. Authenticate with GitHub
1. Open the URL in your browser
2. Enter the code shown in the logs
3. Authorize VS Code

### 4. Connect from VS Code
1. Open VS Code on your local machine
2. Install "Remote - Tunnels" extension if not installed
3. Click Remote Explorer in the sidebar
4. Select "Tunnels" from dropdown
5. Find your tunnel (named after the workspace)
6. Click to connect

## Resource Recommendations

| Use Case | CPU | Memory | GPU |
|----------|-----|--------|-----|
| Light browsing/editing | 4 | 16G | No |
| Python/R analysis | 8 | 32G | No |
| Single-cell analysis | 16 | 64G | Optional |
| Large datasets (1M+ cells) | 32 | 128G | Yes |

## Available Tools in Tunnel

Once connected, you have access to:

**Python Environment:**
- scanpy, anndata, scvi-tools
- pandas, numpy, scipy
- matplotlib, seaborn, plotly
- JupyterLab

**R Environment:**
- Seurat, SingleCellExperiment
- ggplot2, dplyr, tidyr
- Bioconductor packages

**ChIP-seq Tools:**
- bowtie2, samtools
- HOMER (findPeaks, etc.)
- trim_galore, FastQC
- IGVtools

## Monitor Workspace

```bash
# Check if running
runai workspace list | grep vscode

# View status
runai workspace describe vscode-charlotte -p talmo-lab

# Follow logs
runai workspace logs vscode-charlotte -p talmo-lab --follow
```

## Execute Commands in Workspace

```bash
# Open interactive shell
runai workspace exec vscode-charlotte -it -- bash

# Run quick command
runai workspace exec vscode-charlotte -- python --version
runai workspace exec vscode-charlotte -- R --version
```

## Cleanup When Done

**Important**: Delete your tunnel workspace when finished to free cluster resources.

```bash
runai workspace delete vscode-charlotte -p talmo-lab
```

## Troubleshooting

### Tunnel not appearing in VS Code
1. Verify workspace is running: `runai workspace list | grep vscode`
2. Check logs for errors: `runai workspace logs vscode-charlotte -p talmo-lab`
3. Ensure GitHub authentication completed
4. Try refreshing Remote Explorer in VS Code

### Connection drops frequently
- Increase resource allocation
- Check network connectivity
- Workspace may have been evicted - resubmit

### Cannot write to /outputs
- Verify `readwrite` flag in host-path mount
- Check directory permissions on cluster

### Authentication URL not appearing
- Wait a few seconds for container to start
- Check if workspace is still initializing
- Try: `runai workspace describe vscode-charlotte -p talmo-lab`

## Related Commands

- `/submit-chipseq` - Submit ChIP-seq job
- `/submit-singlecell` - Submit single-cell job
- `/monitor-jobs` - Monitor all jobs
- `/cleanup-jobs` - Clean up workspaces
