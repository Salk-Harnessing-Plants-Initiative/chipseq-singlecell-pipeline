# Monitor RunAI Jobs

Monitor RunAI workspaces with status updates and filtering.

## Quick Status Check

```bash
# List all workspaces
runai workspace list

# Filter by type
runai workspace list | grep chipseq
runai workspace list | grep singlecell
runai workspace list | grep vscode
```

## Detailed Status

```bash
# Show detailed workspace info
runai workspace describe <WORKSPACE_NAME> -p talmo-lab
```

## View Logs

```bash
# Follow logs in real-time
runai workspace logs <WORKSPACE_NAME> -p talmo-lab --follow

# Get last N lines
runai workspace logs <WORKSPACE_NAME> -p talmo-lab --tail 100
```

## Using Monitor Script

```bash
# Interactive dashboard with auto-refresh
./scripts/monitor-runai-jobs.sh

# Filter to specific prefix
./scripts/monitor-runai-jobs.sh --prefix chipseq

# Export status to JSON
./scripts/monitor-runai-jobs.sh --export status.json

# Single status check (no refresh)
./scripts/monitor-runai-jobs.sh --once
```

## Dashboard Output

```
═══════════════════════════════════════════════════════════════
ChIP-seq & Single-Cell Pipeline - Job Monitor
Updated: 2025-01-15 14:23:45
═══════════════════════════════════════════════════════════════

SUMMARY:
  Running:    3 jobs
  Succeeded:  12 jobs
  Failed:     1 job
  Pending:    0 jobs
  Total:      16 jobs

RUNNING JOBS:
  NAME                    STATUS    DURATION  MEMORY    CPU
  chipseq-col0-h3k4me3    Running   45m       48/64G    14/16
  singlecell-root-atlas   Running   2h 15m    96/128G   28/32
  vscode-charlotte        Running   3h 30m    24/32G    6/8

RECENTLY COMPLETED:
  chipseq-col0-h3k27me3   Succeeded  1h 23m ago
  singlecell-leaf         Succeeded  2h 05m ago

FAILED JOBS:
  chipseq-test-001        Failed    OOMKilled

Press Ctrl+C to exit. Refreshing every 30s...
```

## Filter by Status

```bash
# Show only running jobs
runai workspace list | grep -E "Running"

# Show only failed jobs
runai workspace list | grep -E "Failed"

# Show only succeeded jobs
runai workspace list | grep -E "Succeeded"

# Count by status
runai workspace list | awk '{print $2}' | sort | uniq -c
```

## Check Specific Job Types

### ChIP-seq Jobs
```bash
runai workspace list | grep chipseq
runai workspace describe chipseq-col0-h3k4me3 -p talmo-lab
runai workspace logs chipseq-col0-h3k4me3 -p talmo-lab --follow
```

### Single-Cell Jobs
```bash
runai workspace list | grep singlecell
runai workspace describe singlecell-root-atlas -p talmo-lab
runai workspace logs singlecell-root-atlas -p talmo-lab --follow
```

### VS Code Tunnels
```bash
runai workspace list | grep vscode
runai workspace describe vscode-charlotte -p talmo-lab
runai workspace logs vscode-charlotte -p talmo-lab --follow
```

## Resource Usage

```bash
# Show resource allocation
runai workspace list -o wide

# Check specific workspace resources
runai workspace describe <NAME> -p talmo-lab | grep -E "CPU|Memory|GPU"
```

## Job Timeline

Track when jobs started and finished:

```bash
# Show workspace events
runai workspace describe <NAME> -p talmo-lab
```

## Wait for Job Completion

```bash
# Simple wait loop
while runai workspace list | grep -q "chipseq-col0.*Running"; do
  echo "Job still running..."
  sleep 60
done
echo "Job completed!"

# With notification
while runai workspace list | grep -q "<JOB_NAME>.*Running"; do
  sleep 60
done && echo "Job done!" | mail -s "RunAI Job Complete" you@email.com
```

## Troubleshooting

### Job stuck in "Pending"
```bash
# Check detailed status
runai workspace describe <NAME> -p talmo-lab

# Common causes:
# - Insufficient cluster resources
# - Invalid image
# - Host path issues
```

### Job shows "Failed"
```bash
# Check logs for error
runai workspace logs <NAME> -p talmo-lab --tail 200

# Common causes:
# - OOMKilled (out of memory)
# - Script error
# - Missing files
```

### Cannot see job in list
```bash
# Check if job exists
runai workspace describe <NAME> -p talmo-lab

# Job may have been deleted or never created
# Check for typos in job name
```

## Related Commands

- `/submit-chipseq` - Submit ChIP-seq job
- `/submit-singlecell` - Submit single-cell job
- `/submit-tunnel` - Submit VS Code tunnel
- `/cleanup-jobs` - Clean up completed jobs
