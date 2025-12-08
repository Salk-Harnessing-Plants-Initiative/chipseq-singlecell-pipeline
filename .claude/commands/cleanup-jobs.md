# Cleanup RunAI Jobs

Safely clean up completed, failed, or orphaned RunAI workspaces.

**WARNING**: This command deletes jobs and associated resources. Always verify what will be deleted before proceeding.

## Safe Command (Recommended)

```bash
# Interactive cleanup with confirmation
./scripts/cleanup-runai.sh
```

## Safety Guardrails

### 1. Always Preview First (Dry Run)

```bash
# See what WOULD be deleted without actually deleting
./scripts/cleanup-runai.sh --dry-run

# Output shows:
# Would delete:
#   chipseq-col0-test (Succeeded, completed 2h ago)
#   singlecell-test-001 (Failed, OOMKilled)
# Total: 2 jobs (0 running, 1 completed, 1 failed)
```

### 2. Never Delete Running Jobs

```bash
# The script ALWAYS excludes running jobs by default
./scripts/cleanup-runai.sh --all

# Running jobs are protected unless --force-running is used
```

### 3. Clean by Category Only

```bash
# Clean ONLY completed jobs (safest)
./scripts/cleanup-runai.sh --completed

# Clean ONLY failed jobs
./scripts/cleanup-runai.sh --failed

# Avoid --all unless you're certain
```

### 4. Age-Based Cleanup

```bash
# Only delete jobs older than 24 hours
./scripts/cleanup-runai.sh --older-than 24h

# Only delete jobs older than 7 days
./scripts/cleanup-runai.sh --older-than 7d
```

## Command Options

| Option | Description |
|--------|-------------|
| `--dry-run` | Preview without deleting |
| `--completed` | Delete only succeeded jobs |
| `--failed` | Delete only failed jobs |
| `--all` | Delete both completed and failed |
| `--older-than <time>` | Only jobs older than threshold |
| `--prefix <name>` | Only jobs matching prefix |
| `--yes` | Skip confirmation prompt |
| `--force-running` | Include running jobs (dangerous!) |

## Interactive Mode (Default)

```bash
./scripts/cleanup-runai.sh

# Prompts:
# Found 12 completed jobs
# Found 3 failed jobs
# Found 2 running jobs (protected)
#
# What would you like to clean up?
# 1) Completed jobs only
# 2) Failed jobs only
# 3) Both completed and failed
# 4) Show detailed list first
# 5) Cancel
# Select option:
```

## Quick Cleanup Commands

### Delete All Completed Jobs
```bash
./scripts/cleanup-runai.sh --completed --yes
```

### Delete All Failed Jobs
```bash
./scripts/cleanup-runai.sh --failed --yes
```

### Delete Specific Job Type
```bash
# Delete all completed ChIP-seq jobs
./scripts/cleanup-runai.sh --completed --prefix chipseq --yes

# Delete all completed single-cell jobs
./scripts/cleanup-runai.sh --completed --prefix singlecell --yes
```

### Delete Old Jobs Only
```bash
# Delete completed jobs older than 1 day
./scripts/cleanup-runai.sh --completed --older-than 24h --yes

# Delete all jobs older than 1 week
./scripts/cleanup-runai.sh --all --older-than 7d --yes
```

## Manual Cleanup (Without Script)

### Delete Single Workspace
```bash
runai workspace delete <WORKSPACE_NAME> -p talmo-lab
```

### Delete Multiple by Filter
```bash
# List what would be deleted first
runai workspace list | grep "chipseq.*Succeeded"

# Then delete (carefully!)
runai workspace list | grep "chipseq.*Succeeded" | awk '{print $1}' | \
  xargs -I {} runai workspace delete {} -p talmo-lab
```

### Delete Failed Jobs Only
```bash
runai workspace list | grep "Failed" | awk '{print $1}' | \
  xargs -I {} runai workspace delete {} -p talmo-lab
```

## Step-by-Step Safe Cleanup Workflow

### Step 1: Check Status
```bash
# See what's running vs completed
runai workspace list

# Count by status
runai workspace list | awk '{print $2}' | sort | uniq -c
```

### Step 2: Verify Outputs Saved
```bash
# Check output directory has expected files
ls /hpi/path/to/outputs/

# For ChIP-seq, look for:
# - *_sorted_rmdup_uniq.bam
# - *_TagDir/
# - *.bedGraph.gz
# - 1_res_read_mapping_stat.txt

# For single-cell, look for:
# - *.h5ad or *.rds
# - plots/figures
```

### Step 3: Dry Run
```bash
./scripts/cleanup-runai.sh --completed --dry-run
```

### Step 4: Execute Cleanup
```bash
./scripts/cleanup-runai.sh --completed --yes
```

### Step 5: Verify
```bash
# Check remaining jobs
runai workspace list

# Should only see running jobs (if any)
```

## Handle Failed Jobs

Before deleting failed jobs, investigate:

```bash
# List failed jobs
runai workspace list | grep Failed

# Check why it failed
runai workspace logs <FAILED_JOB_NAME> -p talmo-lab --tail 100

# Common failures:
# - OOMKilled → Resubmit with more memory
# - Script error → Fix script and resubmit
# - Missing file → Check input paths
```

## Cleanup Protection Features

The cleanup script includes:

1. **No silent deletion** - Always shows what will be deleted
2. **Running job protection** - Never deletes running jobs by default
3. **Confirmation prompts** - Requires explicit confirmation
4. **Dry-run mode** - Test deletions without executing
5. **Age filters** - Only delete jobs older than threshold
6. **Status filters** - Only delete specific status types
7. **Detailed logging** - Records all deletions

## When NOT to Cleanup

Do NOT cleanup if:
- Jobs are still running
- Outputs not yet verified/saved
- Investigating failures
- Within 24h of submission (give time to review)

Safe to cleanup when:
- All jobs completed or failed
- Outputs saved and verified
- Logs reviewed for failed jobs
- Disk space needed

## Troubleshooting

### "Permission denied" errors
```bash
# Check you have delete permissions
runai workspace list  # If this works, you should have permissions
```

### Job won't delete (stuck)
```bash
# Force delete (last resort)
runai workspace delete <NAME> -p talmo-lab --force
```

### Accidentally deleted important job
```bash
# Output files should still exist on cluster
ls /hpi/path/to/outputs/

# Logs may still be retrievable via kubectl
kubectl logs <pod-name> -n runai-talmo-lab

# Resubmit job if needed
/submit-chipseq  # or /submit-singlecell
```

## Related Commands

- `/monitor-jobs` - Check job status before cleanup
- `/submit-chipseq` - Resubmit ChIP-seq jobs
- `/submit-singlecell` - Resubmit single-cell jobs
