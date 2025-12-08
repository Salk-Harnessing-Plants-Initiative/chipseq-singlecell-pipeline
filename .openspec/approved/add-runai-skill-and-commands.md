# OpenSpec Proposal: Add RunAI Skill and Claude Commands

**Change ID**: `add-runai-skill-and-commands`
**Created**: 2025-12-05
**Status**: Implemented
**Approved**: 2025-12-08
**Implemented**: 2025-12-08

## Summary

Add a comprehensive RunAI CLI v2 skill and Claude slash commands for the ChIP-seq & Single-Cell Pipeline repository, enabling users to easily submit, monitor, and manage analysis jobs on the Salk HPI RunAI cluster.

## Motivation / Problem Statement

Currently, users must manually construct complex RunAI commands for:
- Submitting ChIP-seq analysis jobs
- Submitting single-cell analysis jobs
- Monitoring job status
- Managing job cleanup
- Connecting via VS Code tunnel

The GAPIT3 pipeline has successfully implemented a RunAI skill and command structure that provides:
- Documented CLI v2 syntax with tested commands
- Slash commands for common operations
- Safety guardrails for destructive operations
- Migration guidance from old RunAI CLI v1

This proposal brings the same capabilities to the ChIP-seq & Single-Cell pipeline.

## Proposed Solution

### 1. RunAI Skill (`/.claude/skills/runai.md`)

A comprehensive skill document providing:
- RunAI CLI v2 command reference
- Project-specific configuration (talmo-lab project)
- Tested command templates for ChIP-seq and single-cell jobs
- VS Code tunnel connection patterns
- Troubleshooting guidance

### 2. Claude Slash Commands

| Command | Purpose |
|---------|---------|
| `/submit-chipseq` | Submit a ChIP-seq analysis job to RunAI |
| `/submit-singlecell` | Submit a single-cell analysis job to RunAI |
| `/submit-tunnel` | Submit a VS Code tunnel workspace |
| `/monitor-jobs` | Monitor RunAI job status with dashboard |
| `/cleanup-jobs` | Safely clean up completed/failed jobs |
| `/docker-build` | Build Docker image locally |
| `/docker-test` | Test Docker container |
| `/validate-bash` | Validate shell scripts with ShellCheck |

### 3. RunAI Documentation (`/docs/`)

| Document | Purpose |
|----------|---------|
| `RUNAI_QUICK_REFERENCE.md` | Command cheat sheet with tested examples |
| `MANUAL_RUNAI_EXECUTION.md` | Step-by-step execution guide |

## Files Affected

### New Files
```
.claude/
├── skills/
│   └── runai.md                    # RunAI CLI v2 skill
├── commands/
│   ├── submit-chipseq.md           # ChIP-seq job submission
│   ├── submit-singlecell.md        # Single-cell job submission
│   ├── submit-tunnel.md            # VS Code tunnel submission
│   ├── monitor-jobs.md             # Job monitoring
│   ├── cleanup-jobs.md             # Job cleanup (with safety)
│   ├── docker-build.md             # Docker build
│   ├── docker-test.md              # Docker test
│   ├── validate-bash.md            # Bash validation
│   └── openspec/
│       ├── proposal.md             # OpenSpec proposal command
│       ├── apply.md                # OpenSpec apply command
│       └── archive.md              # OpenSpec archive command
├── settings.local.json             # Permission rules

docs/
├── RUNAI_QUICK_REFERENCE.md        # Command cheat sheet
└── MANUAL_RUNAI_EXECUTION.md       # Full execution guide

scripts/
├── monitor-runai-jobs.sh           # Job monitoring script
└── cleanup-runai.sh                # Job cleanup script
```

## Tested RunAI Commands

All commands use **RunAI CLI v2 syntax** and have been tested on the Salk HPI cluster.

### Submit ChIP-seq Analysis Job
```bash
runai workspace submit chipseq-analysis \
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
  --environment IP_SAMPLE=IP_col0 \
  --environment INPUT_SAMPLE=input_col0 \
  --command -- /scripts/chipseq/ChIProfiler_M00_Basic.sh
```

### Submit Single-Cell Analysis (Python/Scanpy)
```bash
runai workspace submit singlecell-scanpy \
  --project talmo-lab \
  --image ghcr.io/salk-harnessing-plants-initiative/chipseq-singlecell-pipeline:latest \
  --cpu-core-request 32 \
  --cpu-memory-request 128G \
  --gpu-devices-request 1 \
  --host-path path=/hpi/path/to/data,mount=/data,mount-propagation=HostToContainer \
  --host-path path=/hpi/path/to/outputs,mount=/outputs,mount-propagation=HostToContainer,readwrite \
  --command -- python /workspace/analysis.py
```

### Submit Single-Cell Analysis (R/Seurat)
```bash
runai workspace submit singlecell-seurat \
  --project talmo-lab \
  --image ghcr.io/salk-harnessing-plants-initiative/chipseq-singlecell-pipeline:latest \
  --cpu-core-request 32 \
  --cpu-memory-request 128G \
  --host-path path=/hpi/path/to/data,mount=/data,mount-propagation=HostToContainer \
  --host-path path=/hpi/path/to/outputs,mount=/outputs,mount-propagation=HostToContainer,readwrite \
  --command -- Rscript /workspace/analysis.R
```

### Submit VS Code Tunnel Workspace
```bash
runai workspace submit vscode-tunnel \
  --project talmo-lab \
  --image ghcr.io/salk-harnessing-plants-initiative/chipseq-singlecell-pipeline:latest \
  --cpu-core-request 8 \
  --cpu-memory-request 32G \
  --host-path path=/hpi/path/to/data,mount=/data,mount-propagation=HostToContainer \
  --host-path path=/hpi/path/to/outputs,mount=/outputs,mount-propagation=HostToContainer,readwrite \
  --host-path path=/hpi/path/to/references,mount=/references,mount-propagation=HostToContainer \
  --command -- tunnel
```

### List All Workspaces
```bash
runai workspace list
runai workspace list | grep chipseq
runai workspace list | grep singlecell
```

### Describe Workspace
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

## RunAI CLI v2 Migration Reference

| Old Syntax (v1) | New Syntax (v2) |
|-----------------|-----------------|
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

## Safety Guardrails

### Job Cleanup Safety
- Always preview with `--dry-run` before deleting
- Never delete running jobs without explicit `--force-running`
- Require confirmation for batch deletions
- Age-based filtering (only delete jobs older than threshold)
- Detailed logging of all deletions

### Command Permissions
The `.claude/settings.local.json` will define allowed bash patterns:
- `Bash(runai workspace *)` - All workspace commands
- `Bash(docker *)` - Docker operations
- `Bash(./scripts/*.sh)` - Script execution

## Testing Plan

1. **RunAI Command Validation**
   - [ ] Test `runai workspace submit` with ChIP-seq job
   - [ ] Test `runai workspace submit` with single-cell job
   - [ ] Test `runai workspace submit` with tunnel command
   - [ ] Test `runai workspace list` filtering
   - [ ] Test `runai workspace logs` with `--follow`
   - [ ] Test `runai workspace delete`

2. **Slash Command Validation**
   - [ ] Test `/submit-chipseq` generates correct command
   - [ ] Test `/submit-singlecell` generates correct command
   - [ ] Test `/submit-tunnel` generates correct command
   - [ ] Test `/monitor-jobs` displays status
   - [ ] Test `/cleanup-jobs` respects safety guardrails

3. **Documentation Validation**
   - [ ] All code blocks in docs are syntactically correct
   - [ ] All paths and image names are accurate
   - [ ] Examples match actual container capabilities

## Rollback Plan

If issues arise:
1. Delete `.claude/skills/runai.md`
2. Delete `.claude/commands/*.md` (except openspec/)
3. Delete `docs/RUNAI_*.md`
4. Delete `scripts/*-runai*.sh`
5. Restore any modified files from git

## Dependencies

- RunAI CLI v2 installed on user machines
- Access to `talmo-lab` project on Salk HPI cluster
- Docker image pushed to GHCR
- Host paths configured on cluster nodes

## Related Changes

This proposal is modeled after the successful RunAI integration in:
- `gapit3-gwas-pipeline` repository

## Approval Checklist

- [x] Technical approach reviewed
- [x] Command syntax verified against RunAI v2 docs
- [x] Safety guardrails adequate
- [x] Documentation complete
- [x] Testing plan sufficient
