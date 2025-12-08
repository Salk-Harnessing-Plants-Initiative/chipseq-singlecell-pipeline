# Tasks: Add RunAI Skill and Claude Commands

**Change ID**: `add-runai-skill-and-commands`

## Task List

### Phase 1: Core Infrastructure

- [ ] **1.1** Create `.claude/` directory structure
  - Create `.claude/skills/` directory
  - Create `.claude/commands/` directory
  - Create `.claude/commands/openspec/` directory

- [ ] **1.2** Create Claude settings file
  - Create `.claude/settings.local.json` with permission rules
  - Allow RunAI workspace commands
  - Allow script execution

### Phase 2: RunAI Skill

- [ ] **2.1** Create RunAI skill document (`.claude/skills/runai.md`)
  - Project configuration (talmo-lab)
  - CLI v2 command reference
  - ChIP-seq job templates
  - Single-cell job templates (Python + R)
  - VS Code tunnel templates
  - Troubleshooting guide
  - Migration reference (v1 → v2)

### Phase 3: Slash Commands

- [ ] **3.1** Create job submission commands
  - `/submit-chipseq` - ChIP-seq analysis submission
  - `/submit-singlecell` - Single-cell analysis submission
  - `/submit-tunnel` - VS Code tunnel workspace

- [ ] **3.2** Create job management commands
  - `/monitor-jobs` - Job status monitoring
  - `/cleanup-jobs` - Safe job cleanup with guardrails

- [ ] **3.3** Create development commands
  - `/docker-build` - Build Docker image
  - `/docker-test` - Test Docker container
  - `/validate-bash` - ShellCheck validation

- [ ] **3.4** Create OpenSpec commands
  - `/openspec:proposal` - Scaffold new proposal
  - `/openspec:apply` - Implement approved change
  - `/openspec:archive` - Archive deployed change

### Phase 4: Helper Scripts

- [ ] **4.1** Create monitoring script (`scripts/monitor-runai-jobs.sh`)
  - Real-time dashboard with status counts
  - Filter by job prefix
  - Export to JSON option
  - Auto-refresh mode

- [ ] **4.2** Create cleanup script (`scripts/cleanup-runai.sh`)
  - Dry-run preview mode
  - Status-based filtering (completed/failed)
  - Age-based filtering
  - Running job protection
  - Confirmation prompts
  - Logging of deletions

### Phase 5: Documentation

- [ ] **5.1** Create quick reference (`docs/RUNAI_QUICK_REFERENCE.md`)
  - Project info section
  - Quick commands for all operations
  - Common workflows
  - Troubleshooting tips
  - CLI v2 migration table

- [ ] **5.2** Create full execution guide (`docs/MANUAL_RUNAI_EXECUTION.md`)
  - Prerequisites
  - Step-by-step workflows
  - ChIP-seq analysis guide
  - Single-cell analysis guide
  - VS Code tunnel guide
  - Monitoring and cleanup
  - Results verification

### Phase 6: Validation

- [ ] **6.1** Test RunAI commands on cluster
  - Submit test ChIP-seq job
  - Submit test single-cell job
  - Submit test tunnel workspace
  - Verify logs and status commands
  - Test cleanup commands

- [ ] **6.2** Test slash commands
  - Verify command expansion
  - Check parameter substitution
  - Validate safety guardrails

- [ ] **6.3** Documentation review
  - Verify all code blocks
  - Check paths and image names
  - Ensure consistency across docs

### Phase 7: Integration

- [ ] **7.1** Update CLAUDE.md
  - Add RunAI skill reference
  - Document new slash commands

- [ ] **7.2** Update README.md
  - Add RunAI quick start section
  - Link to detailed documentation

## Parallelizable Work

The following tasks can be done in parallel:
- Phase 2 (Skill) and Phase 3 (Commands) can proceed simultaneously
- Phase 4 (Scripts) and Phase 5 (Documentation) can proceed simultaneously
- All Phase 3 subtasks (3.1, 3.2, 3.3, 3.4) can be parallelized

## Dependencies

```
Phase 1 → Phase 2, 3, 4
Phase 2, 3, 4, 5 → Phase 6
Phase 6 → Phase 7
```

## Estimated Complexity

| Phase | Complexity | Files |
|-------|------------|-------|
| 1. Infrastructure | Low | 2 |
| 2. RunAI Skill | Medium | 1 |
| 3. Slash Commands | Medium | 11 |
| 4. Helper Scripts | Medium | 2 |
| 5. Documentation | Medium | 2 |
| 6. Validation | Low | 0 |
| 7. Integration | Low | 2 |

**Total new files**: ~20
