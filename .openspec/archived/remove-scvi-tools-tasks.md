# Tasks: Remove scvi-tools Dependency

**Change ID**: `remove-scvi-tools`

## Task List

### Phase 1: Remove Dependencies

- [ ] **1.1** Update pyproject.toml
  - Remove `"scvi-tools"` from dependencies
  - Remove `[tool.uv]` PyTorch index configuration
  - Remove `[tool.uv.sources]` torch configuration

### Phase 2: Fix Container Tool Access

- [ ] **2.1** Update entrypoint.sh for uv best practices
  - Use `uv run` instead of manual venv activation
  - Update help text to show `uv run` examples
  - Remove scvi-tools from available tools list

- [ ] **2.2** Fix VS Code CLI download in Dockerfile
  - Update download URL (current returns 404)
  - Use stable download method

- [ ] **2.3** Update Dockerfile comment
  - Remove PyTorch reference from uv sync comment

### Phase 3: Documentation

- [ ] **3.1** Update README.md
  - Remove scvi-tools from Single-Cell Analysis (Python) section
  - Add note about optional scvi-tools installation

- [ ] **3.2** Update any other references
  - Check CLAUDE.md for scvi-tools mentions
  - Check docs/ for scvi-tools references

### Phase 4: Validation

- [ ] **4.1** Verify Docker build succeeds
  - Push changes to trigger GitHub Actions
  - Monitor build completion
  - Verify image is pushed to GHCR

- [ ] **4.2** Test image on RunAI
  - Pull image from GHCR
  - Verify Python environment works
  - Verify scanpy, anndata, bbknn, harmonypy work

### Phase 5: Commit and Cleanup

- [ ] **5.1** Commit all changes
  - Stage modified files
  - Create descriptive commit message

- [ ] **5.2** Update proposal status
  - Mark as Implemented when build succeeds

## Dependencies

```
Phase 1 → Phase 2 → Phase 3 → Phase 4 → Phase 5
```

All phases are sequential - each depends on the previous.

## Estimated Complexity

| Phase | Complexity | Files |
|-------|------------|-------|
| 1. Remove Dependencies | Low | 1 |
| 2. Fix Docker Build | Low | 1 |
| 3. Documentation | Low | 1-2 |
| 4. Validation | Medium | 0 |
| 5. Commit | Low | 0 |

**Total modified files**: 3-4
