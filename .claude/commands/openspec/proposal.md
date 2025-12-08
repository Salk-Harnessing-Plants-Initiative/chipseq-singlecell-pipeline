---
name: OpenSpec: Proposal
description: Scaffold a new OpenSpec change and validate strictly.
category: OpenSpec
tags: [openspec, change]
---
<!-- OPENSPEC:START -->
**Guardrails**
- Favor straightforward, minimal implementations first and add complexity only when it is requested or clearly required.
- Keep changes tightly scoped to the requested outcome.
- Refer to `.openspec/README.md` for OpenSpec conventions and directory structure.
- Identify any vague or ambiguous details and ask the necessary follow-up questions before editing files.

**Steps**
1. Review `.openspec/README.md` and inspect related code or docs to ground the proposal in current behaviour; note any gaps that require clarification.
2. Choose a unique verb-led `change-id` (e.g., `add-feature-x`, `fix-bug-y`, `update-config-z`).
3. Create proposal file at `.openspec/proposals/<change-id>.md` with:
   - Summary of the change
   - Motivation / problem statement
   - Proposed solution
   - Files affected
   - Testing plan
   - Rollback plan (if applicable)
4. Optionally create `.openspec/proposals/<change-id>-tasks.md` with ordered list of implementation tasks.
5. Validate the proposal is complete and actionable before sharing.

**Proposal Template**
```markdown
# OpenSpec Proposal: <Title>

**Change ID**: `<change-id>`
**Created**: <date>
**Status**: Proposal

## Summary
<1-2 sentence description>

## Motivation / Problem Statement
<Why is this change needed?>

## Proposed Solution
<What will be changed and how?>

## Files Affected
- `path/to/file1` - <description>
- `path/to/file2` - <description>

## Testing Plan
- [ ] Test case 1
- [ ] Test case 2

## Rollback Plan
<How to undo if needed>
```

**Reference**
- View existing proposals: `ls .openspec/proposals/`
- View approved changes: `ls .openspec/approved/`
- View archived changes: `ls .openspec/archived/`
<!-- OPENSPEC:END -->
