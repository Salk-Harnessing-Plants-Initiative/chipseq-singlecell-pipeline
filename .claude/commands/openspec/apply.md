---
name: OpenSpec: Apply
description: Implement an approved OpenSpec change and keep tasks in sync.
category: OpenSpec
tags: [openspec, change]
---
<!-- OPENSPEC:START -->
**Guardrails**
- Only implement changes that have been approved (moved to `.openspec/approved/`).
- Follow the tasks list if provided; update task status as you complete each item.
- Keep changes focused on what's specified in the proposal.
- Test each change before marking tasks complete.

**Steps**
1. Verify the change exists in `.openspec/approved/<change-id>.md`
2. Read the proposal and tasks file (if exists) to understand scope.
3. Implement changes according to the proposal:
   - Follow the tasks list in order if provided
   - Mark each task complete as you finish it
   - Commit logical units of work
4. Update the proposal status to "Implementing" while in progress.
5. When complete, update status to "Implemented" and note any deviations.
6. Run tests specified in the testing plan.
7. Move proposal to `.openspec/archived/` when deployed.

**Task Status Updates**
When a tasks file exists (`.openspec/approved/<change-id>-tasks.md`):
- Mark tasks `[x]` as you complete them
- Add notes for any blockers or changes
- Keep the file updated throughout implementation

**Completion Checklist**
- [ ] All tasks in tasks.md marked complete
- [ ] All files in "Files Affected" have been updated
- [ ] Testing plan executed successfully
- [ ] Proposal status updated to "Implemented"
- [ ] Changes committed with clear messages

**Reference**
- View approved changes: `ls .openspec/approved/`
- Check task status: `cat .openspec/approved/<change-id>-tasks.md`
<!-- OPENSPEC:END -->
