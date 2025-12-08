---
name: OpenSpec: Archive
description: Archive a deployed OpenSpec change and update specs.
category: OpenSpec
tags: [openspec, change]
---
<!-- OPENSPEC:START -->
**Guardrails**
- Only archive changes that have been fully implemented and deployed.
- Verify all tasks are complete before archiving.
- Preserve the proposal for historical reference.

**Steps**
1. Verify the change is complete:
   - All tasks marked done in tasks file
   - Proposal status is "Implemented"
   - Changes have been committed/merged
2. Move files to archived:
   ```bash
   mv .openspec/approved/<change-id>.md .openspec/archived/
   mv .openspec/approved/<change-id>-tasks.md .openspec/archived/  # if exists
   ```
3. Update the proposal with:
   - Final status: "Archived"
   - Completion date
   - Any post-implementation notes
4. Commit the archive operation.

**Archive Naming Convention**
Files in `.openspec/archived/` retain their original names:
- `<change-id>.md` - Main proposal
- `<change-id>-tasks.md` - Tasks list (if applicable)

**Verification Checklist**
- [ ] All implementation tasks complete
- [ ] Changes deployed/merged
- [ ] Proposal moved to archived folder
- [ ] Status updated to "Archived"
- [ ] Completion date recorded

**Reference**
- View archived changes: `ls .openspec/archived/`
- Review archived proposal: `cat .openspec/archived/<change-id>.md`
<!-- OPENSPEC:END -->
