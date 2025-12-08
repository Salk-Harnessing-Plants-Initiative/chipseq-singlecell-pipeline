# OpenSpec

OpenSpec is a structured workflow for managing code changes through proposals, approvals, and implementation.

## Directory Structure

```
.openspec/
├── proposals/    # New change proposals awaiting review
├── approved/     # Approved proposals ready for implementation
├── archived/     # Completed and deployed changes
└── README.md     # This file
```

## Workflow

1. **Proposal** (`/openspec:proposal`) - Create a new change proposal in `proposals/`
2. **Review** - Team reviews and approves the proposal (move to `approved/`)
3. **Apply** (`/openspec:apply`) - Implement the approved change
4. **Archive** (`/openspec:archive`) - Archive completed changes to `archived/`

## Proposal Format

Each proposal should be a markdown file with:
- Clear title and description
- Motivation/problem statement
- Proposed solution
- Files affected
- Testing plan
- Rollback plan (if applicable)

## Usage

```bash
# Create a new proposal
/openspec:proposal

# Implement an approved proposal
/openspec:apply

# Archive a deployed change
/openspec:archive
```
