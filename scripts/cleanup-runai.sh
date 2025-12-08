#!/bin/bash
set -euo pipefail

# =============================================================================
# RunAI Job Cleanup for ChIP-seq & Single-Cell Pipeline
# =============================================================================
# Safely clean up completed, failed, or orphaned RunAI workspaces
#
# Usage:
#   ./scripts/cleanup-runai.sh [OPTIONS]
#
# Options:
#   --dry-run          Preview without deleting
#   --completed        Delete only succeeded jobs
#   --failed           Delete only failed jobs
#   --all              Delete both completed and failed
#   --prefix <name>    Only jobs matching prefix
#   --older-than <age> Only jobs older than threshold (e.g., 24h, 7d)
#   --yes              Skip confirmation prompt
#   --force-running    Include running jobs (DANGEROUS!)
#   --help             Show this help message
# =============================================================================

# Default configuration
PROJECT="${PROJECT:-talmo-lab}"
DRY_RUN=false
DELETE_COMPLETED=false
DELETE_FAILED=false
PREFIX=""
OLDER_THAN=""
SKIP_CONFIRM=false
FORCE_RUNNING=false

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        --completed)
            DELETE_COMPLETED=true
            shift
            ;;
        --failed)
            DELETE_FAILED=true
            shift
            ;;
        --all)
            DELETE_COMPLETED=true
            DELETE_FAILED=true
            shift
            ;;
        --prefix)
            PREFIX="$2"
            shift 2
            ;;
        --older-than)
            OLDER_THAN="$2"
            shift 2
            ;;
        --yes)
            SKIP_CONFIRM=true
            shift
            ;;
        --force-running)
            FORCE_RUNNING=true
            shift
            ;;
        --project)
            PROJECT="$2"
            shift 2
            ;;
        --help|-h)
            echo "RunAI Job Cleanup for ChIP-seq & Single-Cell Pipeline"
            echo ""
            echo "Usage: $0 [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  --dry-run          Preview without deleting"
            echo "  --completed        Delete only succeeded jobs"
            echo "  --failed           Delete only failed jobs"
            echo "  --all              Delete both completed and failed"
            echo "  --prefix <name>    Only jobs matching prefix"
            echo "  --older-than <age> Only jobs older than threshold (e.g., 24h, 7d)"
            echo "  --yes              Skip confirmation prompt"
            echo "  --force-running    Include running jobs (DANGEROUS!)"
            echo "  --project <name>   RunAI project (default: talmo-lab)"
            echo "  --help             Show this help message"
            echo ""
            echo "Examples:"
            echo "  $0 --completed --dry-run          # Preview completed job cleanup"
            echo "  $0 --failed --yes                 # Delete all failed jobs"
            echo "  $0 --all --prefix chipseq         # Delete all chipseq jobs"
            echo "  $0 --completed --older-than 24h   # Delete completed jobs older than 24h"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Interactive mode if no deletion type specified
if ! $DELETE_COMPLETED && ! $DELETE_FAILED; then
    echo -e "${BLUE}RunAI Job Cleanup - Interactive Mode${NC}"
    echo ""

    # Get counts
    workspaces=$(runai workspace list 2>/dev/null | tail -n +2 || true)
    if [[ -n "$PREFIX" ]]; then
        workspaces=$(echo "$workspaces" | grep "$PREFIX" || true)
    fi

    running=$(echo "$workspaces" | grep -c "Running" || echo "0")
    completed=$(echo "$workspaces" | grep -c "Succeeded" || echo "0")
    failed=$(echo "$workspaces" | grep -c "Failed" || echo "0")

    echo "Found:"
    echo -e "  ${GREEN}$completed completed jobs${NC}"
    echo -e "  ${RED}$failed failed jobs${NC}"
    echo -e "  ${YELLOW}$running running jobs (protected)${NC}"
    echo ""

    echo "What would you like to clean up?"
    echo "1) Completed jobs only"
    echo "2) Failed jobs only"
    echo "3) Both completed and failed"
    echo "4) Show detailed list first"
    echo "5) Cancel"
    echo ""
    read -rp "Select option: " choice

    case $choice in
        1) DELETE_COMPLETED=true ;;
        2) DELETE_FAILED=true ;;
        3) DELETE_COMPLETED=true; DELETE_FAILED=true ;;
        4)
            echo ""
            echo "=== ALL WORKSPACES ==="
            echo "$workspaces"
            echo ""
            exit 0
            ;;
        5|*)
            echo "Cancelled."
            exit 0
            ;;
    esac
fi

# Get workspaces to delete
workspaces=$(runai workspace list 2>/dev/null | tail -n +2 || true)

# Apply prefix filter
if [[ -n "$PREFIX" ]]; then
    workspaces=$(echo "$workspaces" | grep "$PREFIX" || true)
fi

# Build list of jobs to delete
to_delete=""

if $DELETE_COMPLETED; then
    completed_jobs=$(echo "$workspaces" | grep "Succeeded" || true)
    to_delete="$to_delete$completed_jobs"$'\n'
fi

if $DELETE_FAILED; then
    failed_jobs=$(echo "$workspaces" | grep "Failed" || true)
    to_delete="$to_delete$failed_jobs"$'\n'
fi

if $FORCE_RUNNING; then
    echo -e "${RED}WARNING: --force-running is enabled. Running jobs will be deleted!${NC}"
    running_jobs=$(echo "$workspaces" | grep "Running" || true)
    to_delete="$to_delete$running_jobs"$'\n'
fi

# Clean up empty lines
to_delete=$(echo "$to_delete" | grep -v '^$' || true)

# Check if there's anything to delete
if [[ -z "$to_delete" ]]; then
    echo -e "${GREEN}No jobs to delete matching the criteria.${NC}"
    exit 0
fi

# Count jobs
job_count=$(echo "$to_delete" | wc -l | tr -d ' ')

# Show what will be deleted
echo ""
echo -e "${BLUE}Jobs to delete ($job_count):${NC}"
echo "$to_delete" | while read -r line; do
    if [[ -n "$line" ]]; then
        name=$(echo "$line" | awk '{print $1}')
        status=$(echo "$line" | awk '{print $2}')
        case $status in
            Succeeded) echo -e "  ${GREEN}✓${NC} $name ($status)" ;;
            Failed)    echo -e "  ${RED}✗${NC} $name ($status)" ;;
            Running)   echo -e "  ${YELLOW}●${NC} $name ($status)" ;;
            *)         echo "  - $name ($status)" ;;
        esac
    fi
done
echo ""

# Dry run mode
if $DRY_RUN; then
    echo -e "${YELLOW}DRY RUN: No jobs were deleted.${NC}"
    echo "Run without --dry-run to actually delete these jobs."
    exit 0
fi

# Confirmation
if ! $SKIP_CONFIRM; then
    echo -e "${RED}WARNING: This will permanently delete $job_count job(s) and their logs.${NC}"
    read -rp "Type 'DELETE' to confirm: " confirm
    if [[ "$confirm" != "DELETE" ]]; then
        echo "Cancelled."
        exit 0
    fi
fi

# Delete jobs
echo ""
echo "Deleting jobs..."
deleted=0
errors=0

echo "$to_delete" | while read -r line; do
    if [[ -n "$line" ]]; then
        name=$(echo "$line" | awk '{print $1}')
        echo -n "  Deleting $name... "
        if runai workspace delete "$name" -p "$PROJECT" 2>/dev/null; then
            echo -e "${GREEN}OK${NC}"
            ((deleted++)) || true
        else
            echo -e "${RED}FAILED${NC}"
            ((errors++)) || true
        fi
    fi
done

echo ""
echo -e "${GREEN}Cleanup complete.${NC}"
echo "  Deleted: $deleted"
if [[ $errors -gt 0 ]]; then
    echo -e "  ${RED}Errors: $errors${NC}"
fi
