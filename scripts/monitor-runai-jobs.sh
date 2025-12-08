#!/bin/bash
set -euo pipefail

# =============================================================================
# RunAI Job Monitor for ChIP-seq & Single-Cell Pipeline
# =============================================================================
# Provides real-time monitoring dashboard for RunAI workspaces
#
# Usage:
#   ./scripts/monitor-runai-jobs.sh [OPTIONS]
#
# Options:
#   --prefix <name>    Filter jobs by prefix (default: all)
#   --project <name>   RunAI project (default: talmo-lab)
#   --watch            Auto-refresh every 30 seconds
#   --once             Single status check (no refresh)
#   --export <file>    Export status to JSON file
#   --help             Show this help message
# =============================================================================

# Default configuration
PROJECT="${PROJECT:-talmo-lab}"
PREFIX="${PREFIX:-}"
WATCH_MODE=false
EXPORT_FILE=""
REFRESH_INTERVAL=30

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --prefix)
            PREFIX="$2"
            shift 2
            ;;
        --project)
            PROJECT="$2"
            shift 2
            ;;
        --watch)
            WATCH_MODE=true
            shift
            ;;
        --once)
            WATCH_MODE=false
            shift
            ;;
        --export)
            EXPORT_FILE="$2"
            shift 2
            ;;
        --help|-h)
            echo "RunAI Job Monitor for ChIP-seq & Single-Cell Pipeline"
            echo ""
            echo "Usage: $0 [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  --prefix <name>    Filter jobs by prefix"
            echo "  --project <name>   RunAI project (default: talmo-lab)"
            echo "  --watch            Auto-refresh every 30 seconds"
            echo "  --once             Single status check (no refresh)"
            echo "  --export <file>    Export status to JSON file"
            echo "  --help             Show this help message"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Function to get workspace list
get_workspaces() {
    if [[ -n "$PREFIX" ]]; then
        runai workspace list 2>/dev/null | grep "$PREFIX" || true
    else
        runai workspace list 2>/dev/null | tail -n +2 || true
    fi
}

# Function to count jobs by status
count_by_status() {
    local status="$1"
    local workspaces="$2"
    echo "$workspaces" | grep -c "$status" || echo "0"
}

# Function to display dashboard
display_dashboard() {
    local workspaces
    workspaces=$(get_workspaces)

    # Count by status
    local running succeeded failed pending
    running=$(count_by_status "Running" "$workspaces")
    succeeded=$(count_by_status "Succeeded" "$workspaces")
    failed=$(count_by_status "Failed" "$workspaces")
    pending=$(count_by_status "Pending" "$workspaces")
    local total=$((running + succeeded + failed + pending))

    # Clear screen if in watch mode
    if $WATCH_MODE; then
        clear
    fi

    # Print header
    echo -e "${CYAN}═══════════════════════════════════════════════════════════════${NC}"
    echo -e "${CYAN}ChIP-seq & Single-Cell Pipeline - Job Monitor${NC}"
    echo -e "${CYAN}Updated: $(date '+%Y-%m-%d %H:%M:%S')${NC}"
    if [[ -n "$PREFIX" ]]; then
        echo -e "${CYAN}Filter: $PREFIX${NC}"
    fi
    echo -e "${CYAN}═══════════════════════════════════════════════════════════════${NC}"
    echo ""

    # Print summary
    echo -e "${BLUE}SUMMARY:${NC}"
    echo -e "  ${GREEN}Running:${NC}    $running jobs"
    echo -e "  ${GREEN}Succeeded:${NC}  $succeeded jobs"
    echo -e "  ${RED}Failed:${NC}     $failed jobs"
    echo -e "  ${YELLOW}Pending:${NC}    $pending jobs"
    echo -e "  Total:      $total jobs"
    echo ""

    # Show running jobs
    if [[ $running -gt 0 ]]; then
        echo -e "${BLUE}RUNNING JOBS:${NC}"
        echo "$workspaces" | grep "Running" | while read -r line; do
            local name status
            name=$(echo "$line" | awk '{print $1}')
            status=$(echo "$line" | awk '{print $2}')
            echo -e "  ${GREEN}●${NC} $name"
        done
        echo ""
    fi

    # Show pending jobs
    if [[ $pending -gt 0 ]]; then
        echo -e "${BLUE}PENDING JOBS:${NC}"
        echo "$workspaces" | grep "Pending" | while read -r line; do
            local name
            name=$(echo "$line" | awk '{print $1}')
            echo -e "  ${YELLOW}○${NC} $name"
        done
        echo ""
    fi

    # Show failed jobs
    if [[ $failed -gt 0 ]]; then
        echo -e "${BLUE}FAILED JOBS:${NC}"
        echo "$workspaces" | grep "Failed" | while read -r line; do
            local name
            name=$(echo "$line" | awk '{print $1}')
            echo -e "  ${RED}✗${NC} $name"
        done
        echo ""
    fi

    # Show recent completed (last 5)
    local completed_list
    completed_list=$(echo "$workspaces" | grep "Succeeded" | head -5)
    if [[ -n "$completed_list" ]]; then
        echo -e "${BLUE}RECENTLY COMPLETED (last 5):${NC}"
        echo "$completed_list" | while read -r line; do
            local name
            name=$(echo "$line" | awk '{print $1}')
            echo -e "  ${GREEN}✓${NC} $name"
        done
        echo ""
    fi

    # Export to JSON if requested
    if [[ -n "$EXPORT_FILE" ]]; then
        cat > "$EXPORT_FILE" << EOF
{
  "timestamp": "$(date -Iseconds)",
  "project": "$PROJECT",
  "prefix": "$PREFIX",
  "summary": {
    "running": $running,
    "succeeded": $succeeded,
    "failed": $failed,
    "pending": $pending,
    "total": $total
  }
}
EOF
        echo -e "${CYAN}Status exported to: $EXPORT_FILE${NC}"
    fi

    # Show refresh message if in watch mode
    if $WATCH_MODE; then
        echo -e "${CYAN}Press Ctrl+C to exit. Refreshing every ${REFRESH_INTERVAL}s...${NC}"
    fi
}

# Main execution
if $WATCH_MODE; then
    while true; do
        display_dashboard
        sleep "$REFRESH_INTERVAL"
    done
else
    display_dashboard
fi
