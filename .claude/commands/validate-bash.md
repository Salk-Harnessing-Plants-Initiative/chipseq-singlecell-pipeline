# Validate Bash Scripts

Validate shell scripts with syntax checking and ShellCheck static analysis.

## Quick Validation

```bash
# Validate all scripts
./scripts/validate-bash.sh

# Or manually with ShellCheck
shellcheck scripts/**/*.sh
```

## Validation Checks

### 1. Bash Syntax Check

```bash
# Check syntax without executing
bash -n scripts/entrypoint.sh
bash -n scripts/chipseq/ChIProfiler_M00_Basic.sh
bash -n scripts/chipseq/ChIProfiler_M01_Stat_PE.sh
bash -n scripts/chipseq/ChIProfiler_M03_CallPeaks.sh
```

### 2. ShellCheck Analysis

```bash
# Run ShellCheck on all scripts
shellcheck scripts/**/*.sh

# With specific exclusions (matching CI)
shellcheck \
  --exclude=SC1090,SC1091,SC2034,SC2046,SC2086,SC2155 \
  scripts/**/*.sh
```

### 3. Shebang Check

```bash
# Verify scripts have proper shebang
head -1 scripts/**/*.sh

# Should show:
# #!/bin/bash
# or
# #!/usr/bin/env bash
```

### 4. Safe Options Check

```bash
# Check for set -e or set -euo pipefail
grep -l "set -e" scripts/**/*.sh
```

## ShellCheck Exclusions

The CI workflow excludes these checks:

| Code | Reason |
|------|--------|
| SC1090 | Can't follow sourced files |
| SC1091 | Source file not found |
| SC2034 | Unused variables (may be external) |
| SC2046 | Word splitting (intentional for parallel) |
| SC2086 | Globbing (intentional) |
| SC2155 | Declare/assign separation |

## Validate Specific Script

```bash
# Validate single script
shellcheck scripts/chipseq/ChIProfiler_M00_Basic.sh

# With verbose output
shellcheck -x scripts/chipseq/ChIProfiler_M00_Basic.sh

# Show only errors (not warnings)
shellcheck -S error scripts/chipseq/ChIProfiler_M00_Basic.sh
```

## Common Issues and Fixes

### Unquoted Variables (SC2086)

```bash
# Bad
cp $file $dest

# Good
cp "$file" "$dest"
```

### Unused Variables (SC2034)

```bash
# If intentionally unused, prefix with _
_unused_var="value"

# Or export if used externally
export MY_VAR="value"
```

### Word Splitting (SC2046)

```bash
# If intentional, disable for line
# shellcheck disable=SC2046
files=$(ls *.txt)
```

### Missing Shebang

```bash
# Add at top of script
#!/bin/bash
set -euo pipefail
```

## CI Validation

The GitHub Actions workflow validates on every push/PR:

```yaml
# .github/workflows/validate-bash-scripts.yml
- name: Validate with ShellCheck
  run: |
    shellcheck \
      --exclude=SC1090,SC1091,SC2034,SC2046,SC2086,SC2155 \
      scripts/**/*.sh
```

## Pre-commit Hook

Add to `.git/hooks/pre-commit`:

```bash
#!/bin/bash
# Validate bash scripts before commit

SCRIPTS=$(git diff --cached --name-only --diff-filter=ACM | grep '\.sh$')
if [ -n "$SCRIPTS" ]; then
    echo "Validating shell scripts..."
    shellcheck $SCRIPTS
    if [ $? -ne 0 ]; then
        echo "ShellCheck found issues. Please fix before committing."
        exit 1
    fi
fi
```

## Validate All Scripts in Repo

```bash
# Find and validate all .sh files
find . -name "*.sh" -type f | xargs shellcheck

# Excluding certain directories
find . -name "*.sh" -type f -not -path "./.git/*" | xargs shellcheck
```

## Expected Output

Successful validation:
```
Checking scripts/entrypoint.sh... OK
Checking scripts/chipseq/ChIProfiler_M00_Basic.sh... OK
Checking scripts/chipseq/ChIProfiler_M01_Stat_PE.sh... OK
Checking scripts/chipseq/ChIProfiler_M03_CallPeaks.sh... OK

All scripts validated successfully!
```

## Troubleshooting

### ShellCheck not installed

```bash
# macOS
brew install shellcheck

# Ubuntu/Debian
apt-get install shellcheck

# Or use Docker
docker run --rm -v "$PWD:/mnt" koalaman/shellcheck:stable scripts/**/*.sh
```

### False positives

Add inline disable comments:
```bash
# shellcheck disable=SC2086
command $unquoted_var
```

Or file-level:
```bash
#!/bin/bash
# shellcheck disable=SC2086,SC2046
```

## Related Commands

- `/docker-build` - Build after validation
- `/docker-test` - Test container
