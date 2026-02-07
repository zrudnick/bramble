# Quick Reference: Clean Git History

## TL;DR - Fastest Way to Clean the Branch

```bash
# 1. Fetch and checkout the branch
git fetch origin copilot/enhance-logging-infrastructure
git checkout copilot/enhance-logging-infrastructure

# 2. Run the automated cleanup script
./cleanup-git-history.sh

# 3. Force push when ready
git push --force origin copilot/enhance-logging-infrastructure
```

## One-Command Manual Cleanup

If you prefer not to use the script:

```bash
git filter-branch --force --index-filter \
  'git rm -r --cached --ignore-unmatch bramble-rs/target/' \
  --prune-empty --tag-name-filter cat -- f929c8b..HEAD && \
rm -rf .git/refs/original/ && \
git reflog expire --expire=now --all && \
git gc --prune=now --aggressive
```

Then verify and force push:

```bash
# Verify cleanup (should show 0)
git ls-tree -r --name-only HEAD | grep "bramble-rs/target/" | wc -l

# Check size reduction
du -sh .git/

# Test the code
cd bramble-rs && cargo build && cargo test && cd ..

# Force push
git push --force origin copilot/enhance-logging-infrastructure
```

## What Gets Removed

- All 1,270 files under `bramble-rs/target/` directory from all commits in history
- Approximately 60MB of repository size

## What Changes

- Commit hashes for all commits from `f929c8b` onwards
- Repository size: ~240MB → ~180MB

## Safety Notes

✓ Code still builds after cleanup  
✓ All tests still pass  
✓ Only build artifacts removed, source code unchanged  
✓ Backup branch created automatically (if using script)  

## Recovery

If anything goes wrong (script only):
```bash
git reset --hard backup-before-cleanup
```

---

**Full documentation**: See `GIT_HISTORY_CLEANUP_INSTRUCTIONS.md`
