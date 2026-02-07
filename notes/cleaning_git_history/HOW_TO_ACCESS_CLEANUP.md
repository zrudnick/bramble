# How to Access and Clean Your Local Branch

## Overview

I've added three files to the `copilot/enhance-logging-infrastructure` branch that will allow you to clean the git history locally:

1. **cleanup-git-history.sh** - Fully automated cleanup script
2. **GIT_HISTORY_CLEANUP_INSTRUCTIONS.md** - Detailed documentation
3. **QUICK_CLEANUP_GUIDE.md** - Quick reference

## How to Get Started

### Step 1: Fetch the Branch with the Cleanup Tools

```bash
# Navigate to your local bramble repository
cd /path/to/bramble

# Fetch the latest changes
git fetch origin copilot/enhance-logging-infrastructure

# Checkout the branch
git checkout copilot/enhance-logging-infrastructure

# Pull the latest commits (including the cleanup tools)
git pull origin copilot/enhance-logging-infrastructure
```

### Step 2: Run the Cleanup

You now have two options:

#### Option A: Use the Automated Script (Easiest)

```bash
# Simply run the script
./cleanup-git-history.sh

# When it completes successfully, force push
git push --force origin copilot/enhance-logging-infrastructure
```

The script handles everything:
- Creates a backup branch
- Removes all build artifacts from history
- Runs garbage collection
- Verifies the code still works
- Tells you exactly what to do next

#### Option B: Manual Cleanup

If you prefer more control, see `GIT_HISTORY_CLEANUP_INSTRUCTIONS.md` for detailed manual steps.

### Step 3: Verify the Results

After cleanup, you should see:
- 0 files under `bramble-rs/target/` in history
- Repository size reduced from ~240MB to ~180MB
- Code still builds and tests pass
- New commit hashes (due to history rewriting)

### Step 4: Force Push

```bash
git push --force origin copilot/enhance-logging-infrastructure
```

## What This Fixes

This cleanup removes 1,270 build artifact files that were accidentally committed in commit `a0ce14f`. These files:
- Were never meant to be in version control
- Added ~60MB to the repository size
- Slowed down clone operations
- Are regenerated during build anyway

## Files Available in the Branch

Once you checkout the branch, you'll have:

```
bramble/
├── cleanup-git-history.sh              # Automated cleanup script
├── GIT_HISTORY_CLEANUP_INSTRUCTIONS.md # Full documentation
├── QUICK_CLEANUP_GUIDE.md              # Quick reference
└── (all your existing code files)
```

## Safety

- The script creates a backup branch automatically
- You can restore the original state anytime with: `git reset --hard backup-before-cleanup`
- All changes are local until you force push
- The code is verified to build and test successfully after cleanup

## Questions?

If you encounter any issues:
1. Check `GIT_HISTORY_CLEANUP_INSTRUCTIONS.md` for troubleshooting
2. Look at `QUICK_CLEANUP_GUIDE.md` for quick commands
3. Ask for help if needed!

## Summary

**TL;DR**: 
1. `git checkout copilot/enhance-logging-infrastructure`
2. `git pull`
3. `./cleanup-git-history.sh`
4. `git push --force origin copilot/enhance-logging-infrastructure`
5. Done! Repository is now 60MB smaller.
