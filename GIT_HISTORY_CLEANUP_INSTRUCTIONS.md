# Git History Cleanup Guide

## Problem

Commit `a0ce14ff9d1dcb924d4dd451e5abbcdbe63b187f` accidentally committed 1,270 build artifact files from the `bramble-rs/target/` directory, bloating the repository by approximately 60MB.

## Solution Options

You have two options to clean the git history:

### Option 1: Automated Script (Recommended)

Run the provided cleanup script:

```bash
./cleanup-git-history.sh
```

The script will:
- ✓ Create a backup branch automatically
- ✓ Remove all files under `bramble-rs/target/` from git history
- ✓ Run garbage collection to reduce repository size
- ✓ Verify that code still builds and tests pass
- ✓ Provide clear next steps

After the script completes successfully, force push:
```bash
git push --force origin copilot/enhance-logging-infrastructure
```

### Option 2: Manual Cleanup

If you prefer to do it manually or the script doesn't work for your setup:

#### Step 1: Fetch the branch

```bash
# Fetch the branch
git fetch origin copilot/enhance-logging-infrastructure
git checkout copilot/enhance-logging-infrastructure
```

#### Step 2: Create a backup

```bash
# Create a backup branch in case something goes wrong
git branch backup-before-cleanup
```

#### Step 3: Remove build artifacts from history

```bash
# Run filter-branch to remove all files under bramble-rs/target/
git filter-branch --force --index-filter \
    'git rm -r --cached --ignore-unmatch bramble-rs/target/' \
    --prune-empty --tag-name-filter cat -- f929c8b..HEAD
```

This will rewrite history starting from commit `f929c8b` (the commit before the artifacts were added).

#### Step 4: Clean up and reduce repository size

```bash
# Remove the backup refs created by filter-branch
rm -rf .git/refs/original/

# Expire all reflog entries and garbage collect
git reflog expire --expire=now --all
git gc --prune=now --aggressive
```

#### Step 5: Verify the cleanup

```bash
# Check that no files remain under bramble-rs/target/
git ls-tree -r --name-only HEAD | grep "bramble-rs/target/" | wc -l
# Should output: 0

# Check repository size
du -sh .git/
# Should be around 180-190MB (down from 240MB)

# Verify code still works
cd bramble-rs
cargo build
cargo test
cd ..
```

#### Step 6: Force push the cleaned branch

```bash
git push --force origin copilot/enhance-logging-infrastructure
```

## Expected Results

After cleanup:

- **Files removed**: 1,270 build artifact files from history
- **Repository size**: Reduced from ~240MB to ~180MB (60MB savings)
- **Commit hashes**: All commits from `f929c8b` onwards will have new hashes
- **Code**: Still builds and all tests pass

## Commit Hash Changes

Due to history rewriting, commit hashes will change:

| Description | Old Hash | New Hash |
|-------------|----------|----------|
| Add logging and progress bar infrastructure | a0ce14f | 5f6ab12* |
| Update gitignore to exclude Rust build artifacts | ccce3ea | db35b53* |
| Address code review feedback: use constant and expect | 48392e6 | 1520735* |
| Improve progress bar refresh rate and update frequency | dca6502 | 4887445* |

*Note: New hashes will be generated when you run the cleanup. The exact hashes may vary.

## Troubleshooting

### If you get authentication errors

Make sure you have push access to the repository and are authenticated:
```bash
git config --list | grep remote.origin.url
# Should show: https://github.com/zrudnick/bramble
```

### If you get "refusing to lose history" errors

Use `--force` flag with git push (this is expected and safe since we intentionally rewrote history):
```bash
git push --force origin copilot/enhance-logging-infrastructure
```

### If something goes wrong

Restore from your backup:
```bash
git reset --hard backup-before-cleanup
```

### If you want to verify what changed

Compare old vs new history:
```bash
# Show commits in your cleaned branch
git log --oneline | head -10

# Show what files changed in each commit
git log --stat --oneline | head -50
```

## Why This Cleanup Is Necessary

1. **Repository size**: The build artifacts add ~60MB of unnecessary data to every clone
2. **Clone performance**: Larger repositories take longer to clone and consume more disk space
3. **Good practice**: Build artifacts should never be committed to version control
4. **CI/CD efficiency**: Smaller repositories mean faster CI/CD pipeline execution

## After Force Push

Once you've force-pushed the cleaned branch:

1. Other collaborators should update their local branches:
   ```bash
   git fetch origin
   git reset --hard origin/copilot/enhance-logging-infrastructure
   ```

2. The repository will be smaller for all future clones

3. The PR will be updated automatically with the new commits

## Questions?

If you encounter any issues or have questions about this cleanup process, please ask!
