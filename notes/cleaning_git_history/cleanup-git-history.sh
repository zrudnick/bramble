#!/bin/bash
#
# Git History Cleanup Script for bramble repository
# This script removes all build artifacts from bramble-rs/target/ directory
# that were accidentally committed in commit a0ce14ff9d1dcb924d4dd451e5abbcdbe63b187f
#

set -e  # Exit on error

echo "=================================================="
echo "Git History Cleanup Script for bramble repository"
echo "=================================================="
echo ""

# Check if we're in the correct repository
if [ ! -d "bramble-rs" ]; then
    echo "ERROR: This script must be run from the root of the bramble repository"
    exit 1
fi

# Check if the branch is correct
CURRENT_BRANCH=$(git rev-parse --abbrev-ref HEAD)
if [ "$CURRENT_BRANCH" != "copilot/enhance-logging-infrastructure" ]; then
    echo "WARNING: Current branch is '$CURRENT_BRANCH', not 'copilot/enhance-logging-infrastructure'"
    read -p "Do you want to continue anyway? (y/N) " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        exit 1
    fi
fi

echo "Step 1: Checking for build artifacts in current HEAD..."
ARTIFACT_COUNT=$(git ls-tree -r --name-only HEAD | grep "bramble-rs/target/" | wc -l)
echo "Found $ARTIFACT_COUNT files under bramble-rs/target/ in current HEAD"

if [ "$ARTIFACT_COUNT" -eq 0 ]; then
    echo ""
    echo "✓ No build artifacts found in current HEAD. History is already clean!"
    echo ""
    echo "Checking if history contains artifacts in older commits..."
    
    # Check the commit where artifacts were added
    if git rev-parse --verify a0ce14f >/dev/null 2>&1; then
        OLD_ARTIFACT_COUNT=$(git ls-tree -r --name-only a0ce14f | grep "bramble-rs/target/" | wc -l)
        echo "Commit a0ce14f has $OLD_ARTIFACT_COUNT files under bramble-rs/target/"
        
        if [ "$OLD_ARTIFACT_COUNT" -eq 0 ]; then
            echo "✓ History is already clean. No action needed!"
            exit 0
        fi
    fi
fi

echo ""
echo "Step 2: Creating backup of current state..."
BACKUP_BRANCH="backup-before-cleanup-$(date +%Y%m%d-%H%M%S)"
git branch "$BACKUP_BRANCH"
echo "✓ Created backup branch: $BACKUP_BRANCH"

echo ""
echo "Step 3: Running git filter-branch to remove bramble-rs/target/ from history..."
echo "This may take a few minutes..."
echo ""

# Find the parent of the commit where artifacts were added
if git rev-parse --verify f929c8b >/dev/null 2>&1; then
    START_COMMIT="f929c8b"
else
    echo "ERROR: Could not find expected commit f929c8b"
    echo "Please manually identify the commit before the artifacts were added"
    exit 1
fi

# Run filter-branch
git filter-branch --force --index-filter \
    'git rm -r --cached --ignore-unmatch bramble-rs/target/' \
    --prune-empty --tag-name-filter cat -- "${START_COMMIT}..HEAD" 2>&1 | tail -20

echo ""
echo "Step 4: Cleaning up old objects and running garbage collection..."
rm -rf .git/refs/original/
git reflog expire --expire=now --all
git gc --prune=now --aggressive

echo ""
echo "Step 5: Verifying cleanup..."
NEW_ARTIFACT_COUNT=$(git ls-tree -r --name-only HEAD | grep "bramble-rs/target/" | wc -l)
echo "Files under bramble-rs/target/ in current HEAD: $NEW_ARTIFACT_COUNT"

if [ "$NEW_ARTIFACT_COUNT" -eq 0 ]; then
    echo "✓ SUCCESS! All build artifacts removed from history"
else
    echo "✗ WARNING: Still found $NEW_ARTIFACT_COUNT files under bramble-rs/target/"
    echo "Cleanup may not have been complete"
fi

echo ""
echo "Step 6: Checking repository size..."
du -sh .git/

echo ""
echo "Step 7: Verifying code still works..."
cd bramble-rs
if cargo build --quiet; then
    echo "✓ Code builds successfully"
else
    echo "✗ WARNING: Build failed. You may need to investigate."
fi

if cargo test --quiet; then
    echo "✓ All tests pass"
else
    echo "✗ WARNING: Some tests failed. You may need to investigate."
fi
cd ..

echo ""
echo "=================================================="
echo "CLEANUP COMPLETE"
echo "=================================================="
echo ""
echo "Your branch now has a cleaned history with all build artifacts removed."
echo "The commit hashes have changed due to history rewriting."
echo ""
echo "NEXT STEPS:"
echo "1. Review the changes with: git log --oneline | head -10"
echo "2. If satisfied, force push to update the remote branch:"
echo "   git push --force origin $CURRENT_BRANCH"
echo ""
echo "If something went wrong, you can restore the original state with:"
echo "   git reset --hard $BACKUP_BRANCH"
echo ""
echo "Backup branch created: $BACKUP_BRANCH"
echo "=================================================="
