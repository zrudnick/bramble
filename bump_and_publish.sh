#!/usr/bin/env bash
# Bump and (optionally) publish the two bramble Rust crates in lockstep:
#
#   bramble-rs   — the I/O-free projection library          (published first)
#   bramble-cli  — the BAM-in/BAM-out CLI (exe: `bramble-rs`) (published second)
#
# Both crates share one version and one git tag (v<version>). bramble-cli depends
# on bramble-rs by `{ path = "../bramble-rs", version = "<version>" }`, so its
# version requirement is bumped too, and the library must reach crates.io before
# the binary is published.
#
# Based on piscem's bump_and_publish.sh, extended for the two-crate split.
set -euo pipefail

die() {
    echo "error: $*" >&2
    exit 1
}

usage() {
    cat <<'EOF'
Usage:
  ./bump_and_publish.sh <version> [--publish] [--dry-run]

Options:
  --publish  After bumping/committing/tagging/pushing, publish BOTH crates to
             crates.io in order: bramble-rs (library) then bramble-cli (binary).
  --dry-run  Show what would be done without modifying tracked files, creating
             commits or tags, pushing, or publishing.
  -h, --help Show this help message.
EOF
}

print_cmd() {
    printf '+'
    printf ' %q' "$@"
    printf '\n'
}

run() {
    print_cmd "$@"
    if [[ "$DRY_RUN" == true ]]; then
        return 0
    fi
    "$@"
}

VERSION=""
PUBLISH=false
DRY_RUN=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        --publish) PUBLISH=true ;;
        --dry-run) DRY_RUN=true ;;
        -h|--help) usage; exit 0 ;;
        -*) die "unknown option: $1" ;;
        *)
            [[ -n "$VERSION" ]] && die "version specified more than once"
            VERSION="$1"
            ;;
    esac
    shift
done

[[ -n "$VERSION" ]] || { usage; exit 1; }

if ! [[ "$VERSION" =~ ^[0-9]+\.[0-9]+\.[0-9]+([+-][0-9A-Za-z.-]+)*$ ]]; then
    die "version must look like X.Y.Z, optionally with prerelease/build suffixes"
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

LIB_DIR="bramble-rs"
CLI_DIR="bramble-cli"
LIB_CARGO="$LIB_DIR/Cargo.toml"
CLI_CARGO="$CLI_DIR/Cargo.toml"
LIB_LOCK="$LIB_DIR/Cargo.lock"
CLI_LOCK="$CLI_DIR/Cargo.lock"
TAG="v${VERSION}"

for f in "$LIB_CARGO" "$CLI_CARGO" "$LIB_LOCK" "$CLI_LOCK"; do
    [[ -f "$f" ]] || die "not found: $f"
done

# --- preconditions -----------------------------------------------------------

CURRENT_VERSION="$(sed -n 's/^version = "\(.*\)"/\1/p' "$LIB_CARGO" | head -1)"
[[ -n "$CURRENT_VERSION" ]] || die "could not determine current version from $LIB_CARGO"

CLI_CURRENT="$(sed -n 's/^version = "\(.*\)"/\1/p' "$CLI_CARGO" | head -1)"
[[ "$CLI_CURRENT" == "$CURRENT_VERSION" ]] || \
    die "crate versions are out of sync: $LIB_DIR=$CURRENT_VERSION, $CLI_DIR=$CLI_CURRENT"

[[ "$CURRENT_VERSION" == "$VERSION" ]] && die "crates are already at version $VERSION"

git rev-parse "$TAG" >/dev/null 2>&1 && die "tag $TAG already exists"
[[ -z "$(git status --porcelain)" ]] || \
    die "working tree is not clean; commit or stash existing changes first"
git remote get-url origin >/dev/null 2>&1 || die "git remote 'origin' is not configured"

# Publish blocker: cargo publish rejects git/path dependencies. The library
# currently pins ksw2rs to a git fork (crates.io ksw2rs 0.1.3 is broken). Refuse
# to attempt publishing until that is a plain crates.io version requirement.
GIT_DEPS="$(sed -n '/^\[dependencies\]/,/^\[/p' "$LIB_CARGO" | grep -nE '\bgit *=' || true)"
if [[ "$PUBLISH" == true && -n "$GIT_DEPS" ]]; then
    echo "$GIT_DEPS" >&2
    die "$LIB_CARGO has git dependencies; crates.io will reject \`cargo publish\`.
       Switch ksw2rs (and any other git dep) to a published crates.io version first.
       Re-run without --publish to bump/commit/tag/push only."
fi

echo "Current version : $CURRENT_VERSION"
echo "New version     : $VERSION"
echo "Tag             : $TAG"
echo "Crates          : $LIB_DIR, $CLI_DIR"
echo "Publish         : $PUBLISH"
echo "Dry-run         : $DRY_RUN"
echo

# --- helpers -----------------------------------------------------------------

# Bump the top-of-file `version = "..."` of a Cargo.toml (the [package] version).
bump_manifest_version() {
    local cargo="$1"
    sed -i.bak "1,/^version = /s/^version = \".*\"/version = \"${VERSION}\"/" "$cargo"
    rm -f "${cargo}.bak"
}

# Bump the version requirement of a named [[package]] in a Cargo.lock, scoped to
# that package's stanza (name line .. following blank line) so other packages of
# the same version are untouched.
bump_lock_pkg() {
    local lock="$1" pkg="$2"
    sed -i.bak "/^name = \"${pkg}\"$/,/^$/s/^version = \".*\"/version = \"${VERSION}\"/" "$lock"
    rm -f "${lock}.bak"
}

# Bump bramble-cli's dependency requirement on bramble-rs:
#   bramble-rs = { path = "../bramble-rs", version = "X" }  ->  version = "<VERSION>"
bump_cli_dep_on_lib() {
    sed -i.bak -E "s#(bramble-rs *= *\{[^}]*version *= *\")[^\"]*(\")#\1${VERSION}\2#" "$CLI_CARGO"
    rm -f "${CLI_CARGO}.bak"
}

backup() { cp "$1" "$1.release.bak"; }
restore_all() {
    for f in "$LIB_CARGO" "$CLI_CARGO" "$LIB_LOCK" "$CLI_LOCK"; do
        [[ -f "$f.release.bak" ]] && cp "$f.release.bak" "$f"
    done
}
clear_backups() {
    for f in "$LIB_CARGO" "$CLI_CARGO" "$LIB_LOCK" "$CLI_LOCK"; do
        rm -f "$f.release.bak"
    done
}

COMMIT_CREATED=false
cleanup() {
    local status=$?
    if [[ "$status" -ne 0 && "$DRY_RUN" == false && "$COMMIT_CREATED" == false ]]; then
        restore_all
        echo "restored manifests/lockfiles after failure" >&2
    fi
    clear_backups
    return "$status"
}
trap cleanup EXIT

# package (no compile) sanity check for one crate dir
package_check() {
    local dir="$1"
    local tmp
    tmp="$(mktemp -d "${TMPDIR:-/tmp}/bramble-release-check.XXXXXX")"
    ( cd "$dir" && CARGO_TARGET_DIR="$tmp" cargo package --allow-dirty --no-verify >/dev/null )
    rm -rf "$tmp"
}

# --- preflight (pre-bump) ----------------------------------------------------
#
# Only the library is package-checked here. bramble-cli depends on bramble-rs by
# `{ path, version }`; on publish the path is stripped, so `cargo package` for the
# cli resolves bramble-rs from the registry and FAILS until the matching library
# version is published. The cli is therefore validated by its own `cargo publish`
# verify step, after the library has been published below.

echo "Preflight (pre-bump): cargo package --no-verify ($LIB_DIR)"
if [[ "$DRY_RUN" == true ]]; then
    echo "Dry-run: would run package_check on $LIB_DIR"
else
    package_check "$LIB_DIR"
fi

# --- apply bumps -------------------------------------------------------------

echo
echo "Bumping versions $CURRENT_VERSION -> $VERSION"
if [[ "$DRY_RUN" == true ]]; then
    echo "Dry-run: would rewrite:"
    echo "  $LIB_CARGO [package].version"
    echo "  $CLI_CARGO [package].version + bramble-rs dependency version"
    echo "  $LIB_LOCK  (bramble-rs)"
    echo "  $CLI_LOCK  (bramble-rs, bramble-cli)"
else
    for f in "$LIB_CARGO" "$CLI_CARGO" "$LIB_LOCK" "$CLI_LOCK"; do backup "$f"; done

    bump_manifest_version "$LIB_CARGO"
    bump_manifest_version "$CLI_CARGO"
    bump_cli_dep_on_lib

    bump_lock_pkg "$LIB_LOCK" "bramble-rs"
    bump_lock_pkg "$CLI_LOCK" "bramble-rs"
    bump_lock_pkg "$CLI_LOCK" "bramble-cli"

    # verify
    [[ "$(sed -n 's/^version = "\(.*\)"/\1/p' "$LIB_CARGO" | head -1)" == "$VERSION" ]] \
        || die "library manifest bump failed"
    [[ "$(sed -n 's/^version = "\(.*\)"/\1/p' "$CLI_CARGO" | head -1)" == "$VERSION" ]] \
        || die "cli manifest bump failed"
    grep -qE "bramble-rs *= *\{[^}]*version *= *\"${VERSION}\"" "$CLI_CARGO" \
        || die "cli dependency-on-library bump failed"
fi

# --- post-bump validation ----------------------------------------------------

echo
echo "Post-bump validation: cargo package --no-verify ($LIB_DIR)"
if [[ "$DRY_RUN" == true ]]; then
    echo "Dry-run: would re-run package_check on $LIB_DIR"
else
    package_check "$LIB_DIR"
fi

# --- commit / tag / push -----------------------------------------------------

run git add "$LIB_CARGO" "$CLI_CARGO" "$LIB_LOCK" "$CLI_LOCK"
run git commit -m "chore(release): bump bramble crates to v${VERSION}"
[[ "$DRY_RUN" == false ]] && COMMIT_CREATED=true

run git tag -a "$TAG" -m "Release ${VERSION}"
run git push origin HEAD
run git push origin "$TAG"

# --- publish (library first, then binary) ------------------------------------

if [[ "$PUBLISH" == true ]]; then
    echo
    echo "Publishing $LIB_DIR to crates.io"
    run bash -c "cd '$LIB_DIR' && cargo publish"

    # `cargo publish` above blocks until bramble-rs ${VERSION} is live on the
    # registry index before returning, so the cli publish can resolve it. Poll
    # the sparse index as a belt-and-suspenders guard against index lag.
    echo
    echo "Confirming bramble-rs ${VERSION} is visible on the crates.io index..."
    if [[ "$DRY_RUN" == false ]]; then
        for _ in $(seq 1 30); do
            if curl -fsS -A "bramble-release ($(git config user.email))" \
                 "https://index.crates.io/br/am/bramble-rs" 2>/dev/null \
                 | grep -q "\"vers\":\"${VERSION}\""; then
                break
            fi
            sleep 10
        done
    fi

    echo "Publishing $CLI_DIR to crates.io"
    run bash -c "cd '$CLI_DIR' && cargo publish"
else
    echo
    echo "Skipping crates.io publish; re-run with --publish to publish both crates for v${VERSION}"
fi

echo
if [[ "$DRY_RUN" == true ]]; then
    echo "Dry-run complete"
else
    echo "Release bump complete for v${VERSION}"
fi
