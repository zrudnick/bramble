#!/bin/bash
# This script ensures cgranges submodule has the necessary meson.build file

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CGRANGES_DIR="$SCRIPT_DIR/subprojects/cgranges"
MESON_BUILD_SRC="$SCRIPT_DIR/subprojects/packagefiles/cgranges/meson.build"
MESON_BUILD_DST="$CGRANGES_DIR/meson.build"

# Initialize submodule if not already done
if [ ! -f "$CGRANGES_DIR/cgranges.c" ]; then
    echo "Initializing cgranges submodule..."
    git -C "$SCRIPT_DIR" submodule update --init "$CGRANGES_DIR"
fi

# Copy meson.build if it doesn't exist
if [ ! -f "$MESON_BUILD_DST" ]; then
    echo "Copying meson.build to cgranges submodule..."
    cp "$MESON_BUILD_SRC" "$MESON_BUILD_DST"
fi

echo "cgranges setup complete"