#!/usr/bin/env bash

# Setup cgranges submodule first
./setup_cgranges.sh

meson setup build --buildtype=debug #--buildtype=release
ninja -C build -v
