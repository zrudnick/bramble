#!/usr/bin/env bash

# Setup cgranges submodule first
./setup_cgranges.sh

meson setup build --buildtype=release #--buildtype=release 
ninja -C build -v
