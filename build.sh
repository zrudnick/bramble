#!/usr/bin/env bash

meson setup build --buildtype=debug #--buildtype=release 
ninja -C build -v
