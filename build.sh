#!/usr/bin/env bash

meson setup build --buildtype=release #--buildtype=debug
ninja -C build -v
