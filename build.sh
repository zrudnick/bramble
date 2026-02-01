#!/usr/bin/env bash

meson setup build --buildtype=release #--buildtype=release 
ninja -C build -v
