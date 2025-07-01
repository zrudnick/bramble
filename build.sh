#!/usr/bin/env bash

cd htslib
./configure
make -j8
cd ..
meson setup builddir
cd builddir
ninja
