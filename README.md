# Bramble: projection of spliced genomic alignments into transcriptomic space 🌿🫐


For the `meson` build system. Currently, system dependencies are required for 
`zlib`, `bz2`, and `lzma`, though these could be built as dependencies later
if wanted (an attempt will be made using wrapdb to build `zlib` and `bz2` if
not found on the system).

Also, to build you will need to have `meson`, `ninja`, and `autotools` installed.

To build you can either run `build.sh` directly, or do the following:


```
> meson setup build
```

then change to the build directory and build it

```
> cd build
> ninja
```

This should create the `bramble` executable!


This has been tested with `meson` 1.7 and 1.8, but likely also works with older versions.




