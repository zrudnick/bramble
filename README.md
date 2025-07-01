# Bramble: genome to transcriptome coordinate conversion 🌿🫐


For the `meson` build system. Currently, system dependencies are required for 
`zlib`, `bz2`, `lzma`, `curl`, and `libdeflate`, though these could be versioned 
later if wanted (an attempt will be made using wrapdb to build `zlib` and `bz2` if
not found on the system).


First, directly build htslib (in-source):

```
> cd htslib
> ./configure && make -j8
```

Then, back out to the root directory

```
> meson setup builddir
```

The build it

```
> cd builddir
> ninja
```

This should create the `bramble` executable!




