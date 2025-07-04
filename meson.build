project('bramble', 'cpp', default_options : ['cpp_std=c++17'])

cc = meson.get_compiler('cpp')

# Compiler and linker options
base_flags = [
  '-Wall', '-Wextra', '-fsigned-char',
  '-D_FILE_OFFSET_BITS=64', '-D_LARGEFILE_SOURCE',
  '-fno-strict-aliasing', '-fno-rtti',
]

if get_option('buildtype') == 'debug'
  base_flags += ['-DDEBUG', '-D_DEBUG', '-DGDEBUG']
endif

add_project_arguments(base_flags, language: 'cpp')

# External static libraries
htslib_dep = dependency('', fallback: ['htslib', 'htslib_dep'])
quill_dep = dependency('quill', fallback: ['quill', 'quill_dep'])
zlib_dep = dependency('zlib', fallback: ['zlib', 'zlib_dep'])
bz2_dep = dependency('bzip2', fallback: ['bzip2', 'bzip2_dep'])
lzma_dep = dependency('liblzma', required: true)
curl_dep = dependency('libcurl', required: true)
libdeflate_dep = dependency('libdeflate', required: true)

libs = [libdeflate_dep, curl_dep, zlib_dep, bz2_dep, lzma_dep, cc.find_library('m', required: false)]

# Paths to external directories
# htslib_dir = 'htslib'
gclib_dir = 'gclib'
#htslib_sp = subproject('htslib')
htslib_inc = subproject('htslib').get_variable('htslib_inc')
htslib_lib = subproject('htslib').get_variable('htslib_lib')
gclib_dep = declare_dependency(
  include_directories: gclib_dir,
  dependencies: [htslib_dep]
)

gclib_target = custom_target(
  'gclib-target',
  output: 'gclib.stat',
  command: ['sh', '-c', 'touch gclib.stat'],
  depends: htslib_lib
)

# Include directories
inc_dirs = [htslib_inc] + include_directories('include', gclib_dir)

# Add pthread if available
if not host_machine.system().startswith('windows')
  libs += [dependency('threads')]
endif

# Source files
gclib_srcs = files(
  'gclib/GBase.cpp', 'gclib/GArgs.cpp', 'gclib/GStr.cpp', 'gclib/GSam.cpp',
  'gclib/gdna.cpp', 'gclib/codons.cpp', 'gclib/GFastaIndex.cpp',
  'gclib/GFaSeqGet.cpp', 'gclib/gff.cpp', 'gclib/GThreads.cpp'
)

main_srcs = files('src/bramble.cpp', 'src/bam.cpp', 'src/reads.cpp', 'src/tree.cpp')

# Optional sanitizers
if get_option('memcheck')
  add_project_arguments([
    '-fsanitize=address',
    '-fsanitize=undefined',
    '-fno-omit-frame-pointer'
  ], language: 'cpp')
  libs += [cc.find_library('asan', required: false), cc.find_library('ubsan', required: false)]
endif

if get_option('tsan')
  add_project_arguments([
    '-fsanitize=thread',
    '-fsanitize=undefined',
    '-fno-omit-frame-pointer'
  ], language: 'cpp')
  libs += [cc.find_library('tsan', required: false), cc.find_library('ubsan', required: false)]
endif

# Optional: disable threads
if get_option('nothreads')
  add_project_arguments('-DNOTHREADS', language: 'cpp')
else
  libs += [dependency('threads')]
endif

# Optional: statically link libstdc++
if get_option('static_cpp')
  add_project_link_arguments(['-static-libstdc++', '-static-libgcc'], language: 'cpp')
endif


executable('bramble', gclib_srcs + main_srcs + gclib_target,
  include_directories: inc_dirs,
  dependencies: [htslib_dep, quill_dep, gclib_dep] + libs,
  #link_args : ['-fsanitize=address']
)
