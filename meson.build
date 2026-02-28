project('bramble', 'cpp', default_options : ['cpp_std=c++17'])

cc = meson.get_compiler('cpp')

# Compiler and linker options
base_flags = [
  '-Wall', '-Wextra', '-fsigned-char',
  '-D_FILE_OFFSET_BITS=64', '-D_LARGEFILE_SOURCE',
  '-fno-strict-aliasing'
]

if get_option('buildtype') == 'debug'
  base_flags += ['-DDEBUG', '-D_DEBUG', '-DGDEBUG']
endif

add_project_arguments(base_flags, language: 'cpp')

# External static libraries
htslib_dep = dependency('', fallback: ['htslib', 'htslib_dep'])
quill_dep = dependency('quill', fallback: ['quill', 'quill_dep'])
# unordered_dense_dep = dependency('unordered_dense', fallback: ['unordered_dense', 'unordered_dense_dep'])
zlib_dep = dependency('zlib', fallback: ['zlib', 'zlib_dep'])
# Use cc.find_library as fallback for bzip2 when dependency() fails
bz2_dep = dependency('bzip2', required: false)
if not bz2_dep.found()
  bz2_dep = cc.find_library('bz2', required: true)
endif
# Use cc.find_library as fallback for liblzma when dependency() fails
lzma_dep = dependency('liblzma', required: false)
if not lzma_dep.found()
  lzma_dep = cc.find_library('lzma', required: true)
endif
lzma_dep = dependency('liblzma', required: true)
curl_dep = dependency('libcurl', required: true)
libdeflate_dep = dependency('libdeflate', fallback: ['libdeflate', 'libdeflate_dep'], required: true)

libs = [libdeflate_dep, zlib_dep, curl_dep, bz2_dep, lzma_dep, 
        cc.find_library('m', required: false), cc.find_library('stdc++fs', required: false)]

# Paths to external directories
# htslib_dir = 'htslib'
gclib_dir = 'gclib'
ksw2_dir = 'ksw2'
#htslib_sp = subproject('htslib')
htslib_inc = subproject('htslib').get_variable('htslib_inc')
htslib_lib = subproject('htslib').get_variable('htslib_lib')
cgranges_dep = subproject('cgranges').get_variable('cgranges_dep')
unordered_dense_dep = subproject('unordered_dense').get_variable('unordered_dense_dep')
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
inc_dirs = [htslib_inc] + include_directories('include', gclib_dir) + include_directories('include', ksw2_dir)

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

ksw2_srcs = files('ksw2/ksw2_extz2_sse.cpp')

main_srcs = files('src/bramble.cpp', 'src/bam.cpp', 'src/threads.cpp',
'src/core.cpp', 'src/mates.cpp', 'src/evaluate.cpp', 'src/g2t.cpp')

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

executable('bramble', gclib_srcs + ksw2_srcs + main_srcs + gclib_target,
  include_directories: inc_dirs,
  dependencies: [htslib_dep, quill_dep, gclib_dep, cgranges_dep, unordered_dense_dep] + libs,
)
