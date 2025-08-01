project('libdeflate', 'c',
  version : '1.24',
  license : 'MIT',
  default_options : ['c_std=c99']
)

# Core library sources
libdeflate_sources = files(
  'lib/adler32.c',
  'lib/crc32.c',
  'lib/deflate_compress.c',
  'lib/deflate_decompress.c',
  'lib/gzip_compress.c',
  'lib/gzip_decompress.c',
  'lib/utils.c',
  'lib/zlib_compress.c',
  'lib/zlib_decompress.c',
)

# Add CPU-specific sources
host_cpu = host_machine.cpu_family()
if host_cpu == 'x86_64' or host_cpu == 'x86'
  libdeflate_sources += files('lib/x86/cpu_features.c')
elif host_cpu == 'aarch64' or host_cpu == 'arm'
  libdeflate_sources += files('lib/arm/cpu_features.c')
endif

# Include directory
inc_dir = include_directories('.')

# Compile arguments for optimization
libdeflate_args = []
cc = meson.get_compiler('c')

if cc.get_id() == 'gcc' or cc.get_id() == 'clang'
  # Basic optimizations
  foreach arg : ['-O3', '-fomit-frame-pointer']
    if cc.has_argument(arg)
      libdeflate_args += arg
    endif
  endforeach
  
  # CPU-specific optimizations
  if host_cpu == 'x86_64' or host_cpu == 'x86'
    foreach arg : ['-msse2', '-mpclmul']
      if cc.has_argument(arg)
        libdeflate_args += arg
      endif
    endforeach
  endif
endif

# Build the library
libdeflate_lib = static_library('deflate',
  libdeflate_sources,
  include_directories : inc_dir,
  c_args : libdeflate_args,
)

# Create dependency
libdeflate_dep = declare_dependency(
  link_with : libdeflate_lib,
  include_directories : inc_dir,
)
