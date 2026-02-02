project('cgranges', 'c', version: '0.1.1')

cgranges_inc = include_directories('.', 'cpp')

cgranges_lib = static_library('cgranges',
  'cgranges.c',
  include_directories: cgranges_inc
)

cgranges_dep = declare_dependency(
  include_directories: cgranges_inc,
  link_with: cgranges_lib
)
