project('cgranges', 'c',
  default_options : ['warning_level=1', 'c_std=c99']
)

cgranges_inc = include_directories('.')

cgranges_lib = static_library('cgranges',
  'cgranges.c',
  include_directories : cgranges_inc,
)

cgranges_dep = declare_dependency(
  include_directories : cgranges_inc,
  link_with : cgranges_lib,
)
