project('cgranges', 'cpp',
  version : '0.1.1',
  license : 'MIT'
)

# For C++ header-only version
cgranges_inc = include_directories('cpp')

cgranges_dep = declare_dependency(
  include_directories : cgranges_inc
)

meson.override_dependency('cgranges', cgranges_dep)
