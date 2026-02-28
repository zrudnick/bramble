project('ksw2', 'cpp', default_options: ['cpp_std=c++17'])

ksw2_inc = include_directories('.')

ksw2_lib = static_library('ksw2', 'ksw2_extz2_sse.cpp',
  include_directories: ksw2_inc)

ksw2_dep = declare_dependency(
  include_directories: ksw2_inc,
  link_with: ksw2_lib)
