project('skfmm', 'cpp', version: '2024.03.18')

py = import('python').find_installation()

includes = include_directories('skfmm')

srcs = [
  'skfmm/fmm.cpp',
  'skfmm/heap.cpp',
  'skfmm/base_marcher.cpp',
  'skfmm/distance_marcher.cpp',
  'skfmm/travel_time_marcher.cpp',
  'skfmm/extension_velocity_marcher.cpp',
]

py.install_sources(
  ['skfmm/__init__.py', 'skfmm/pfmm.py'],
  subdir: 'skfmm',
  pure: false,
)

py.extension_module(
    'cfmm',
    srcs,
    install: true,
    subdir: 'skfmm',
    include_directories: includes,
)


py.extension_module(
    'pheap',
    ['skfmm/pheap.cpp', 'skfmm/heap.cpp'],
    install: true,
    subdir: 'skfmm',
    include_directories: includes,
)
