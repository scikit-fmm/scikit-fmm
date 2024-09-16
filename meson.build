project('skfmm', 'cpp', version: '2024.09.16')

py = import('python').find_installation()

incdir_numpy = meson.get_external_property('numpy-include-dir', 'not-given')
if incdir_numpy == 'not-given'
  incdir_numpy = run_command(py,
    [
      '-c',
      '''
import os
import numpy as np
try:
  incdir = os.path.relpath(np.get_include())
except Exception:
  incdir = np.get_include()
print(incdir)
'''
    ],
    check: true
  ).stdout().strip()
endif

inc_np = include_directories(incdir_numpy)
np_dep = declare_dependency(include_directories: inc_np)
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
  ['skfmm/__init__.py',
   'skfmm/pfmm.py',
   'skfmm/heap.py'],
  subdir: 'skfmm',
  pure: false,
)

py.extension_module(
    'cfmm',
    srcs,
    install: true,
    subdir: 'skfmm',
    include_directories: [includes, inc_np],
    dependencies: [np_dep],
)


py.extension_module(
    'pheap',
    ['skfmm/pheap.cpp',
     'skfmm/heap.cpp'],
    install: true,
    subdir: 'skfmm',
    include_directories: includes,
)
