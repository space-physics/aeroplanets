project('trans++', 'cpp', 'fortran',
  default_options : ['default_library=static', 'buildtype=release', 'warning_level=0'])

blas = dependency('blas')

boost_hdr = dependency('boost')
boost_regex = dependency('boost', modules : 'regex')
boost_fs = dependency('boost', modules : 'filesystem')
boost_sys = dependency('boost', modules : 'system')

subdir('src')
