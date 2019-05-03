subdir('math')
subdir('tinyxpath')
subdir('cppscixml')
subdir('elecimpact')

eflux = library('eflux', 'eflux/eflux.cpp')

subdir('photo')
#-------------
subdir('observ')
subdir('emit')
#----------
phflux = library('phflux', 'phflux/phflux.cpp')

dgphiv = library('dgphiv', 'proton/fortran/dgphiv.f', 'proton/fortran/expogk.f')

proton = library('proton', 'proton/proton.cpp',
  link_with : [phflux, dgphiv],
  dependencies : blas)

subdir('planet')
#--- indp

subdir('species')

subdir('chem')

protocos = library('protocos', 'protoncosmics/protocos.cpp')

subdir('atmo')
#--- Executables.

aero1d_exe = executable('aero1d', 'main/main.cpp',
  link_with : [atmo, mathfun],
  dependencies : boost_sys)

aerochem_exe = executable('aerochem', 'chemsimu/main.cpp',
  link_with : atmo,
  dependencies : boost_sys)

aerobranch = executable('aerobranch', 'computebranchphot/main.cpp',
  link_with : atmo,
  dependencies : boost_sys)
