iri = library('iri', 'Earth/iris12.f', 'Earth/irif12.f', 'Earth/msis90.f', 'Earth/cira86_pour_iri12.f')

earth = library('earth', 'Earth/irimsis.cpp',
  include_directories : '..',
  link_with : iri)

venus = library('venus', 'Venus/venssubs_V05.f', 'Venus/pvatmos.f')

titan = library('titan', 'Titan/titandata.cpp',
  include_directories : '..')

mars = library('mars', 'Mars/datamarinerviking.cpp', 'Mars/marsatmotim.cpp', 'Mars/marsbinary.cpp',
  include_directories : '..')

planet = library('planet', 'anonym.cpp', 'earth.cpp', 'mars.cpp', 'planet.cpp', 'test.cpp', 'titan.cpp', 'venus.cpp',
  include_directories : '..',
  link_with : [earth, venus, titan, mars])
