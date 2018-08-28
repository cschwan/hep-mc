project(
    'hep-mc',
    'cpp',
    version : '0.8',
    license : 'GPL3',
    default_options : 'cpp_std=c++11',
    meson_version : '>=0.40'
)

# must be before `doc` because doxygen depends on the headers
subdir('include')

if get_option('doxygen')
    subdir('doc')
endif

pkg = import('pkgconfig')
pkg.generate(
    description : 'A C++ Template Library for Monte Carlo Integration',
    name : meson.project_name(),
    version : meson.project_version()
)

hep_mc_dep = declare_dependency(include_directories : incdir)

if get_option('mpi')
    mpi_dep = dependency('mpi', language : 'cpp')
endif
    
if get_option('examples')
    subdir('examples')
endif

subdir('tests')
