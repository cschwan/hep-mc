project(
	'hep-mc',
	'cpp',
	version : '0.6',
	license : 'GPL-3',
	default_options : 'cpp_std=c++11',
	meson_version : '>=0.42'
)

mpidep = dependency('mpi', language : 'cpp')
incdir = include_directories('include')

if get_option('enable-doxygen')
	subdir('doc')
endif

if get_option('enable-examples')
	subdir('examples')
endif

subdir('include')
subdir('tests')

hep_mc_dep = declare_dependency(include_directories : incdir)
