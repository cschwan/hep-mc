project(
	'hep-mc',
	'cpp',
	version : '0.7',
	license : 'GPL3',
	default_options : 'cpp_std=c++11',
	meson_version : '>=0.42'
)

# must be before `doc` because doxygen depends on the headers
subdir('include')

if get_option('doxygen')
	subdir('doc')
endif

if get_option('examples')
	subdir('examples')
endif

subdir('tests')

hep_mc_dep = declare_dependency(include_directories : incdir)
