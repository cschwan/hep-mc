config = configuration_data()
config.set('srcdir', meson.current_source_dir())
config.set('top_srcdir', meson.source_root())
config.set('builddir', meson.current_build_dir())
config.set('top_builddir', meson.build_root())
config.set('PACKAGE_NAME', meson.project_name())
config.set('PACKAGE_VERSION', meson.project_version())

doxygen = find_program('doxygen')

doxyfile = configure_file(
	configuration : config,
	input : 'Doxyfile.in',
	output : 'Doxyfile'
)

datadir = join_paths(get_option('datadir'), 'doc', meson.project_name())

custom_target(
	'html',
	build_by_default : true,
	command : [ doxygen, doxyfile ],
	depend_files : headers,
	input : doxyfile,
	install : true,
	install_dir : datadir,
	output : 'html'
)
