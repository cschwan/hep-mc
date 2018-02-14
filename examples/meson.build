executable(
	'distributions_example',
	'distributions_example.cpp',
	include_directories : incdir
)

executable(
	'helper_example',
	'helper_example.cpp',
	include_directories : incdir
)

executable(
	'multi_channel_disabled_channels',
	'multi_channel_disabled_channels.cpp',
	include_directories : incdir
)

executable(
	'multi_channel_example',
	'multi_channel_example.cpp',
	include_directories : incdir
)

executable(
	'vegas_example',
	'vegas_example.cpp',
	include_directories : incdir
)

executable(
	'vegas_read_pdf',
	'vegas_read_pdf.cpp',
	include_directories : incdir
)

executable(
	'vegas_write_pdf',
	'vegas_write_pdf.cpp',
	include_directories : incdir
)

executable(
	'vegas_stop_after_precision',
	'vegas_stop_after_precision.cpp',
	include_directories : incdir
)

if get_option('mpi')

mpi = dependency('mpi', language : 'cpp')

executable(
	'mpi_vegas_example',
	'mpi_vegas_example.cpp',
	include_directories : incdir,
	dependencies : mpi
)

executable(
	'mpi_multi_channel_example',
	'multi_channel_example.cpp',
	cpp_args : '-DUSE_MPI',
	include_directories : incdir,
	dependencies : mpi
)

endif
