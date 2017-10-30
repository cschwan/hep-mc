libgenz = static_library(
	'genz',
	'genz_integrand.hpp',
	'genz_integrand.cpp',
	include_directories : incdir
)

libgtest = static_library(
	'gtest',
	'gtest/gtest.h',
	'gtest-all.cpp',
	cpp_args : '-DGTEST_HAS_PTHREAD=0'
)

libgtest_main = static_library(
	'gtest-main',
	'gtest/gtest.h',
	'gtest_main.cpp',
	cpp_args : '-DGTEST_HAS_PTHREAD=0'
)

test_discrete_distribution = executable(
	'test_discrete_distribution',
	'test_discrete_distribution.cpp',
	cpp_args : '-DGTEST_HAS_THREAD=0',
	include_directories : incdir,
	link_with : [ libgtest, libgtest_main ]
)

test_distributions = executable(
	'test_distributions',
	'test_distributions.cpp',
	cpp_args : '-DGTEST_HAS_THREAD=0',
	include_directories : incdir,
	link_with : [ libgtest, libgtest_main ]
)

test_genz_integrand = executable(
	'test_genz_integrand',
	'test_genz_integrand.cpp',
	cpp_args : '-DGTEST_HAS_THREAD=0',
	include_directories : incdir,
	link_with : [ libgenz, libgtest, libgtest_main ]
)

test_mc_helper = executable(
	'test_mc_helper',
	'test_mc_helper.cpp',
	cpp_args : '-DGTEST_HAS_THREAD=0',
	include_directories : incdir,
	link_with : [ libgtest, libgtest_main ]
)

test_mc_point = executable(
	'test_mc_point',
	'test_mc_point.cpp',
	cpp_args : '-DGTEST_HAS_THREAD=0',
	include_directories : incdir,
	link_with : [ libgtest, libgtest_main ]
)

test_mc_result = executable(
	'test_mc_result',
	'test_mc_result.cpp',
	cpp_args : '-DGTEST_HAS_THREAD=0',
	include_directories : incdir,
	link_with : [ libgtest, libgtest_main ]
)

test_numerical_results = executable(
	'test_numerical_results',
	'test_numerical_results.cpp',
	cpp_args : '-DGTEST_HAS_THREAD=0',
	include_directories : incdir,
	link_with : [ libgtest, libgtest_main ]
)

test_vegas_icdf_bug = executable(
	'test_vegas_icdf_bug',
	'test_vegas_icdf_bug.cpp',
	cpp_args : '-DGTEST_HAS_THREAD=0',
	include_directories : incdir,
	link_with : [ libgtest, libgtest_main ]
)

test_vegas_pdf = executable(
	'test_vegas_pdf',
	'test_vegas_pdf.cpp',
	cpp_args : '-DGTEST_HAS_THREAD=0',
	include_directories : incdir,
	link_with : [ libgtest, libgtest_main ]
)

if get_option('enable-mpi')
	libgtest_mpi_main = static_library(
		'gtest-mpi-main',
		'gtest/gtest.h',
		'gtest_mpi_main.cpp',
		cpp_args : '-DGTEST_HAS_PTHREAD=0',
		dependencies : mpidep
	)

	test_mpi_numerical_results = executable(
		'test_mpi_numerical_results',
		'test_numerical_results.cpp',
		cpp_args : [ '-DGTEST_HAS_THREAD=0', '-DHEP_USE_MPI' ],
		dependencies : mpidep,
		include_directories : incdir,
		link_with : [ libgtest, libgtest_mpi_main ]
	)

	config = configuration_data()
	config.set('abs_builddir', meson.current_build_dir())
	shell_script = configure_file(input : 'test_mpi_numerical_results.sh.in',
		output : 'test_mpi_numerical_results.sh', configuration : config)

	test('test_mpi_numerical_results', shell_script)
endif

test('test_discrete_distribution', test_discrete_distribution)
test('test_distributions', test_distributions)
test('test_genz_integrand', test_genz_integrand, timeout : 120)
test('test_mc_helper', test_mc_helper)
test('test_mc_point', test_mc_point)
test('test_mc_result', test_mc_result)
test('test_numerical_results', test_numerical_results)
test('test_vegas_icdf_bug', test_vegas_icdf_bug)
test('test_vegas_pdf', test_vegas_pdf)
