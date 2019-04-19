catch_subproj = subproject('catch2')
catch_dep = catch_subproj.get_variable('catch2_dep')

# make it a static library because `catch.cpp` contains `main`
libcatch = static_library('catch', 'catch.cpp', dependencies : catch_dep)

libcatch_dep = declare_dependency(dependencies : catch_dep, link_with : libcatch)

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

programs = [
    'test_distributions',
    'test_mc_helper',
    'test_numerical_results',
    'test_vegas_icdf_bug',
    'test_vegas_pdf'
]

foreach p : programs
    exe = executable(
        p,
        p + '.cpp',
        cpp_args : '-DGTEST_HAS_THREAD=0',
        include_directories : incdir,
        link_with : [ libgtest, libgtest_main ]
    )

    test(p, exe)
endforeach

tests = [
    'test_discrete_distribution',
    'test_mc_point',
    'test_mc_result',
    'test_multi_channel_chkpt',
    'test_plain_chkpt',
    'test_vegas_chkpt'
]

foreach t : tests
    exe = executable(t, t + '.cpp', dependencies : [ libcatch_dep, hep_mc_dep ])
    test(t, exe)
endforeach

test_genz_integrand = executable(
    'test_genz_integrand',
    'test_genz_integrand.cpp',
    cpp_args : '-DGTEST_HAS_THREAD=0',
    include_directories : incdir,
    link_with : [ libgenz, libgtest, libgtest_main ]
)

test('test_genz_integrand', test_genz_integrand, timeout : 120)

if get_option('mpi')
    libgtest_mpi_main = static_library(
        'gtest-mpi-main',
        'gtest/gtest.h',
        'gtest_mpi_main.cpp',
        cpp_args : '-DGTEST_HAS_PTHREAD=0',
        dependencies : mpi_dep
    )

    test_mpi_numerical_results = executable(
        'test_mpi_numerical_results',
        'test_numerical_results.cpp',
        cpp_args : [ '-DGTEST_HAS_THREAD=0', '-DHEP_USE_MPI' ],
        dependencies : mpi_dep,
        include_directories : incdir,
        link_with : [ libgtest, libgtest_mpi_main ]
    )

    mpiexec = find_program('mpiexec')
    test('test_mpi_numerical_results', mpiexec,
        args : [ '-np', '2', test_mpi_numerical_results])
endif
