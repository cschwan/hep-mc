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
    'test_mc_helper',
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
    'test_distributions_with_plain',
    'test_mc_point',
    'test_mc_result',
    'test_multi_channel',
    'test_multi_channel_chkpt',
    'test_plain',
    'test_plain_chkpt',
    'test_vegas',
    'test_vegas_chkpt',
    'test_vegas_pdf'
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
    # make it a static library because `catch.cpp` contains `main`
    libcatch_mpi = static_library('catch_mpi', 'catch_mpi.cpp', dependencies : catch_dep)
    libcatch_mpi_dep = declare_dependency(dependencies : catch_dep, link_with : libcatch_mpi)

    mpi_tests = [
        'test_distributions_with_plain',
        'test_multi_channel',
        'test_plain',
        'test_vegas'
    ]

    mpiexec = find_program('mpiexec')

    foreach t : mpi_tests
        exe = executable(t + '_mpi', t + '.cpp', cpp_args : [ '-DHEP_USE_MPI' ],
            dependencies : [ hep_mc_dep, libcatch_mpi_dep, mpi_dep ])

        test(t + '_mpi', mpiexec, args : [ '-np', '2', exe ])
    endforeach
endif
