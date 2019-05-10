catch_dep = dependency('catch2', fallback : [ 'catch2', 'catch2_dep' ])

# make it a static library because `catch.cpp` contains `main`
libcatch = static_library('catch', 'catch.cpp', dependencies : catch_dep)
libcatch_dep = declare_dependency(dependencies : catch_dep, link_with : libcatch)

tests = [
    'test_discrete_distribution',
    'test_mc_helper',
    'test_mc_point',
    'test_mc_result',
    'test_multi_channel',
    'test_multi_channel_chkpt',
    'test_non_finite_integrand',
    'test_plain',
    'test_plain_chkpt',
    'test_plain_with_distributions',
    'test_plain_with_genz_integrands',
    'test_vegas',
    'test_vegas_chkpt',
    'test_vegas_pdf',
    'test_vegas_with_genz_integrands'
]

foreach t : tests
    exe = executable(t, t + '.cpp', dependencies : [ libcatch_dep, hep_mc_dep ])
    test(t, exe, timeout : 60)
endforeach

if get_option('mpi')
    # make it a static library because `catch.cpp` contains `main`
    libcatch_mpi = static_library('catch_mpi', 'catch_mpi.cpp', dependencies : catch_dep)
    libcatch_mpi_dep = declare_dependency(dependencies : catch_dep, link_with : libcatch_mpi)

    mpi_tests = [
        'test_multi_channel',
        'test_plain',
        'test_plain_with_distributions',
        'test_vegas'
    ]

    mpiexec = find_program('mpiexec')

    foreach t : mpi_tests
        exe = executable(t + '_mpi', t + '.cpp', cpp_args : [ '-DHEP_USE_MPI' ],
            dependencies : [ hep_mc_dep, libcatch_mpi_dep, mpi_dep ])

        test(t + '_mpi', mpiexec, args : [ '-np', '2', exe ])
    endforeach
endif
