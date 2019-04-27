examples = [
    'checkpoints',
    'distributions_example',
    'helper_example',
    'multi_channel_disabled_channels',
    'multi_channel_example',
    'two_dimensional_distribution',
    'vegas_example',
    'vegas_read_pdf',
    'vegas_write_pdf',
    'vegas_stop_after_precision',
]

foreach e : examples
    exe = executable(e, e + '.cpp', dependencies : hep_mc_dep)
endforeach

if get_option('mpi')
    executable(
        'mpi_vegas_example',
        'mpi_vegas_example.cpp',
        include_directories : incdir,
        dependencies : [ hep_mc_dep, mpi_dep ]
    )
    
    executable(
        'mpi_multi_channel_example',
        'multi_channel_example.cpp',
        cpp_args : '-DUSE_MPI',
        dependencies : [ hep_mc_dep, mpi_dep ]
    )
endif
