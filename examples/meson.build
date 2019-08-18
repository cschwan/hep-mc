examples = [
    'checkpoints',
    'checkpoints2',
    'distributions',
    'helpers',
    'multi_channel_disabled_channels',
    'multi_channel',
    'two_dimensional_distribution',
    'vegas',
    'vegas_stop_after_precision',
]

foreach e : examples
    exe = executable(e, e + '.cpp', dependencies : hep_mc_dep)
endforeach

if get_option('mpi')
    executable(
        'mpi_vegas',
        'mpi_vegas.cpp',
        include_directories : incdir,
        dependencies : [ hep_mc_dep, mpi_dep ]
    )
    
    executable(
        'mpi_multi_channel',
        'multi_channel.cpp',
        cpp_args : '-DUSE_MPI',
        dependencies : [ hep_mc_dep, mpi_dep ]
    )
endif
