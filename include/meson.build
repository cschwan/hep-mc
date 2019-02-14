incdir = include_directories('.')

headers1 = [
    'hep/mc/accumulator.hpp',
    'hep/mc/accumulator_fwd.hpp',
    'hep/mc/chkpt.hpp',
    'hep/mc/discrete_distribution.hpp',
    'hep/mc/distribution_parameters.hpp',
    'hep/mc/distribution_result.hpp',
    'hep/mc/generator_helper.hpp',
    'hep/mc/global_configuration.hpp',
    'hep/mc/integrand.hpp',
    'hep/mc/mc_helper.hpp',
    'hep/mc/mc_point.hpp',
    'hep/mc/mc_result.hpp',
    'hep/mc/mpi_helper.hpp',
    'hep/mc/mpi_multi_channel.hpp',
    'hep/mc/mpi_multi_channel_callback.hpp',
    'hep/mc/mpi_plain.hpp',
    'hep/mc/mpi_plain_callback.hpp',
    'hep/mc/mpi_vegas.hpp',
    'hep/mc/mpi_vegas_callback.hpp',
    'hep/mc/multi_channel.hpp',
    'hep/mc/multi_channel_callback.hpp',
    'hep/mc/multi_channel_integrand.hpp',
    'hep/mc/multi_channel_map.hpp',
    'hep/mc/multi_channel_max_difference.hpp',
    'hep/mc/multi_channel_point.hpp',
    'hep/mc/multi_channel_result.hpp',
    'hep/mc/multi_channel_weight_info.hpp',
    'hep/mc/plain.hpp',
    'hep/mc/plain_callback.hpp',
    'hep/mc/plain_chkpt.hpp',
    'hep/mc/plain_result.hpp',
    'hep/mc/projector.hpp',
    'hep/mc/vegas.hpp',
    'hep/mc/vegas_callback.hpp',
    'hep/mc/vegas_pdf.hpp',
    'hep/mc/vegas_point.hpp',
    'hep/mc/vegas_result.hpp',
]

headers2 = [
    'hep/mc.hpp',
    'hep/mc-mpi.hpp',
]

headers = files(headers1, headers2)

install_headers(headers1, subdir : 'hep/mc')
install_headers(headers2, subdir : 'hep')
