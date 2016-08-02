nobase_include_HEADERS = \
	hep/mc/internal/accumulator.hpp \
	hep/mc/internal/accumulator_fwd.hpp \
	hep/mc/internal/discrete_distribution.hpp \
	hep/mc/internal/generator_helper.hpp \
	hep/mc/internal/mpi_helper.hpp \
	hep/mc/distribution_parameters.hpp \
	hep/mc/distribution_result.hpp \
	hep/mc/global_configuration.hpp \
	hep/mc/integrand.hpp \
	hep/mc/mc_helper.hpp \
	hep/mc/mc_point.hpp \
	hep/mc/mc_result.hpp \
	hep/mc/multi_channel.hpp \
	hep/mc/multi_channel_callback.hpp \
	hep/mc/multi_channel_integrand.hpp \
	hep/mc/multi_channel_map.hpp \
	hep/mc/multi_channel_max_difference.hpp \
	hep/mc/multi_channel_point.hpp \
	hep/mc/multi_channel_result.hpp \
	hep/mc/multi_channel_weight_info.hpp \
	hep/mc/mpi_multi_channel.hpp \
	hep/mc/mpi_multi_channel_callback.hpp \
	hep/mc/mpi_plain.hpp \
	hep/mc/mpi_vegas.hpp \
	hep/mc/mpi_vegas_callback.hpp \
	hep/mc/plain.hpp \
	hep/mc/plain_result.hpp \
	hep/mc/projector.hpp \
	hep/mc/vegas.hpp \
	hep/mc/vegas_callback.hpp \
	hep/mc/vegas_pdf.hpp \
	hep/mc/vegas_point.hpp \
	hep/mc/vegas_result.hpp \
	hep/mc.hpp \
	hep/mc-mpi.hpp

# create a timestamp file so we can depend on it instead of depending on every
# single header
stamp-headers: $(nobase_include_HEADERS)
	echo "timestamp for headers" > $@

# build this by default
all-local: stamp-headers

MOSTLYCLEANFILES = stamp-headers
