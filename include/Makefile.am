nobase_include_HEADERS = \
	hep/mc/default_vegas_parallelizer.hpp \
	hep/mc/mc_point.hpp \
	hep/mc/mc_result.hpp \
	hep/mc/plain.hpp \
	hep/mc/vegas.hpp \
	hep/mc.hpp

# create a timestamp file so we can depend on it instead of depending on every
# single header
stamp-headers: $(nobase_include_HEADERS)
	echo "timestamp for headers" > $@

# build this by default
all-local: stamp-headers

MOSTLYCLEANFILES = stamp-headers
