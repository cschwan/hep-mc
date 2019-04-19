#include "hep/mc/multi_channel_integrand.hpp"
#include "hep/mc/multi_channel.hpp"

#include "catch2/catch.hpp"

#include <cstddef>
#include <sstream>
#include <vector>

template <typename T>
T linear_function(hep::mc_point<T> const& point, hep::projector<T>& projector)
{
    T const x = point.point().at(0);
    T const v = x;

    projector.add(0, x, v);

    return v;
}

template <typename T>
T map(
    std::size_t channel,
    std::vector<T> const& random_numbers,
    std::vector<T>& coordinates,
    std::vector<std::size_t> const& enabled_channels,
    std::vector<T>& densities,
    hep::multi_channel_map
) {
    assert( channel == 0 );
    assert( enabled_channels.at(0) == 0 );

    coordinates.at(0) = random_numbers.at(0);
    densities.at(0) = T(1.0);

    return T(1.0);
}

TEMPLATE_TEST_CASE("multi_channel_chkpt<T> empty serialization", "", float, double, long double)
{
    using T = TestType;

    // create a checkpoint using five iterations with PLAIN
    auto chkpt1 = hep::make_multi_channel_chkpt<T>();

    chkpt1.channels(1);

    CHECK( chkpt1.beta()       == T(0.25) );
    CHECK( chkpt1.min_weight() == T() );

    // serialize checkpoint
    std::ostringstream stream1;
    chkpt1.serialize(stream1);

    // unserialize checkpoint
    auto in = std::istringstream(stream1.str());
    auto const chkpt2 = decltype (chkpt1) {in};

    // serialize previous checkpoint
    std::ostringstream stream2;
    chkpt2.serialize(stream2);

    // string representation should be the same
    CHECK( stream1.str() == stream2.str() );
}

TEMPLATE_TEST_CASE("multi_channel_chkpt<T> serialization", "", float, double, long double)
{
    using T = TestType;

    // create a checkpoint using five iterations with PLAIN
    auto const chkpt1 = hep::multi_channel(
        hep::make_multi_channel_integrand<T>(
            linear_function<T>,
            1,
            map<T>,
            1,
            1,
            hep::make_dist_params<T>(10, T(0.0), T(1.0), "distribution name")
        ),
        std::vector<std::size_t>(5, 1000)
    );

    // serialize checkpoint
    std::ostringstream stream1;
    chkpt1.serialize(stream1);

    // unserialize checkpoint
    auto in = std::istringstream(stream1.str());
    auto const chkpt2 = decltype (chkpt1) {in};

    // serialize previous checkpoint
    std::ostringstream stream2;
    chkpt2.serialize(stream2);

    // string representation should be the same
    CHECK( stream1.str() == stream2.str() );
}
