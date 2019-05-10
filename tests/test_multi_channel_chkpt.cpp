#include "hep/mc/multi_channel_integrand.hpp"
#include "hep/mc/multi_channel.hpp"

#include <catch2/catch.hpp>

#include <cstddef>
#include <sstream>
#include <vector>

template <typename T>
T linear_function(hep::multi_channel_point<T> const& point, hep::projector<T>& projector)
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
    hep::multi_channel_map /*action*/
) {
    REQUIRE( channel == 0 );
    REQUIRE( enabled_channels.at(0) == 0 );

    coordinates.at(0) = random_numbers.at(0);
    densities.at(0) = T(1.0);

    return T(1.0);
}

TEMPLATE_TEST_CASE("empty serialization", "[multi_channel_chkpt]", float, double, long double)
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
    auto const chkpt2 = hep::make_multi_channel_chkpt<T, std::mt19937>(in);

    // serialize previous checkpoint
    std::ostringstream stream2;
    chkpt2.serialize(stream2);

    // string representation should be the same
    CHECK( stream1.str() == stream2.str() );
}

TEMPLATE_TEST_CASE("serialization", "[multi_channel_chkpt]", float, double, long double)
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
    auto const chkpt2 = hep::make_multi_channel_chkpt<T, std::mt19937>(in);

    // serialize previous checkpoint
    std::ostringstream stream2;
    chkpt2.serialize(stream2);

    // string representation should be the same
    CHECK( stream1.str() == stream2.str() );
}

TEMPLATE_TEST_CASE("empty stream construction", "[multi_channel_chkpt]", float, double, long double)
{
    using T = TestType;

    // empty stream
    std::istringstream in;

    // construct using an empty stream
    auto const chkpt1 = hep::make_multi_channel_chkpt<T, std::mt19937>(in);
    // construct using the default values
    auto const chkpt2 = hep::make_multi_channel_chkpt<T>();

    std::ostringstream out1;
    chkpt1.serialize(out1);
    std::ostringstream out2;
    chkpt2.serialize(out2);

    CHECK( out1.str() == out2.str() );
}

TEMPLATE_TEST_CASE("weight construction", "[multi_channel_chkpt]", float, double, long double)
{
    using T = TestType;

    std::vector<T> const weights = { T(0.125), T(0.125), T(0.75) };

    auto const chkpt = hep::make_multi_channel_chkpt(weights, T(0.125), T(0.25));

    CHECK( chkpt.beta() == T(0.25) );
    CHECK( chkpt.min_weight() == T(0.125) );
    CHECK( chkpt.channel_weights().at(0) == T(0.125) );
    CHECK( chkpt.channel_weights().at(1) == T(0.125) );
    CHECK( chkpt.channel_weights().at(2) == T(0.75) );
}

TEMPLATE_TEST_CASE("series", "[multi_channel_chkpt]", float, double, long double)
{
    using T = TestType;

    // check if we get the same checkpoint when we perform five iterations in one go, and when we
    // perfom first one iteration, then resume the resulting checkpoint with three iterations, and
    // finally resume it with one iteration

    auto const chkptA = hep::multi_channel(
        hep::make_multi_channel_integrand<T>(
            linear_function<T>,
            1,
            map<T>,
            1,
            1,
            hep::make_dist_params<T>(10, T(0.0), T(1.0), "distribution #1")
        ),
        std::vector<std::size_t>(5, 1000)
    );

    auto const tmp_chkpt1 = hep::multi_channel(
        hep::make_multi_channel_integrand<T>(
            linear_function<T>,
            1,
            map<T>,
            1,
            1,
            hep::make_dist_params<T>(10, T(0.0), T(1.0), "distribution #1")
        ),
        std::vector<std::size_t>(1, 1000)
    );

    auto const tmp_chkpt2 = hep::multi_channel(
        hep::make_multi_channel_integrand<T>(
            linear_function<T>,
            1,
            map<T>,
            1,
            1,
            hep::make_dist_params<T>(10, T(0.0), T(1.0), "distribution #1")
        ),
        std::vector<std::size_t>(3, 1000),
        tmp_chkpt1
    );

    auto const chkptB = hep::multi_channel(
        hep::make_multi_channel_integrand<T>(
            linear_function<T>,
            1,
            map<T>,
            1,
            1,
            hep::make_dist_params<T>(10, T(0.0), T(1.0), "distribution #1")
        ),
        std::vector<std::size_t>(1, 1000),
        tmp_chkpt2
    );

    // Check if both checkpoints are the same
    std::ostringstream stream1;
    std::ostringstream stream2;

    chkptA.serialize(stream1);
    chkptB.serialize(stream2);

    CHECK( stream1.str() == stream2.str() );
}
