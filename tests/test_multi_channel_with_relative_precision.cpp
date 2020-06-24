#ifndef HEP_USE_MPI
#include "hep/mc.hpp"
#else
#include "hep/mc-mpi.hpp"
#endif

#include <catch2/catch.hpp>

#include <algorithm>
#include <cstddef>
#include <limits>
#include <random>
#include <vector>

template <typename T>
T function(hep::multi_channel_point<T> const& point)
{
    T const x = point.point().at(0);
    T const y = point.point().at(1);
    T const f = T(3.0) / T(2.0) * (x * x + y * y);

    // check that channel zero is always disabled
    CHECK( point.channel() != 0 );
    CHECK( point.point().at(0) == point.coordinates().at(0) );

    return f;
}

template <typename T>
T densities(
    std::size_t /*channel*/,
    std::vector<T> const& random_numbers,
    std::vector<T>& coordinates,
    std::vector<std::size_t> const& enabled_channels,
    std::vector<T>& densities,
    hep::multi_channel_map action
) {
    if (action == hep::multi_channel_map::calculate_densities)
    {
        for (std::size_t channel : enabled_channels)
        {
            CHECK( channel != 0 );

            densities[channel] = T(1.0);
        }

        return T(1.0);
    }

    std::copy(random_numbers.begin(), random_numbers.end(), coordinates.begin());

    return T(1.0);
}

TEMPLATE_TEST_CASE("multi_channel integration", "", float, double, long double)
{
    using T = TestType;

    // make sure that one random number suffices for all floating point types - otherwise the test
    // results depend on the numeric type
    REQUIRE( std::numeric_limits<std::mt19937_64::result_type>::digits >=
        std::numeric_limits<T>::digits );

    // check with unnormalized weights
    std::vector<T> const weights = { T(), T(1.0), T(1.0), T(1.0) };

#ifndef HEP_USE_MPI
    auto const results = hep::multi_channel(
#else
    auto const results = hep::mpi_multi_channel(
        MPI_COMM_WORLD,
#endif
        hep::make_multi_channel_integrand<T, T>(function<T>, 2, densities<T>, 2, 4),
        std::vector<std::size_t>(100, 10000),
        hep::make_multi_channel_chkpt<T>(weights, T(0.01), T(0.25), std::mt19937_64()),
#ifndef HEP_USE_MPI
        hep::callback<hep::multi_channel_chkpt<T>>{hep::callback_mode::verbose, "", T(0.001)}
#else
        hep::mpi_callback<hep::multi_channel_chkpt<T>>{hep::callback_mode::verbose, "", T(0.001)}
#endif
    ).results();

    for (auto const& result : results)
    {
        CHECK( result.calls()          == 10000 );
        CHECK( result.non_zero_calls() == 10000 );
        CHECK( result.finite_calls()   == 10000 );
    }

    auto const& result = hep::accumulate<hep::weighted_with_variance>(results.begin(),
        results.end());

    CHECK( result.error() / result.value() < T(0.001) );
    CHECK( results.size() == 41 );
}
