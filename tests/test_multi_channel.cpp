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

TEMPLATE_TEST_CASE("multi_channel integration", "", float, double /*, long double*/)
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
        hep::make_multi_channel_integrand<T>(function<T>, 2, densities<T>, 2, 4),
        std::vector<std::size_t>(5, 10000),
        hep::make_multi_channel_chkpt<T>(weights, T(0.01), T(0.25), std::mt19937_64())
    ).results();

    for (auto const& result : results)
    {
        CHECK( result.calls()          == 10000 );
        CHECK( result.non_zero_calls() == 10000 );
        CHECK( result.finite_calls()   == 10000 );
    }

    CHECK_THAT( results.at(0).value() , Catch::WithinULP(T(1.000889500971818506844e+00), 128) );
    CHECK_THAT( results.at(1).value() , Catch::WithinULP(T(1.007360827957602174439e+00), 128) );
    CHECK_THAT( results.at(2).value() , Catch::WithinULP(T(9.993263442793500384736e-01), 128) );
    CHECK_THAT( results.at(3).value() , Catch::WithinULP(T(9.991312981071902387764e-01), 128) );
    CHECK_THAT( results.at(4).value() , Catch::WithinULP(T(1.008099287068812917106e+00), 128) );

    CHECK_THAT( results.at(0).error() , Catch::WithinULP(T(6.259819425011487338953e-03), 128) );
    CHECK_THAT( results.at(1).error() , Catch::WithinULP(T(6.326650126493282871129e-03), 128) );
    CHECK_THAT( results.at(2).error() , Catch::WithinULP(T(6.296140552521163745576e-03), 128) );
    CHECK_THAT( results.at(3).error() , Catch::WithinULP(T(6.332853475743208759244e-03), 128) );
    CHECK_THAT( results.at(4).error() , Catch::WithinULP(T(6.348087527478404719798e-03), 128) );
}
