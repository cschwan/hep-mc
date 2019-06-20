#ifndef HEP_USE_MPI
#include "hep/mc.hpp"
#else
#include "hep/mc-mpi.hpp"
#endif

#include <catch2/catch.hpp>

#include <cstddef>
#include <limits>
#include <random>
#include <vector>

template <typename T>
T function(hep::mc_point<T> const& point)
{
    T const x = point.point().at(0);
    T const y = point.point().at(1);
    T const f = T(3.0) / T(2.0) * (x * x + y * y);

    return f;
}

TEMPLATE_TEST_CASE("stop after precision", "", float, double, long double)
{
    using T = TestType;

    // make sure that one random number suffices for all floating point types - otherwise the test
    // results depend on the numeric type
    REQUIRE( std::numeric_limits<std::mt19937_64::result_type>::digits >=
        std::numeric_limits<T>::digits );

#ifndef HEP_USE_MPI
    auto const results = hep::vegas(
#else
    auto const results = hep::mpi_vegas(
        MPI_COMM_WORLD,
#endif
        hep::make_integrand<T>(function<T>, 2),
        std::vector<std::size_t>(100, 10000),
        hep::make_vegas_chkpt<T>(8, T(1.5), std::mt19937_64()),
#ifndef HEP_USE_MPI
        hep::callback<hep::vegas_chkpt<T>>{hep::callback_mode::verbose, "", T(0.001)}
#else
        hep::mpi_callback<hep::vegas_chkpt<T>>{hep::callback_mode::verbose, "", T(0.001)}
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
    CHECK( results.size() == 7 );
}
