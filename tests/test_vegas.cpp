#ifndef HEP_USE_MPI
#include "hep/mc.hpp"
#else
#include "hep/mc-mpi.hpp"
#endif

#include "catch2/catch.hpp"

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

TEMPLATE_TEST_CASE("vegas integration", "", float, double /*, long double*/)
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
        std::vector<std::size_t>(5, 10000),
        hep::make_vegas_chkpt<T>(8, T(1.5), std::mt19937_64())
    ).results();

    for (auto const& result : results)
    {
        CHECK( result.calls()          == 10000 );
        CHECK( result.non_zero_calls() == 10000 );
        CHECK( result.finite_calls()   == 10000 );
    }

    CHECK_THAT( results.at(0).value() , Catch::WithinULP(T(1.00174442996819830119e+00), 256) );
    CHECK_THAT( results.at(1).value() , Catch::WithinULP(T(1.00237281522120332683e+00), 256) );
    CHECK_THAT( results.at(2).value() , Catch::WithinULP(T(1.00096173351515851154e+00), 256) );
    CHECK_THAT( results.at(3).value() , Catch::WithinULP(T(9.98733376891225725115e-01), 256) );
    CHECK_THAT( results.at(4).value() , Catch::WithinULP(T(9.98414891829362159009e-01), 256) );

    CHECK_THAT( results.at(0).error() , Catch::WithinULP(T(6.33329853753390638773e-03), 256) );
    CHECK_THAT( results.at(1).error() , Catch::WithinULP(T(2.51428830094928175573e-03), 256) );
    CHECK_THAT( results.at(2).error() , Catch::WithinULP(T(2.21815471345527884809e-03), 256) );
    CHECK_THAT( results.at(3).error() , Catch::WithinULP(T(2.20202856630089260614e-03), 256) );
    CHECK_THAT( results.at(4).error() , Catch::WithinULP(T(2.22194654802565982583e-03), 256) );
}
