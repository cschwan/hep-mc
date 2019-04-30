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

TEMPLATE_TEST_CASE("plain integration", "", float, double /*, long double*/)
{
    using T = TestType;

    // make sure that one random number suffices for all floating point types - otherwise the test
    // results depend on the numeric type
    REQUIRE( std::numeric_limits<std::mt19937_64::result_type>::digits >=
        std::numeric_limits<T>::digits );

#ifndef HEP_USE_MPI
    auto const results = hep::plain(
#else
    auto const results = hep::mpi_plain(
        MPI_COMM_WORLD,
#endif
        hep::make_integrand<T>(function<T>, 2),
        std::vector<std::size_t>(5, 10000),
        hep::make_plain_chkpt<T>(std::mt19937_64())
    ).results();

    for (auto const& result : results)
    {
        CHECK( result.calls()          == 10000 );
        CHECK( result.non_zero_calls() == 10000 );
        CHECK( result.finite_calls()   == 10000 );
    }

    CHECK_THAT( results.at(0).value() , Catch::WithinULP(T(1.001744429968198301186e+00), 64) );
    CHECK_THAT( results.at(1).value() , Catch::WithinULP(T(9.977324216697471673037e-01), 64) );
    CHECK_THAT( results.at(2).value() , Catch::WithinULP(T(1.006695216628880117397e+00), 64) );
    CHECK_THAT( results.at(3).value() , Catch::WithinULP(T(1.000871034592982816843e+00), 64) );
    CHECK_THAT( results.at(4).value() , Catch::WithinULP(T(1.003271572625603427396e+00), 64) );

    CHECK_THAT( results.at(0).error() , Catch::WithinULP(T(6.333298537533906387727e-03), 64) );
    CHECK_THAT( results.at(1).error() , Catch::WithinULP(T(6.255383837302876383591e-03), 64) );
    CHECK_THAT( results.at(2).error() , Catch::WithinULP(T(6.372858507247217516485e-03), 64) );
    CHECK_THAT( results.at(3).error() , Catch::WithinULP(T(6.319768554589118109626e-03), 64) );
    CHECK_THAT( results.at(4).error() , Catch::WithinULP(T(6.347542820469781915043e-03), 64) );
}
