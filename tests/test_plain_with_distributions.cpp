#ifndef HEP_USE_MPI
#include "hep/mc.hpp"
#else
#include "hep/mc-mpi.hpp"
#endif

#include <catch2/catch.hpp>

#include <cmath>
#include <cstddef>
#include <limits>
#include <random>

template <typename T>
T linear_function(hep::mc_point<T> const& point, hep::projector<T>& projector)
{
    T const x = point.point().at(0);
    T const f = T(2.0) * x;

    projector.add(0, x, f);

    return f;
}

TEMPLATE_TEST_CASE("distribution accumulation", "", float, double /*, long double*/)
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
        hep::make_integrand<T>(
            linear_function<T>,
            1,
            hep::make_dist_params<T>(10, T(), T(1.0), "distribution #1")
        ),
        std::vector<std::size_t>(1, 100000),
        hep::make_plain_chkpt<T>(std::mt19937_64())
    ).results();

    auto const distribution = results.front().distributions().at(0).results();

    CHECK( distribution.at(0).non_zero_calls() == 10005 );
    CHECK( distribution.at(1).non_zero_calls() ==  9873 );
    CHECK( distribution.at(2).non_zero_calls() == 10076 );
    CHECK( distribution.at(3).non_zero_calls() ==  9920 );
    CHECK( distribution.at(4).non_zero_calls() ==  9864 );
    CHECK( distribution.at(5).non_zero_calls() == 10072 );
    CHECK( distribution.at(6).non_zero_calls() == 10115 );
    CHECK( distribution.at(7).non_zero_calls() == 10115 );
    CHECK( distribution.at(8).non_zero_calls() ==  9883 );
    CHECK( distribution.at(9).non_zero_calls() == 10077 );

    CHECK_THAT( distribution.at(0).value() , Catch::WithinULP(T(1.001907801327140122254e-01), 4) );
    CHECK_THAT( distribution.at(1).value() , Catch::WithinULP(T(2.960630313431077687274e-01), 4) );
    CHECK_THAT( distribution.at(2).value() , Catch::WithinULP(T(5.038947110406520502762e-01), 4) );
    CHECK_THAT( distribution.at(3).value() , Catch::WithinULP(T(6.946591244586934042228e-01), 4) );
    CHECK_THAT( distribution.at(4).value() , Catch::WithinULP(T(8.879290492453523290779e-01), 4) );
    CHECK_THAT( distribution.at(5).value() , Catch::WithinULP(T(1.106813171796360488561e+00), 4) );
    CHECK_THAT( distribution.at(6).value() , Catch::WithinULP(T(1.313818502654181215100e+00), 4) );
    CHECK_THAT( distribution.at(7).value() , Catch::WithinULP(T(1.516945172991019164833e+00), 4) );
    CHECK_THAT( distribution.at(8).value() , Catch::WithinULP(T(1.680209014305308981053e+00), 4) );
    CHECK_THAT( distribution.at(9).value() , Catch::WithinULP(T(1.914334532635264938269e+00), 4) );

    CHECK_THAT( distribution.at(0).error() , Catch::WithinULP(T(1.113557100141976586154e-03), 32) );
    CHECK_THAT( distribution.at(1).error() , Catch::WithinULP(T(2.886177922370993690687e-03), 32) );
    CHECK_THAT( distribution.at(2).error() , Catch::WithinULP(T(4.795659298348139786893e-03), 32) );
    CHECK_THAT( distribution.at(3).error() , Catch::WithinULP(T(6.644698084012381649117e-03), 32) );
    CHECK_THAT( distribution.at(4).error() , Catch::WithinULP(T(8.507262662791624736008e-03), 32) );
    CHECK_THAT( distribution.at(5).error() , Catch::WithinULP(T(1.047466673437628753468e-02), 32) );
    CHECK_THAT( distribution.at(6).error() , Catch::WithinULP(T(1.239863093572293159486e-02), 32) );
    CHECK_THAT( distribution.at(7).error() , Catch::WithinULP(T(1.431138887864306891215e-02), 32) );
    CHECK_THAT( distribution.at(8).error() , Catch::WithinULP(T(1.605465629957762490140e-02), 32) );
    CHECK_THAT( distribution.at(9).error() , Catch::WithinULP(T(1.809318175022538258698e-02), 32) );
}
