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
T integrand(hep::mc_point<T> const& point, hep::projector<T>& projector)
{
    // two-dimensional distribution
    T const x = point.point().at(0);
    T const y = point.point().at(1);
    T const f = T(4.0) * x * y;

    projector.add(0, x, y, f);

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
        hep::make_integrand<T, T>(
            integrand<T>,
            2,
            hep::distribution_parameters<T>(4, 6, T(), T(1.0), T(-0.25), T(1.25), "distribution #1")
        ),
        std::vector<std::size_t>(1, 100000),
        hep::make_plain_chkpt<T>(std::mt19937_64())
    ).results();

    auto const& bins = results.front().distributions().at(0).results();

    REQUIRE( bins.size() == 24 );

    CHECK( bins.at( 0).non_zero_calls() == 0 );
    CHECK( bins.at( 1).non_zero_calls() == 0 );
    CHECK( bins.at( 2).non_zero_calls() == 0 );
    CHECK( bins.at( 3).non_zero_calls() == 0 );
    CHECK( bins.at(20).non_zero_calls() == 0 );
    CHECK( bins.at(21).non_zero_calls() == 0 );
    CHECK( bins.at(22).non_zero_calls() == 0 );
    CHECK( bins.at(23).non_zero_calls() == 0 );

    CHECK( bins.at( 0).value() == T() );
    CHECK( bins.at( 1).value() == T() );
    CHECK( bins.at( 2).value() == T() );
    CHECK( bins.at( 3).value() == T() );
    CHECK( bins.at(20).value() == T() );
    CHECK( bins.at(21).value() == T() );
    CHECK( bins.at(22).value() == T() );
    CHECK( bins.at(23).value() == T() );

    CHECK( bins.at( 0).error() == T() );
    CHECK( bins.at( 1).error() == T() );
    CHECK( bins.at( 2).error() == T() );
    CHECK( bins.at( 3).error() == T() );
    CHECK( bins.at(20).error() == T() );
    CHECK( bins.at(21).error() == T() );
    CHECK( bins.at(22).error() == T() );
    CHECK( bins.at(23).error() == T() );

    CHECK( bins.at( 4).non_zero_calls() == 6275 );
    CHECK( bins.at( 5).non_zero_calls() == 6224 );
    CHECK( bins.at( 6).non_zero_calls() == 6232 );
    CHECK( bins.at( 7).non_zero_calls() == 6231 );
    CHECK( bins.at( 8).non_zero_calls() == 6305 );
    CHECK( bins.at( 9).non_zero_calls() == 6192 );
    CHECK( bins.at(10).non_zero_calls() == 6283 );
    CHECK( bins.at(11).non_zero_calls() == 6137 );
    CHECK( bins.at(12).non_zero_calls() == 6159 );
    CHECK( bins.at(13).non_zero_calls() == 6279 );
    CHECK( bins.at(14).non_zero_calls() == 6293 );
    CHECK( bins.at(15).non_zero_calls() == 6348 );
    CHECK( bins.at(16).non_zero_calls() == 6256 );
    CHECK( bins.at(17).non_zero_calls() == 6121 );
    CHECK( bins.at(18).non_zero_calls() == 6366 );
    CHECK( bins.at(19).non_zero_calls() == 6299 );

    CHECK_THAT( bins.at( 4).value() , Catch::WithinULP(T(6.273918338116714768038e-02), 4) );
    CHECK_THAT( bins.at( 5).value() , Catch::WithinULP(T(1.867342319479317085991e-01), 4) );
    CHECK_THAT( bins.at( 6).value() , Catch::WithinULP(T(3.138616415792370977086e-01), 4) );
    CHECK_THAT( bins.at( 7).value() , Catch::WithinULP(T(4.384800514103403067221e-01), 4) );
    CHECK_THAT( bins.at( 8).value() , Catch::WithinULP(T(1.873282568236345351528e-01), 4) );
    CHECK_THAT( bins.at( 9).value() , Catch::WithinULP(T(5.577452223219004951598e-01), 4) );
    CHECK_THAT( bins.at(10).value() , Catch::WithinULP(T(9.433639157449172221248e-01), 4) );
    CHECK_THAT( bins.at(11).value() , Catch::WithinULP(T(1.287612637922220693455e+00), 4) );
    CHECK_THAT( bins.at(12).value() , Catch::WithinULP(T(3.085258313637559095103e-01), 4) );
    CHECK_THAT( bins.at(13).value() , Catch::WithinULP(T(9.437871147027528395961e-01), 4) );
    CHECK_THAT( bins.at(14).value() , Catch::WithinULP(T(1.573752362633750532107e+00), 4) );
    CHECK_THAT( bins.at(15).value() , Catch::WithinULP(T(2.219259903350687465964e+00), 4) );
    CHECK_THAT( bins.at(16).value() , Catch::WithinULP(T(4.373832443213468370473e-01), 4) );
    CHECK_THAT( bins.at(17).value() , Catch::WithinULP(T(1.283771627755091892685e+00), 4) );
    CHECK_THAT( bins.at(18).value() , Catch::WithinULP(T(2.224626313677139543851e+00), 4) );
    CHECK_THAT( bins.at(19).value() , Catch::WithinULP(T(3.082318732307455763902e+00), 4) );

    CHECK_THAT( bins.at( 4).error() , Catch::WithinULP(T(1.040068422453591732379e-03), 32) );
    CHECK_THAT( bins.at( 5).error() , Catch::WithinULP(T(2.724169446037201852213e-03), 32) );
    CHECK_THAT( bins.at( 6).error() , Catch::WithinULP(T(4.508790882942817602586e-03), 32) );
    CHECK_THAT( bins.at( 7).error() , Catch::WithinULP(T(6.290108729694922768934e-03), 32) );
    CHECK_THAT( bins.at( 8).error() , Catch::WithinULP(T(2.722541915010574052165e-03), 32) );
    CHECK_THAT( bins.at( 9).error() , Catch::WithinULP(T(7.131866145076577196537e-03), 32) );
    CHECK_THAT( bins.at(10).error() , Catch::WithinULP(T(1.183446478125283318626e-02), 32) );
    CHECK_THAT( bins.at(11).error() , Catch::WithinULP(T(1.629297564703037392253e-02), 32) );
    CHECK_THAT( bins.at(12).error() , Catch::WithinULP(T(4.448947366842394582378e-03), 32) );
    CHECK_THAT( bins.at(13).error() , Catch::WithinULP(T(1.184523282620369712076e-02), 32) );
    CHECK_THAT( bins.at(14).error() , Catch::WithinULP(T(1.948394607802795694613e-02), 32) );
    CHECK_THAT( bins.at(15).error() , Catch::WithinULP(T(2.725102411867845891823e-02), 32) );
    CHECK_THAT( bins.at(16).error() , Catch::WithinULP(T(6.258049791224178284770e-03), 32) );
    CHECK_THAT( bins.at(17).error() , Catch::WithinULP(T(1.627300481355644373461e-02), 32) );
    CHECK_THAT( bins.at(18).error() , Catch::WithinULP(T(2.727638627454345610897e-02), 32) );
    CHECK_THAT( bins.at(19).error() , Catch::WithinULP(T(3.786197424909581959234e-02), 32) );

    auto const& x = mid_points_x(results.front().distributions().at(0));
    auto const& y = mid_points_y(results.front().distributions().at(0));

    CHECK( x.at( 0) == T(0.125) ); CHECK( y.at( 0) == T(-0.125) );
    CHECK( x.at( 1) == T(0.375) ); CHECK( y.at( 1) == T(-0.125) );
    CHECK( x.at( 2) == T(0.625) ); CHECK( y.at( 2) == T(-0.125) );
    CHECK( x.at( 3) == T(0.875) ); CHECK( y.at( 3) == T(-0.125) );
    CHECK( x.at( 4) == T(0.125) ); CHECK( y.at( 4) == T( 0.125) );
    CHECK( x.at( 5) == T(0.375) ); CHECK( y.at( 5) == T( 0.125) );
    CHECK( x.at( 6) == T(0.625) ); CHECK( y.at( 6) == T( 0.125) );
    CHECK( x.at( 7) == T(0.875) ); CHECK( y.at( 7) == T( 0.125) );
    CHECK( x.at( 8) == T(0.125) ); CHECK( y.at( 8) == T( 0.375) );
    CHECK( x.at( 9) == T(0.375) ); CHECK( y.at( 9) == T( 0.375) );
    CHECK( x.at(10) == T(0.625) ); CHECK( y.at(10) == T( 0.375) );
    CHECK( x.at(11) == T(0.875) ); CHECK( y.at(11) == T( 0.375) );
    CHECK( x.at(12) == T(0.125) ); CHECK( y.at(12) == T( 0.625) );
    CHECK( x.at(13) == T(0.375) ); CHECK( y.at(13) == T( 0.625) );
    CHECK( x.at(14) == T(0.625) ); CHECK( y.at(14) == T( 0.625) );
    CHECK( x.at(15) == T(0.875) ); CHECK( y.at(15) == T( 0.625) );
    CHECK( x.at(16) == T(0.125) ); CHECK( y.at(16) == T( 0.875) );
    CHECK( x.at(17) == T(0.375) ); CHECK( y.at(17) == T( 0.875) );
    CHECK( x.at(18) == T(0.625) ); CHECK( y.at(18) == T( 0.875) );
    CHECK( x.at(19) == T(0.875) ); CHECK( y.at(19) == T( 0.875) );
    CHECK( x.at(20) == T(0.125) ); CHECK( y.at(20) == T( 1.125) );
    CHECK( x.at(21) == T(0.375) ); CHECK( y.at(21) == T( 1.125) );
    CHECK( x.at(22) == T(0.625) ); CHECK( y.at(22) == T( 1.125) );
    CHECK( x.at(23) == T(0.875) ); CHECK( y.at(23) == T( 1.125) );
}
