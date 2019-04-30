#include "hep/mc/integrand.hpp"
#include "hep/mc/plain.hpp"
#include "hep/mc/vegas.hpp"

#include "genz_integrand.hpp"

#include <catch2/catch.hpp>

#include <cmath>
#include <cstddef>
#include <tuple>
#include <vector>

TEST_CASE("plain with genz integrands", "")
{
    using std::fabs;
    using T = double;

    // these parameters were picked from the diploma thesis "High-Dimensional Numerical Integration
    // on Parallel Computers" from Rudolf Schuerer, chapter 2. First parameter is the number of
    // dimensions, second parameter is the maximally allowed deviation from the reference result in
    // multiples of the error the monte carlo integrator returns.
    auto const data = GENERATE( values({
//        std::make_tuple(  2, 2, 110, 3, genz::oscillatory),
        std::make_tuple(  3, 2, 110, 3, genz::oscillatory),
        std::make_tuple(  5, 1, 110, 3, genz::oscillatory),
        std::make_tuple( 10, 1, 110, 3, genz::oscillatory),
        std::make_tuple( 15, 2, 110, 3, genz::oscillatory),
        std::make_tuple( 20, 2, 110, 3, genz::oscillatory),
        std::make_tuple( 25, 2, 110, 3, genz::oscillatory),
        std::make_tuple( 30, 1, 110, 3, genz::oscillatory),
        std::make_tuple( 40, 1, 110, 3, genz::oscillatory),
        std::make_tuple( 60, 1, 110, 3, genz::oscillatory),
        std::make_tuple( 80, 1, 110, 3, genz::oscillatory),
        std::make_tuple(100, 1, 110, 3, genz::oscillatory),
        std::make_tuple(  2, 2, 600, 4, genz::product_peak),
        std::make_tuple(  3, 1, 600, 4, genz::product_peak),
        std::make_tuple(  5, 1, 600, 4, genz::product_peak),
        std::make_tuple( 10, 2, 600, 4, genz::product_peak),
        std::make_tuple( 15, 1, 600, 4, genz::product_peak),
        std::make_tuple( 20, 2, 600, 4, genz::product_peak),
        std::make_tuple( 25, 3, 600, 4, genz::product_peak),
        std::make_tuple( 30, 3, 600, 4, genz::product_peak),
//        std::make_tuple( 40, 1, 600, 4, genz::product_peak),   // TOO SMALL
//        std::make_tuple( 60, 1, 600, 4, genz::product_peak),   // TOO SMALL
        std::make_tuple( 80, 1, 600, 4, genz::product_peak),
        std::make_tuple(100, 1, 600, 4, genz::product_peak),
        std::make_tuple(  2, 3, 600, 4, genz::corner_peak),
        std::make_tuple(  3, 2, 600, 4, genz::corner_peak),
        std::make_tuple(  5, 2, 600, 4, genz::corner_peak),
        std::make_tuple( 10, 1, 600, 4, genz::corner_peak),
        std::make_tuple( 15, 2, 600, 4, genz::corner_peak),
        std::make_tuple( 20, 2, 600, 4, genz::corner_peak),
//        std::make_tuple( 25, 1, 600, 4, genz::corner_peak),
//        std::make_tuple( 30, 1, 600, 4, genz::corner_peak),
//        std::make_tuple( 40, 1, 600, 4, genz::corner_peak),
//        std::make_tuple( 60, 1, 600, 4, genz::corner_peak),
//        std::make_tuple( 80, 1, 600, 4, genz::corner_peak),
//        std::make_tuple(100, 1, 600, 4, genz::corner_peak),
        std::make_tuple(  2, 2, 100, 2, genz::gaussian),
        std::make_tuple(  3, 1, 100, 2, genz::gaussian),
        std::make_tuple(  5, 2, 100, 2, genz::gaussian),
        std::make_tuple( 10, 2, 100, 2, genz::gaussian),
        std::make_tuple( 15, 1, 100, 2, genz::gaussian),
        std::make_tuple( 20, 2, 100, 2, genz::gaussian),
        std::make_tuple( 25, 2, 100, 2, genz::gaussian),
        std::make_tuple( 30, 3, 100, 2, genz::gaussian),
        std::make_tuple( 40, 1, 100, 2, genz::gaussian),
        std::make_tuple( 60, 2, 100, 2, genz::gaussian),
        std::make_tuple( 80, 1, 100, 2, genz::gaussian),
        std::make_tuple(100, 1, 100, 2, genz::gaussian),
        std::make_tuple(  2, 2, 150, 4, genz::c0_function),
        std::make_tuple(  3, 1, 150, 4, genz::c0_function),
        std::make_tuple(  5, 1, 150, 4, genz::c0_function),
        std::make_tuple( 10, 2, 150, 4, genz::c0_function),
        std::make_tuple( 15, 1, 150, 4, genz::c0_function),
        std::make_tuple( 20, 1, 150, 4, genz::c0_function),
        std::make_tuple( 25, 2, 150, 4, genz::c0_function),
        std::make_tuple( 30, 3, 150, 4, genz::c0_function),
        std::make_tuple( 40, 1, 150, 4, genz::c0_function),
        std::make_tuple( 60, 2, 150, 4, genz::c0_function),
        std::make_tuple( 80, 1, 150, 4, genz::c0_function),
        std::make_tuple(100, 1, 150, 4, genz::c0_function),
        std::make_tuple(  2, 3, 100, 4, genz::discontinuous),
        std::make_tuple(  3, 1, 100, 4, genz::discontinuous),
        std::make_tuple(  5, 3, 100, 4, genz::discontinuous),
        std::make_tuple( 10, 2, 100, 4, genz::discontinuous),
        std::make_tuple( 15, 1, 100, 4, genz::discontinuous),
        std::make_tuple( 20, 1, 100, 4, genz::discontinuous),
        std::make_tuple( 25, 4, 100, 4, genz::discontinuous),  // PLAIN = 1
        std::make_tuple( 30, 1, 100, 4, genz::discontinuous),
        std::make_tuple( 40, 1, 100, 4, genz::discontinuous),
        std::make_tuple( 60, 1, 100, 4, genz::discontinuous),
        std::make_tuple( 80, 2, 100, 4, genz::discontinuous),
        std::make_tuple(100, 2, 100, 4, genz::discontinuous)
    }) );

    int dimension = std::get<0>(data);
    int deviation = std::get<1>(data);
    T diff        = std::get<2>(data);
    T comp        = std::get<3>(data) * T(0.5);
    auto type     = std::get<4>(data);

    INFO( "dimension=" << dimension );
    INFO( "deviation=" << deviation );
    INFO( "diff=" << diff );
    INFO( "comp=" << comp );
    INFO( "type=" << type );

    auto params    = genz::parameters<T>(dimension, diff, comp, T(0.1));
    auto integrand = genz::integrand<T>(type, params.affective(), params.unaffective());

    auto const result = hep::vegas(
        hep::make_integrand<T>(integrand, dimension),
        std::vector<std::size_t>(1, 100000)
    ).results().front();

    // approximation should lie with the error interval
    CHECK( fabs(integrand.reference_result() - result.value()) <= deviation * result.error() );

    if (result.value() != T())
    {
        // relative error should not be larger than 15%
        CHECK( result.error() / fabs(result.value()) < T(0.15) );
    }
}
