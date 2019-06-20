#include "hep/mc.hpp"

#include <catch2/catch.hpp>

#include <sstream>
#include <string>

TEMPLATE_TEST_CASE("distribution parameters", "", float, double, long double)
{
    using T = TestType;

    hep::distribution_parameters<T> params1{8, T(1.0), (5.0), "distribution #1"};

    CHECK( params1.bins_x() == 8 );
    CHECK( params1.bins_y() == 1 );
    CHECK( params1.name() == "distribution #1" );
    CHECK( params1.x_min() == T(1.0) );
    CHECK( params1.y_min() == T() );
    CHECK( params1.bin_size_x() == T(0.5) );
    CHECK( params1.bin_size_y() == T(1.0) );

    hep::distribution_parameters<T> params2{8, 2, T(1.0), (5.0), T(2.0), T(4.0), "distribution #2"};

    CHECK( params2.bins_x() == 8 );
    CHECK( params2.bins_y() == 2 );
    CHECK( params2.name() == "distribution #2" );
    CHECK( params2.x_min() == T(1.0) );
    CHECK( params2.y_min() == T(2.0) );
    CHECK( params2.bin_size_x() == T(0.5) );
    CHECK( params2.bin_size_y() == T(1.0) );
}

TEMPLATE_TEST_CASE("serialization/deserialization", "", float, double, long double)
{
    using T = TestType;

    hep::distribution_parameters<T> tmp1{8, T(1.0), (5.0), "distribution #1"};

    std::stringstream stream1;
    tmp1.serialize(stream1);
    hep::distribution_parameters<T> params1{stream1};

    CHECK( params1.bins_x() == 8 );
    CHECK( params1.bins_y() == 1 );
    CHECK( params1.name() == "distribution #1" );
    CHECK( params1.x_min() == T(1.0) );
    CHECK( params1.y_min() == T() );
    CHECK( params1.bin_size_x() == T(0.5) );
    CHECK( params1.bin_size_y() == T(1.0) );

    hep::distribution_parameters<T> tmp2{8, 2, T(1.0), (5.0), T(2.0), T(4.0), "distribution #2"};

    std::stringstream stream2;
    tmp2.serialize(stream2);
    hep::distribution_parameters<T> params2{stream2};

    CHECK( params2.bins_x() == 8 );
    CHECK( params2.bins_y() == 2 );
    CHECK( params2.name() == "distribution #2" );
    CHECK( params2.x_min() == T(1.0) );
    CHECK( params2.y_min() == T(2.0) );
    CHECK( params2.bin_size_x() == T(0.5) );
    CHECK( params2.bin_size_y() == T(1.0) );
}
