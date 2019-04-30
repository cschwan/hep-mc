#include "hep/mc/mc_point.hpp"

#include <catch2/catch.hpp>

#include <vector>

TEMPLATE_TEST_CASE("mc_point<T> creation", "", float, double, long double)
{
    using T = TestType;

    std::vector<T> vector{ T(1.0), T(2.0), T(3.0), T(4.0) };

    hep::mc_point<T> point1{vector};

    // check default argument value
    CHECK( point1.weight() == T(1.0) );

    CHECK( point1.point().at(0) == T(1.0) );
    CHECK( point1.point().at(1) == T(2.0) );
    CHECK( point1.point().at(2) == T(3.0) );
    CHECK( point1.point().at(3) == T(4.0) );

    hep::mc_point<T> point2{vector, T(2.0)};

    // check custom weight instead of the default argument
    CHECK( point2.weight() == T(2.0) );
}
