#include "hep/mc/mc_result.hpp"

#include <catch2/catch.hpp>

#include <sstream>

TEMPLATE_TEST_CASE("mc_result<T> construction", "", float, double, long double)
{
    using T = TestType;

    hep::mc_result<T> result{100, 99, 98, T(100.0), T(100.0)};

    CHECK( result.calls()          == 100u );
    CHECK( result.non_zero_calls() ==  99u );
    CHECK( result.finite_calls()   ==  98u );
    CHECK( result.sum()            == T(100.0) );
    CHECK( result.sum_of_squares() == T(100.0) );

    CHECK( result.value()          == T(1.0) );
    CHECK( result.variance()       == T() );
    CHECK( result.error()          == T() );
}

TEMPLATE_TEST_CASE("mc_result<T> serialization", "", float, double, long double)
{
    using T = TestType;

    std::stringstream stream;

    // write a result to a stream read it back in
    hep::mc_result<T>{100, 99, 98, T(100.0), T(100.0)}.serialize(stream);
    hep::mc_result<T> result{stream};

    // Is it the same result?

    CHECK( result.calls()          == 100u );
    CHECK( result.non_zero_calls() ==  99u );
    CHECK( result.finite_calls()   ==  98u );
    CHECK( result.sum()            == T(100.0) );
    CHECK( result.sum_of_squares() == T(100.0) );

    CHECK( result.value()          == T(1.0) );
    CHECK( result.variance()       == T() );
    CHECK( result.error()          == T() );
}

TEMPLATE_TEST_CASE("create_result function", "", float, double, long double)
{
    using T = TestType;

    auto result = hep::create_result(100, 99, 98, T(1.0), T());

    CHECK( result.calls()          == 100u );
    CHECK( result.non_zero_calls() ==  99u );
    CHECK( result.finite_calls()   ==  98u );
    CHECK( result.sum()            == T(100.0) );
    CHECK( result.sum_of_squares() == T(100.0) );

    CHECK( result.value()          == T(1.0) );
    CHECK( result.variance()       == T() );
    CHECK( result.error()          == T() );
}
