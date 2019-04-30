#include "hep/mc/distribution_parameters.hpp"
#include "hep/mc/distribution_result.hpp"
#include "hep/mc/mc_helper.hpp"

#include <catch2/catch.hpp>

#include <cmath>
#include <limits>
#include <vector>

TEMPLATE_TEST_CASE("accumulate_result with one result", "", float, double /*, long double*/)
{
    using T = TestType;

    std::vector<hep::mc_result<T>> vector{hep::mc_result<T>(1000, 999, 998, T(1000.0), T(2000.0))};

    auto result1 = hep::accumulate<hep::weighted_with_variance>(vector.begin(), vector.end());

    CHECK( result1.calls()          == 1000 );
    CHECK( result1.non_zero_calls() ==  999 );
    CHECK( result1.finite_calls()   ==  998 );

    CHECK_THAT( result1.value() , Catch::WithinULP(T(1.0), 4) );
    CHECK_THAT( result1.error() , Catch::WithinULP(vector.front().error(), 4) );

    auto result2 = hep::accumulate<hep::weighted_equally>(vector.begin(), vector.end());

    CHECK( result2.calls()          == 1000 );
    CHECK( result2.non_zero_calls() ==  999 );
    CHECK( result2.finite_calls()   ==  998 );

    CHECK_THAT( result2.value() , Catch::WithinULP(T(1.0), 4) );
    CHECK_THAT( result2.error() , Catch::WithinULP(vector.front().error(), 4) );
}

TEMPLATE_TEST_CASE("accumulate_result with two results", "", float, double /*, long double*/)
{
    using std::sqrt;
    using T = TestType;

    std::vector<hep::mc_result<T>> vector{
        hep::mc_result<T>(100, 99, 98, T(100.0), T(10000.0)),
        hep::mc_result<T>(100, 99, 98, T(200.0), T(40000.0))
    };

    auto result1 = hep::accumulate<hep::weighted_with_variance>(vector.begin(), vector.end());

    CHECK( result1.calls()          == 200 );
    CHECK( result1.non_zero_calls() == 198 );
    CHECK( result1.finite_calls()   == 196 );

    CHECK_THAT( result1.sum()            , Catch::WithinULP(T(240.0), 4) );
    CHECK_THAT( result1.sum_of_squares() , Catch::WithinULP(T(32128.0), 4) );

    auto result2 = hep::accumulate<hep::weighted_equally>(vector.begin(), vector.end());

    CHECK( result2.calls()          == 200 );
    CHECK( result2.non_zero_calls() == 198 );
    CHECK( result2.finite_calls()   == 196 );

    CHECK_THAT( result2.sum()            , Catch::WithinULP(T(300.0), 4) );
    CHECK_THAT( result2.sum_of_squares() , Catch::WithinULP(T(10400.0), 4) );
}

TEMPLATE_TEST_CASE("chi_square_dof with zero results", "", float, double /*, long double*/)
{
    using T = TestType;

    std::vector<hep::mc_result<T>> results;

    T result1 = hep::chi_square_dof<hep::weighted_with_variance>(results.begin(), results.end());

    CHECK( result1 == T() );

    T result2 = hep::chi_square_dof<hep::weighted_equally>(results.begin(), results.end());

    CHECK( result2 == T() );
}

TEMPLATE_TEST_CASE("chi_square_dof with one result", "", float, double /*, long double*/)
{
    using T = TestType;

    std::vector<hep::mc_result<T>> vector{hep::mc_result<T>(100, 99, 98, T(100.0), T(10000.0))};

    T result1 = hep::chi_square_dof<hep::weighted_with_variance>(vector.begin(), vector.end());

    CHECK( result1 == std::numeric_limits<T>::infinity() );

    T result2 = hep::chi_square_dof<hep::weighted_equally>(vector.begin(), vector.end());

    CHECK( result2 == std::numeric_limits<T>::infinity() );
}

TEMPLATE_TEST_CASE("chi_square_dof with two results", "", float, double /*, long double*/)
{
    using T = TestType;

    std::vector<hep::mc_result<T>> vector{
        hep::mc_result<T>(100, 99, 98, T(100.0), T(10000.0)),
        hep::mc_result<T>(100, 99, 98, T(200.0), T(40000.0))
    };

    T result1 = hep::chi_square_dof<hep::weighted_with_variance>(vector.begin(), vector.end());

    CHECK_THAT( result1, Catch::WithinULP(T(0.2), 4) );

    T result2 = hep::chi_square_dof<hep::weighted_equally>(vector.begin(), vector.end());

    CHECK_THAT( result2, Catch::WithinULP(T(0.3125), 4) );
}
