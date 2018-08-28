#include "gtest/gtest.h"

#include "hep/mc/distribution_parameters.hpp"
#include "hep/mc/distribution_result.hpp"
#include "hep/mc/mc_helper.hpp"

#include <cmath>
#include <limits>
#include <vector>

typedef testing::Types<float, double, long double> MyT;
template <typename T> class McHelper : public testing::Test { };
TYPED_TEST_CASE(McHelper, MyT);

TYPED_TEST(McHelper, CumulativeResult0WithOne)
{
    typedef TypeParam T;

    std::vector<hep::mc_result<T>> one_result = {
        hep::mc_result<T>(1000, 1000, 1000, T(10.0), T(20.0))
    };
    auto result = hep::accumulate<hep::weighted_with_variance>(
        one_result.begin(), one_result.end());

    EXPECT_EQ( one_result.front().calls() , result.calls() );
    EXPECT_EQ( one_result.front().non_zero_calls() , result.calls() );
    EXPECT_EQ( one_result.front().finite_calls() , result.calls() );
    EXPECT_NEAR( one_result.front().value() , result.value() , T(1e-10) );
    EXPECT_NEAR( one_result.front().error() , result.error() , T(1e-10) );
}

TYPED_TEST(McHelper, CumulativeResult1WithOne)
{
    typedef TypeParam T;

    std::vector<hep::mc_result<T>> one_result = {
        hep::mc_result<T>(1000, 1000, 1000, T(10.0), T(20.0))
    };
    auto result = hep::accumulate<hep::weighted_equally>(one_result.begin(),
        one_result.end());

    EXPECT_EQ( one_result.front().calls() , result.calls() );
    EXPECT_EQ( one_result.front().non_zero_calls() , result.calls() );
    EXPECT_EQ( one_result.front().finite_calls() , result.calls() );
    EXPECT_NEAR( one_result.front().value() , result.value() , T(1e-10) );
    EXPECT_NEAR( one_result.front().error() , result.error() , T(1e-10) );
}

TYPED_TEST(McHelper, CumulativeResult0WithTwo)
{
    typedef TypeParam T;

    // two times the same result
    std::vector<hep::mc_result<T>> two_results = {
        hep::mc_result<T>(100, 100, 100, T(100.0), T(200.0)),
        hep::mc_result<T>(100, 100, 100, T(100.0), T(200.0))
    };
    auto result = hep::accumulate<hep::weighted_with_variance>(
        two_results.begin(), two_results.end());

    // should give two times the calls
    EXPECT_EQ( 2 * two_results.front().calls() , result.calls() );
    EXPECT_EQ( 2 * two_results.front().non_zero_calls() , result.calls() );
    EXPECT_EQ( 2 * two_results.front().finite_calls() , result.calls() );
    // the same result and
    EXPECT_NEAR( two_results.front().value() , result.value() , T(1e-10) );
    // an error reduces by 1/sqrt(2)
    EXPECT_NEAR( two_results.front().error() / std::sqrt(T(2.0)) ,
        result.error() , T(1e-9) );
}

TYPED_TEST(McHelper, CumulativeResult1WithTwo)
{
    typedef TypeParam T;

    // two times the same result
    std::vector<hep::mc_result<T>> two_results = {
        hep::mc_result<T>(100, 100, 100, T(100.0), T(200.0)),
        hep::mc_result<T>(100, 100, 100, T(100.0), T(200.0))
    };
    auto result = hep::accumulate<hep::weighted_equally>(two_results.begin(),
        two_results.end());

    // should give two times the calls
    EXPECT_EQ( 2 * two_results.front().calls() , result.calls() );
    EXPECT_EQ( 2 * two_results.front().non_zero_calls() , result.calls() );
    EXPECT_EQ( 2 * two_results.front().finite_calls() , result.calls() );
    // the same result and
    EXPECT_NEAR( two_results.front().value() , result.value() , T(1e-10) );
    // zero error
    EXPECT_EQ( T() , result.error() );
}

TYPED_TEST(McHelper, CumulativeResult1WithThree)
{
    typedef TypeParam T;

    // two times the same result
    std::vector<hep::mc_result<T>> two_results = {
        hep::mc_result<T>(100, 100, 100, T(100.0), T(200.0)),
        hep::mc_result<T>(100, 100, 100, T(100.0), T(200.0)),
        hep::mc_result<T>(100, 100, 100, T(100.0), T(200.0)),
    };
    auto result = hep::accumulate<hep::weighted_equally>(two_results.begin(),
        two_results.end());

    // should give three times the calls
    EXPECT_EQ( 3 * two_results.front().calls() , result.calls() );
    EXPECT_EQ( 3 * two_results.front().non_zero_calls() , result.calls() );
    EXPECT_EQ( 3 * two_results.front().finite_calls() , result.calls() );
    // the same result and
    EXPECT_NEAR( two_results.front().value() , result.value() , T(1e-10) );
    // zero error
    EXPECT_EQ( T() , result.error() );
}

TYPED_TEST(McHelper, ChiSquareDof0WithZero)
{
    typedef TypeParam T;

    std::vector<hep::mc_result<T>> zero_results;
    T result = hep::chi_square_dof<hep::weighted_with_variance>(
        zero_results.begin(), zero_results.end());

    EXPECT_NEAR( T() , result , T(1e-10) );
}

TYPED_TEST(McHelper, ChiSquareDof1WithZero)
{
    typedef TypeParam T;

    std::vector<hep::mc_result<T>> zero_results;
    T result = hep::chi_square_dof<hep::weighted_equally>(zero_results.begin(),
        zero_results.end());

    EXPECT_NEAR( T() , result , T(1e-10) );
}

TYPED_TEST(McHelper, ChiSquareDof0WithOne)
{
    typedef TypeParam T;

    std::vector<hep::mc_result<T>> one_result = {
        hep::mc_result<T>(1000, 1000, 1000, T(10.0), T(20.0))
    };
    T result = hep::chi_square_dof<hep::weighted_with_variance>(
        one_result.begin(), one_result.end());

    EXPECT_EQ( std::numeric_limits<T>::infinity() , result );
}

TYPED_TEST(McHelper, ChiSquareDof1WithOne)
{
    typedef TypeParam T;

    std::vector<hep::mc_result<T>> one_result = {
        hep::mc_result<T>(1000, 1000, 1000, T(10.0), T(20.0))
    };
    T result = hep::chi_square_dof<hep::weighted_equally>(one_result.begin(),
        one_result.end());

    EXPECT_EQ( std::numeric_limits<T>::infinity() , result );
}

TYPED_TEST(McHelper, ChiSquareDof0WithTwo)
{
    typedef TypeParam T;

    // two times the same result
    std::vector<hep::mc_result<T>> two_results = {
        hep::mc_result<T>(100, 100, 100, T(100.0), T(200.0)),
        hep::mc_result<T>(100, 100, 100, T(100.0), T(200.0))
    };
    T result = hep::chi_square_dof<hep::weighted_with_variance>(
        two_results.begin(), two_results.end());

    EXPECT_NEAR( T() , result , T(1e-10) );
}

TYPED_TEST(McHelper, ChiSquareDof1WithTwo)
{
    typedef TypeParam T;

    // two times the same result
    std::vector<hep::mc_result<T>> two_results = {
        hep::mc_result<T>(100, 100, 100, T(100.0), T(200.0)),
        hep::mc_result<T>(100, 100, 100, T(100.0), T(200.0))
    };
    T result = hep::chi_square_dof<hep::weighted_equally>(two_results.begin(),
        two_results.end());

    EXPECT_NEAR( T() , result , T(1e-10) );
}

TYPED_TEST(McHelper, DistributionAccumulator)
{
    using T = TypeParam;
    using std::sqrt;

    std::vector<hep::distribution_result<T>> it0 = {
        hep::distribution_result<T>{
            hep::make_dist_params(3, T(), T(3.0)),
            std::vector<hep::mc_result<T>>{
                hep::create_result<T>(100, 100, 100, T(1.0), T(0.1)),
                hep::create_result<T>(100, 100, 100, T(2.0), T(0.2)),
                hep::create_result<T>(100, 100, 100, T(3.0), T(0.3))
            }
        },
        hep::distribution_result<T>{
            hep::make_dist_params(2, T(), T(2.0)),
            std::vector<hep::mc_result<T>>{
                hep::create_result<T>(100, 100, 100, T(4.0), T(0.4)),
                hep::create_result<T>(100, 100, 100, T(8.0), T(0.8))
            }
        }
    };

    auto const ir = hep::create_result(100, 100, 100, T(5.0), T(0.5));

    std::vector<hep::plain_result<T>> results = {
        hep::plain_result<T>{it0, ir.non_zero_calls(), ir.finite_calls(), ir.calls(), ir.sum(), ir.sum_of_squares()},
        hep::plain_result<T>{it0, ir.non_zero_calls(), ir.finite_calls(), ir.calls(), ir.sum(), ir.sum_of_squares()}
    };

    auto const result = hep::accumulate<hep::weighted_with_variance>(
        results.begin(), results.end());

    // check integrated result
    EXPECT_NEAR( T(5.0), result.value(), T(1e-10) );
    EXPECT_NEAR( sqrt(T(0.125)), result.error(), T(1e-8) );
    EXPECT_EQ( result.calls(), 200u );
    EXPECT_EQ( result.non_zero_calls(), 200u );
    EXPECT_EQ( result.finite_calls(), 200u );

    // check parameters of both distributions
    EXPECT_EQ( result.distributions().at(0).parameters().bins_x(), 3u );
    EXPECT_EQ( result.distributions().at(0).parameters().x_min(), T() );
//    EXPECT_EQ( result.distributions().at(0).parameters().x_max(), T(3.0) );
    EXPECT_EQ( result.distributions().at(1).parameters().bins_x(), 2u );
    EXPECT_EQ( result.distributions().at(1).parameters().x_min(), T() );
//    EXPECT_EQ( result.distributions().at(1).parameters().x_max(), T(2.0) );

    // check first distribution
    EXPECT_NEAR( T(1.0), result.distributions().at(0).results().at(0).value(), T(1e-10) );
    EXPECT_NEAR( sqrt(T(0.005)), result.distributions().at(0).results().at(0).error(), T(1e-10) );
    EXPECT_EQ( 200u, result.distributions().at(0).results().at(0).calls() );
    EXPECT_EQ( 200u, result.distributions().at(0).results().at(0).non_zero_calls() );
    EXPECT_EQ( 200u, result.distributions().at(0).results().at(0).finite_calls() );
    EXPECT_NEAR( T(2.0), result.distributions().at(0).results().at(1).value(), T(1e-10) );
    EXPECT_NEAR( sqrt(T(0.02)), result.distributions().at(0).results().at(1).error(), T(1e-10) );
    EXPECT_EQ( 200u, result.distributions().at(0).results().at(1).calls() );
    EXPECT_EQ( 200u, result.distributions().at(0).results().at(1).non_zero_calls() );
    EXPECT_EQ( 200u, result.distributions().at(0).results().at(1).finite_calls() );
    EXPECT_NEAR( T(3.0), result.distributions().at(0).results().at(2).value(), T(1e-10) );
    EXPECT_NEAR( sqrt(T(0.045)), result.distributions().at(0).results().at(2).error(), T(1e-10) );
    EXPECT_EQ( 200u, result.distributions().at(0).results().at(2).calls() );
    EXPECT_EQ( 200u, result.distributions().at(0).results().at(2).non_zero_calls() );
    EXPECT_EQ( 200u, result.distributions().at(0).results().at(2).finite_calls() );

    // check second distribution
    EXPECT_NEAR( T(4.0), result.distributions().at(1).results().at(0).value(), T(1e-10) );
    EXPECT_NEAR( sqrt(T(0.08)), result.distributions().at(1).results().at(0).error(), T(1e-10) );
    EXPECT_EQ( 200u, result.distributions().at(1).results().at(0).calls() );
    EXPECT_EQ( 200u, result.distributions().at(1).results().at(0).non_zero_calls() );
    EXPECT_EQ( 200u, result.distributions().at(1).results().at(0).finite_calls() );
    EXPECT_NEAR( T(8.0), result.distributions().at(1).results().at(1).value(), T(1e-6) );
    EXPECT_NEAR( sqrt(T(0.32)), result.distributions().at(1).results().at(1).error(), T(1e-10) );
    EXPECT_EQ( 200u, result.distributions().at(1).results().at(1).calls() );
    EXPECT_EQ( 200u, result.distributions().at(1).results().at(1).non_zero_calls() );
    EXPECT_EQ( 200u, result.distributions().at(1).results().at(1).finite_calls() );
}
