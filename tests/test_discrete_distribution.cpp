#include "hep/mc/discrete_distribution.hpp"

#include "catch2/catch.hpp"

#include <cstddef>
#include <limits>
#include <random>
#include <vector>

TEMPLATE_TEST_CASE("discrete_distribution with no bins", "", float, double, long double)
{
    using T = TestType;

    std::mt19937_64 rng;
    std::vector<T> weights;

    hep::discrete_distribution<std::size_t, T> distribution{weights.begin(), weights.end()};

    CHECK( distribution(rng) == 0u );
    CHECK( distribution(rng) == 0u );
    CHECK( distribution(rng) == 0u );
    CHECK( distribution(rng) == 0u );
    CHECK( distribution(rng) == 0u );
}

TEMPLATE_TEST_CASE("discrete_distribution with one bin", "", float, double, long double)
{
    using T = TestType;

    std::mt19937 rng;
    std::vector<T> weights{T(1.0)};

    hep::discrete_distribution<std::size_t, T> distribution{weights.begin(), weights.end()};

    CHECK( distribution(rng) == 0u );
    CHECK( distribution(rng) == 0u );
    CHECK( distribution(rng) == 0u );
    CHECK( distribution(rng) == 0u );
    CHECK( distribution(rng) == 0u );
}

TEMPLATE_TEST_CASE("discrete_distribution with two bins", "", float, double, long double)
{
    using T = TestType;

    // make sure that one random number suffices for all floating point types - otherwise the test
    // results depend on the numeric type
    REQUIRE( std::numeric_limits<std::mt19937_64::result_type>::digits >=
        std::numeric_limits<T>::digits );

    std::mt19937_64 rng;
    std::vector<T> weights{T(1.0), T(2.0)};

    hep::discrete_distribution<std::size_t, T> distribution{weights.begin(), weights.end()};

    CHECK( distribution(rng) == 1u );
    CHECK( distribution(rng) == 0u );
    CHECK( distribution(rng) == 1u );
    CHECK( distribution(rng) == 1u );
    CHECK( distribution(rng) == 0u );
    CHECK( distribution(rng) == 1u );
    CHECK( distribution(rng) == 0u );
    CHECK( distribution(rng) == 0u );
    CHECK( distribution(rng) == 1u );
    CHECK( distribution(rng) == 1u );
    CHECK( distribution(rng) == 0u );
    CHECK( distribution(rng) == 1u );
    CHECK( distribution(rng) == 0u );
    CHECK( distribution(rng) == 1u );
    CHECK( distribution(rng) == 1u );
    CHECK( distribution(rng) == 1u );
    CHECK( distribution(rng) == 1u );
    CHECK( distribution(rng) == 1u );
    CHECK( distribution(rng) == 1u );
    CHECK( distribution(rng) == 0u );
}
