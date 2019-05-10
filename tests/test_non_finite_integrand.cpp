#include "hep/mc/distribution_parameters.hpp"
#include "hep/mc/integrand.hpp"
#include "hep/mc/plain.hpp"

#include <catch2/catch.hpp>

#include <cstddef>
#include <limits>
#include <vector>

template <typename T>
T function(hep::mc_point<T> const& point, hep::projector<T>& projector)
{
    static bool return_nan = true;

    REQUIRE( std::numeric_limits<T>::has_quiet_NaN );
    REQUIRE( std::numeric_limits<T>::has_infinity );

    T value;

    if (return_nan)
    {
        value = std::numeric_limits<T>::quiet_NaN();
    }
    else
    {
        value = std::numeric_limits<T>::infinity();
    }

    return_nan = !return_nan;

    REQUIRE( point.point().size() <= 2 );

    if (point.point().size() == 1)
    {
        projector.add(0, point.point().at(0), value);
    }
    else if (point.point().size() == 2)
    {
        projector.add(0, point.point().at(0), point.point().at(1), value);
    }

    return value;
}

template <typename T>
T function_no_dist(hep::mc_point<T> const& /*point*/)
{
    static bool return_nan = true;

    REQUIRE( std::numeric_limits<T>::has_quiet_NaN );
    REQUIRE( std::numeric_limits<T>::has_infinity );

    T value;

    if (return_nan)
    {
        value = std::numeric_limits<T>::quiet_NaN();
    }
    else
    {
        value = std::numeric_limits<T>::infinity();
    }

    return_nan = !return_nan;

    return value;
}

TEMPLATE_TEST_CASE("non-finite 1D integrand w/o distributions", "", float, double, long double)
{
    using T = TestType;

    // random numbers don't really matter here
    auto const results = hep::plain(
        hep::make_integrand<T>(function_no_dist<T>, 1),
        std::vector<std::size_t>(5, 10000)
    ).results();

    for (auto const& result : results)
    {
        CHECK( result.calls()          == 10000 );
        CHECK( result.non_zero_calls() == 10000 ); // NaN or Inf are non-zero
        CHECK( result.finite_calls()   == 0 );

        CHECK( result.value() == T() );
        CHECK( result.error() == T() );

        CHECK( result.distributions().empty() );
    }
}

TEMPLATE_TEST_CASE("non-finite 1D integrand w/ distributions", "", float, double, long double)
{
    using T = TestType;

    // random numbers don't really matter here
    auto const results = hep::plain(
        hep::make_integrand<T>(
            function<T>,
            1,
            hep::distribution_parameters<T>{2, T(), T(1.0), "distribution #1"}
        ),
        std::vector<std::size_t>(5, 10000)
    ).results();

    for (auto const& result : results)
    {
        CHECK( result.calls()          == 10000 );
        CHECK( result.non_zero_calls() == 10000 ); // NaN or Inf are non-zero
        CHECK( result.finite_calls()   == 0 );

        CHECK( result.value() == T() );
        CHECK( result.error() == T() );

        for (auto const& bin : result.distributions().at(0).results())
        {
            CHECK( bin.calls()          == 10000 );
//            CHECK( bin.non_zero_calls() == 10000 ); // NaN or Inf are non-zero
            CHECK( bin.finite_calls()   == 0 );

            CHECK( bin.value() == T() );
            CHECK( bin.error() == T() );
        }
    }
}

TEMPLATE_TEST_CASE("non-finite 2D integrand w/ distributions", "", float, double, long double)
{
    using T = TestType;

    // random numbers don't really matter here
    auto const results = hep::plain(
        hep::make_integrand<T>(
            function<T>,
            2,
            hep::distribution_parameters<T>{2, 2, T(), T(1.0), T(), T(1.0), "distribution #1"}
        ),
        std::vector<std::size_t>(5, 10000)
    ).results();

    for (auto const& result : results)
    {
        CHECK( result.calls()          == 10000 );
        CHECK( result.non_zero_calls() == 10000 ); // NaN or Inf are non-zero
        CHECK( result.finite_calls()   == 0 );

        CHECK( result.value() == T() );
        CHECK( result.error() == T() );

        for (auto const& bin : result.distributions().at(0).results())
        {
            CHECK( bin.calls()          == 10000 );
//            CHECK( bin.non_zero_calls() == 10000 ); // NaN or Inf are non-zero
            CHECK( bin.finite_calls()   == 0 );

            CHECK( bin.value() == T() );
            CHECK( bin.error() == T() );
        }
    }
}
