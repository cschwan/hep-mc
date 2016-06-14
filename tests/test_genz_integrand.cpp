#include "gtest/gtest.h"
#include "genz_integrand.hpp"

#include "hep/mc/integrand.hpp"
#include "hep/mc/plain.hpp"
#include "hep/mc/vegas.hpp"

#include <cstddef>
#include <tuple>
#include <vector>

// these parameters were picked from the diploma thesis "High-Dimensional
// Numerical Integration on Parallel Computers" from Rudolf Schuerer, chapter 2.
// First parameter is the number of dimensions, second parameter is the
// maximally allowed deviation from the reference result in multiples of the
// error the monte carlo integrator returns.
auto test_data = {
	std::make_tuple(  2, 2, 110, 3, genz::oscillatory),
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
//	std::make_tuple( 40, 1, 600, 4, genz::product_peak),   // TOO SMALL
//	std::make_tuple( 60, 1, 600, 4, genz::product_peak),   // TOO SMALL
	std::make_tuple( 80, 1, 600, 4, genz::product_peak),
	std::make_tuple(100, 1, 600, 4, genz::product_peak),
	std::make_tuple(  2, 3, 600, 4, genz::corner_peak),
	std::make_tuple(  3, 2, 600, 4, genz::corner_peak),
	std::make_tuple(  5, 2, 600, 4, genz::corner_peak),
	std::make_tuple( 10, 1, 600, 4, genz::corner_peak),
	std::make_tuple( 15, 2, 600, 4, genz::corner_peak),
	std::make_tuple( 20, 2, 600, 4, genz::corner_peak),
//	std::make_tuple( 25, 1, 600, 4, genz::corner_peak),
//	std::make_tuple( 30, 1, 600, 4, genz::corner_peak),
//	std::make_tuple( 40, 1, 600, 4, genz::corner_peak),
//	std::make_tuple( 60, 1, 600, 4, genz::corner_peak),
//	std::make_tuple( 80, 1, 600, 4, genz::corner_peak),
//	std::make_tuple(100, 1, 600, 4, genz::corner_peak),
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
	std::make_tuple(100, 2, 100, 4, genz::discontinuous),
};

using test_data_type = std::tuple<int, int, int, int, genz::integrand_type>;

template <typename T>
void check_plain_integrator(test_data_type data)
{
	int dimension = std::get<0>(data);
	int deviation = std::get<1>(data);
	T diff        = std::get<2>(data);
	T comp        = std::get<3>(data) * T(0.5);
	auto type     = std::get<4>(data);
	T limit       = T(0.1);
	int calls     = 100000;

	auto params    = genz::parameters<T>(dimension, diff, comp, limit);
	auto integrand = genz::integrand<T>(type, params.affective(),
		params.unaffective());

	auto result = hep::plain(
		hep::make_integrand<T>(integrand, dimension),
		calls
	);

	// approximation should lie with the error interval
	EXPECT_NEAR( integrand.reference_result(), result.value(), deviation *
		result.error() );

	// relative error should not be larger than 5%
	if (result.value() != T())
	{
		EXPECT_GT( T(0.15) , result.error() / result.value() );
	}
}

template <typename T>
void check_vegas_integrator(test_data_type data)
{
	int dimension = std::get<0>(data);
	int deviation = std::get<1>(data);
	T diff        = std::get<2>(data);
	T comp        = std::get<3>(data) * T(0.5);
	auto type     = std::get<4>(data);
	T limit       = T(0.1);
	int calls     = 100000;

	auto params    = genz::parameters<T>(dimension, diff, comp, limit);
	auto integrand = genz::integrand<T>(type, params.affective(),
		params.unaffective());

	auto results = hep::vegas(
		hep::make_integrand<T>(integrand, dimension),
		std::vector<std::size_t>(10, calls / 10)
	);
	auto result = hep::cumulative_result1(results.begin(), results.end());

	// approximation should lie with the error interval
	EXPECT_NEAR( integrand.reference_result(), result.value(), deviation *
		result.error() );

	// relative error should not be larger than 5%
	if (result.value() != T())
	{
		EXPECT_GT( T(0.07) , result.error() / result.value() );
	}
}

template <typename T>
struct GenzTest : public ::testing::TestWithParam<test_data_type>
{
	using TypeParam = T;
};

using GenzTestWithDouble = GenzTest<double>;

TEST_P(GenzTestWithDouble, CheckPlainIntegrator)
{
	check_plain_integrator<TypeParam>(GetParam());
}

TEST_P(GenzTestWithDouble, CheckVegasIntegrator)
{
	check_vegas_integrator<TypeParam>(GetParam());
}

INSTANTIATE_TEST_CASE_P(GenzTestSuite, GenzTestWithDouble,
	::testing::ValuesIn(test_data.begin(), test_data.end()));
