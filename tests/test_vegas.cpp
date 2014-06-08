#include <gtest/gtest.h>

#include <hep/mc/vegas.hpp>

#include <cmath>
#include <cstddef>
#include <vector>

template <typename T>
class integrand
{
public:
	integrand(std::size_t type)
		: type(type)
	{
	}

	T operator()(hep::vegas_point<T> const& sample) const
	{
		std::vector<T> const& x = sample.point;

		switch (type)
		{
		case 0:
			return T(1.0);

		case 1:
			return std::acos(T(-1.0)) / T(2.0) *
				std::sin(std::acos(T(-1.0)) * x[0]);

		case 2:
			return std::exp(x[0]);

		case 3:
			return x[0] * x[0];

		case 4:
			return x[0] * x[0] + x[1] * x[1] + x[2] * x[2] + x[3] * x[3];

		case 5: // not square-integrable
			return T(1.0) / std::sqrt(x[0]);

		case 6:
			return T(3.0) / T(2.0) * (T(1.0) -
				(T(2.0) * x[0] - T(1.0)) * (T(2.0) * x[0] - T(1.0)));

		case 7:
			return T(2.0) * x[0];

		default:
			return T(0.0);
		}
	}

private:
	std::size_t type;
};

typedef testing::Types<float, double, long double> MyT;
template <typename T> class Vegas : public testing::Test { };
TYPED_TEST_CASE(Vegas, MyT);

TYPED_TEST(Vegas, CheckAdjustGrid)
{
	typedef TypeParam T;

	// uniform distribution
	hep::piecewise_constant_pdf<T> old_grid(1, 4);

	// with this data it should remain a uniform distribution
	std::vector<T> adjustment_data = { T(1.0), T(1.0), T(1.0), T(1.0), T(), T() };
	auto new_grid = hep::vegas_adjust_grid(T(1.5), old_grid, adjustment_data);

	EXPECT_EQ( old_grid(0, 0), new_grid(0, 0) );
	EXPECT_EQ( old_grid(0, 1), new_grid(0, 1) );
	EXPECT_EQ( old_grid(0, 2), new_grid(0, 2) );
	EXPECT_EQ( old_grid(0, 3), new_grid(0, 3) );

	// this should make the first bin(s) smaller
	adjustment_data = { T(4.0), T(1.0), T(1.0), T(1.0), T(), T() };
	new_grid = hep::vegas_adjust_grid(T(1.0), old_grid, adjustment_data);

	EXPECT_FLOAT_EQ( T(2.072390895598839e-01), new_grid(0, 0) );
	EXPECT_FLOAT_EQ( T(4.303460252668945e-01), new_grid(0, 1) );
	EXPECT_FLOAT_EQ( T(7.047478748905420e-01), new_grid(0, 2) );
	EXPECT_FLOAT_EQ( T(1.000000000000000e+00), new_grid(0, 3) );

	new_grid = hep::vegas_adjust_grid(T(1.5), old_grid, adjustment_data);

	EXPECT_FLOAT_EQ( T(1.904424247771927e-01), new_grid(0, 0) );
	EXPECT_FLOAT_EQ( T(4.002750242439110e-01), new_grid(0, 1) );
	EXPECT_FLOAT_EQ( T(6.761486310730158e-01), new_grid(0, 2) );
	EXPECT_FLOAT_EQ( T(1.000000000000000e+00), new_grid(0, 3) );

	new_grid = hep::vegas_adjust_grid(T(2.0), old_grid, adjustment_data);

	EXPECT_FLOAT_EQ( T(1.760695619039624e-01), new_grid(0, 0) );
	EXPECT_FLOAT_EQ( T(3.727972665008580e-01), new_grid(0, 1) );
	EXPECT_FLOAT_EQ( T(6.426226726635829e-01), new_grid(0, 2) );
	EXPECT_FLOAT_EQ( T(1.000000000000000e+00), new_grid(0, 3) );
}

TYPED_TEST(Vegas, IntegrateSquareFunction)
{
	typedef TypeParam T;

	std::size_t const iterations = 10;
	std::size_t const steps_per_iteration = 10000;
	std::size_t const bins = 100;

	auto result = hep::vegas<T>(
		1,
		std::vector<std::size_t>(iterations, steps_per_iteration),
		integrand<T>(6),
		bins
	);

	// check if the result is within 3-sigma range
	EXPECT_NEAR( T(1.0) , result.back().value , T(3.0) * result.back().error );
}
