#include "gtest/gtest.h"

#include "hep/mc/mpi_vegas.hpp"

#include <cmath>
#include <cstddef>
#include <vector>

#include <mpi.h>

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

TYPED_TEST(Vegas, IntegrateSquareFunction)
{
	typedef TypeParam T;

	std::size_t const iterations = 10;
	std::size_t const steps_per_iteration = 10000;
	std::size_t const bins = 100;

	auto result = hep::mpi_vegas<T>(
		MPI_COMM_WORLD,
		1,
		std::vector<std::size_t>(iterations, steps_per_iteration),
		integrand<T>(6),
		bins
	);

	// check if the result is within 3-sigma range
	EXPECT_NEAR( T(1.0) , result.back().value , T(3.0) * result.back().error );
}
