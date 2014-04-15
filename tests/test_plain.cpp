#define CATCH_CONFIG_MAIN
#include <catch.hpp>

#include <hep/mc/plain.hpp>

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

	T operator()(hep::mc_point<T> const& point) const
	{
		std::vector<T> const& x = point.point;

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

		default:
			return T(0.0);
		}
	}

private:
	std::size_t type;
};

typedef double T;

TEST_CASE("Test PLAIN with constant function", "[plain]")
{
	std::size_t steps = 10;
	hep::mc_result<T> result = hep::plain<T>(1, steps, integrand<T>(0));

	CHECK( result.value == T(1.0) );
	CHECK( result.error == T(0.0) );
}

TEST_CASE("Test PLAIN with sine function", "[plain]")
{
	std::size_t steps = 1000;
	hep::mc_result<T> result = hep::plain<T>(1, steps, integrand<T>(1));

	CHECK( result.value <= T(1.0) + result.error );
	CHECK( result.value >= T(1.0) - result.error );
}

TEST_CASE("Test PLAIN with exponential function", "[plain]")
{
	std::size_t steps = 1000;
	hep::mc_result<T> result = hep::plain<T>(1, steps, integrand<T>(2));
	T reference = std::exp(T(1.0)) - T(1.0);

	CHECK( result.value <= reference + result.error );
	CHECK( result.value >= reference - result.error );
}

TEST_CASE("Test PLAIN with square function", "[plain]")
{
	std::size_t steps = 10000;
	hep::mc_result<T> result = hep::plain<T>(1, steps, integrand<T>(3));
	T reference = T(1.0) / T(3.0);

	CHECK( result.value <= reference + result.error );
	CHECK( result.value >= reference - result.error );
}

TEST_CASE("Test PLAIN with 4D-square function", "[plain]")
{
	std::size_t steps = 10000;
	hep::mc_result<T> result = hep::plain<T>(4, steps, integrand<T>(4));
	T reference = T(4.0) / T(3.0);

	CHECK( result.value <= reference + T(2.0) * result.error );
	CHECK( result.value >= reference - T(2.0) * result.error );
}

TEST_CASE("Test PLAIN with reverse square-root function", "[plain]")
{
	std::size_t steps = 1000;
	hep::mc_result<T> result = hep::plain<T>(1, steps, integrand<T>(5));
	T reference = T(2.0);

	CHECK( result.value <= reference + T(2e-2) );
	CHECK( result.value >= reference - T(2e-2) );
}
