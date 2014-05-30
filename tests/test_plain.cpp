#include <gtest/gtest.h>

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

typedef testing::Types<float, double, long double> MyT;
template <typename T> class Plain : public testing::Test { };
TYPED_TEST_CASE(Plain, MyT);

TYPED_TEST(Plain, IntegrateConstantFunction)
{
	typedef TypeParam T;

	hep::mc_result<T> result = hep::plain<T>(1, 10, integrand<T>(0));

	EXPECT_EQ( T(1.0) , result.value );
	EXPECT_EQ( T()    , result.error );
}

TYPED_TEST(Plain, IntegrateSineFunction)
{
	typedef TypeParam T;

	hep::mc_result<T> result = hep::plain<T>(1, 1000, integrand<T>(1));

	EXPECT_NEAR( T(1.0) , result.value , result.error );
}

TYPED_TEST(Plain, IntegrateExponentialFunction)
{
	typedef TypeParam T;

	hep::mc_result<T> result = hep::plain<T>(1, 1000, integrand<T>(2));

	EXPECT_NEAR( std::exp(T(1.0)) - T(1.0) , result.value , result.error );
}

TYPED_TEST(Plain, IntegrateSquareFunction)
{
	typedef TypeParam T;

	hep::mc_result<T> result = hep::plain<T>(1, 10000, integrand<T>(3));

	EXPECT_NEAR( T(1.0) / T(3.0) , result.value , result.error );
}

TYPED_TEST(Plain, Integrate4DSquareFunction)
{
	typedef TypeParam T;

	hep::mc_result<T> result = hep::plain<T>(4, 10000, integrand<T>(4));

	EXPECT_NEAR( T(4.0) / T(3.0), result.value , T(2.0) * result.error );
}

TYPED_TEST(Plain, IntegrateReverseSqureRootFunction)
{
	typedef TypeParam T;

	hep::mc_result<T> result = hep::plain<T>(1, 10000, integrand<T>(5));

	// not square-integrable therefore error is not a meaningful quantity
	EXPECT_NEAR( T(2.0) , result.value , T(1e-2) );
}
