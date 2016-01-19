#include "gtest/gtest.h"

#include "hep/mc/mc_point.hpp"

#include <vector>

typedef testing::Types<float, double, long double> MyT;
template <typename T> class McPoint : public testing::Test { };
TYPED_TEST_CASE(McPoint, MyT);

TYPED_TEST(McPoint, ConstructorAndMemberVariables)
{
	typedef TypeParam T;

	std::vector<T> point{ T(0.1), T(0.9), T(0.133), T(0.4) };
	hep::mc_point<T> result(point);

	EXPECT_EQ( T(1.0) , result.weight() );

	EXPECT_NEAR( T(0.1)   , result.point()[0] , T(1e-10) );
	EXPECT_NEAR( T(0.9)   , result.point()[1] , T(1e-10) );
	EXPECT_NEAR( T(0.133) , result.point()[2] , T(1e-10) );
	EXPECT_NEAR( T(0.4)   , result.point()[3] , T(1e-10) );
}
