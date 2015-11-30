#include "gtest/gtest.h"

#include "hep/mc/mc_result.hpp"

typedef testing::Types<float, double, long double> MyT;
template <typename T> class McResult : public testing::Test { };
TYPED_TEST_CASE(McResult, MyT);

TYPED_TEST(McResult, ConstructorAndMemberVariables)
{
	typedef TypeParam T;

	hep::mc_result<T> result1(100, T(100.0), T(100.0));

	EXPECT_EQ( 100U , result1.calls() );
	EXPECT_NEAR( T(1.0) , result1.value() , T(1e-10) );
	EXPECT_NEAR( T()    , result1.error() , T(1e-10) );
}

TYPED_TEST(McResult, CreateResult)
{
	typedef TypeParam T;

	auto const result = hep::create_result(100, T(1.0), T());

	EXPECT_EQ( 100U , result.calls() );
	EXPECT_NEAR( T(1.0) , result.value() , T(1e-10) );
	EXPECT_NEAR( T()    , result.error() , T(1e-10) );
}
