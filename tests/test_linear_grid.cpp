#include <gtest/gtest.h>

#include <hep/mc/linear_grid.hpp>

#include <sstream>

typedef testing::Types<float, double, long double> MyT;
template <typename T> class LinearGrid : public testing::Test { };
TYPED_TEST_CASE(LinearGrid, MyT);

TYPED_TEST(LinearGrid, ConstructorAndMembers)
{
	typedef TypeParam T;

	hep::linear_grid<T> grid(2, 3);

	ASSERT_EQ( 3U , grid.bins() );
	ASSERT_EQ( 2U , grid.dimensions() );

	ASSERT_EQ( T() , grid(0, 0) );
	ASSERT_EQ( T() , grid(0, 1) );
	ASSERT_EQ( T() , grid(0, 2) );
	ASSERT_EQ( T() , grid(1, 0) );
	ASSERT_EQ( T() , grid(1, 1) );
	ASSERT_EQ( T() , grid(1, 2) );
}

TYPED_TEST(LinearGrid, SetAndGetAndCopyConstructor)
{
	typedef TypeParam T;

	hep::linear_grid<T> grid1(2, 3);

	grid1(0, 0) = T(0.0);
	grid1(0, 1) = T(0.1);
	grid1(0, 2) = T(0.2);
	grid1(1, 0) = T(1.0);
	grid1(1, 1) = T(1.1);
	grid1(1, 2) = T(1.2);

	ASSERT_EQ( T(0.0) , grid1(0, 0) );
	ASSERT_EQ( T(0.1) , grid1(0, 1) );
	ASSERT_EQ( T(0.2) , grid1(0, 2) );
	ASSERT_EQ( T(1.0) , grid1(1, 0) );
	ASSERT_EQ( T(1.1) , grid1(1, 1) );
	ASSERT_EQ( T(1.2) , grid1(1, 2) );

	hep::linear_grid<T> grid2(grid1);

	ASSERT_EQ( T(0.0) , grid2(0, 0) );
	ASSERT_EQ( T(0.1) , grid2(0, 1) );
	ASSERT_EQ( T(0.2) , grid2(0, 2) );
	ASSERT_EQ( T(1.0) , grid2(1, 0) );
	ASSERT_EQ( T(1.1) , grid2(1, 1) );
	ASSERT_EQ( T(1.2) , grid2(1, 2) );
}

TYPED_TEST(LinearGrid, StreamOperators)
{
	typedef TypeParam T;

	hep::linear_grid<T> grid1(2, 3);

	grid1(0, 0) = T(0.0);
	grid1(0, 1) = T(0.1);
	grid1(0, 2) = T(0.2);
	grid1(1, 0) = T(1.0);
	grid1(1, 1) = T(1.1);
	grid1(1, 2) = T(1.2);

	std::stringstream iostream;

	iostream << grid1;
	hep::linear_grid<T> grid2(2, 3);
	iostream >> grid2;

	ASSERT_NEAR( T(0.0) , grid2(0, 0) , T(1e-10) );
	ASSERT_NEAR( T(0.1) , grid2(0, 1) , T(1e-10) );
	ASSERT_NEAR( T(0.2) , grid2(0, 2) , T(1e-10) );
	ASSERT_NEAR( T(1.0) , grid2(1, 0) , T(1e-10) );
	ASSERT_NEAR( T(1.1) , grid2(1, 1) , T(1e-10) );
	ASSERT_NEAR( T(1.2) , grid2(1, 2) , T(1e-10) );
}
