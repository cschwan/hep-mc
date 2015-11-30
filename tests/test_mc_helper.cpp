#include "gtest/gtest.h"

#include "hep/mc/mc_helper.hpp"

#include <limits>
#include <vector>

typedef testing::Types<float, double, long double> MyT;
template <typename T> class McHelper : public testing::Test { };
TYPED_TEST_CASE(McHelper, MyT);

TYPED_TEST(McHelper, CumulativeResult0WithOne)
{
	typedef TypeParam T;

	std::vector<hep::mc_result<T>> one_result = {
		hep::mc_result<T>(1000, T(10.0), T(20.0))
	};
	auto result = hep::cumulative_result0(one_result.begin(), one_result.end());

	EXPECT_EQ( one_result.front().calls() , result.calls() );
	EXPECT_NEAR( one_result.front().value() , result.value() , T(1e-10) );
	EXPECT_NEAR( one_result.front().error() , result.error() , T(1e-10) );
}

TYPED_TEST(McHelper, CumulativeResult1WithOne)
{
	typedef TypeParam T;

	std::vector<hep::mc_result<T>> one_result = {
		hep::mc_result<T>(1000, T(10.0), T(20.0))
	};
	auto result = hep::cumulative_result1(one_result.begin(), one_result.end());

	EXPECT_EQ( one_result.front().calls() , result.calls() );
	EXPECT_NEAR( one_result.front().value() , result.value() , T(1e-10) );
	EXPECT_NEAR( one_result.front().error() , result.error() , T(1e-10) );
}

TYPED_TEST(McHelper, CumulativeResult0WithTwo)
{
	typedef TypeParam T;

	// two times the same result
	std::vector<hep::mc_result<T>> two_results = {
		hep::mc_result<T>(100, T(100.0), T(200.0)),
		hep::mc_result<T>(100, T(100.0), T(200.0))
	};
	auto result = hep::cumulative_result0(two_results.begin(),
		two_results.end());

	// should give two times the calls
	EXPECT_EQ( 2 * two_results.front().calls() , result.calls() );
	// the same result and
	EXPECT_NEAR( two_results.front().value() , result.value() , T(1e-10) );
	// an error reduces by 1/sqrt(2)
	EXPECT_NEAR( two_results.front().error() / std::sqrt(T(2.0)) ,
		result.error() , T(1e-9) );
}

TYPED_TEST(McHelper, CumulativeResult1WithTwo)
{
	typedef TypeParam T;

	// two times the same result
	std::vector<hep::mc_result<T>> two_results = {
		hep::mc_result<T>(100, T(100.0), T(200.0)),
		hep::mc_result<T>(100, T(100.0), T(200.0))
	};
	auto result = hep::cumulative_result1(two_results.begin(),
		two_results.end());

	// should give two times the calls
	EXPECT_EQ( 2 * two_results.front().calls() , result.calls() );
	// the same result and
	EXPECT_NEAR( two_results.front().value() , result.value() , T(1e-10) );
	// zero error
	EXPECT_EQ( T() , result.error() );
}

TYPED_TEST(McHelper, CumulativeResult1WithThree)
{
	typedef TypeParam T;

	// two times the same result
	std::vector<hep::mc_result<T>> two_results = {
		hep::mc_result<T>(100, T(100.0), T(200.0)),
		hep::mc_result<T>(100, T(100.0), T(200.0)),
		hep::mc_result<T>(100, T(100.0), T(200.0)),
	};
	auto result = hep::cumulative_result1(two_results.begin(),
		two_results.end());

	// should give three times the calls
	EXPECT_EQ( 3 * two_results.front().calls() , result.calls() );
	// the same result and
	EXPECT_NEAR( two_results.front().value() , result.value() , T(1e-10) );
	// zero error
	EXPECT_EQ( T() , result.error() );
}

TYPED_TEST(McHelper, ChiSquareDof0WithZero)
{
	typedef TypeParam T;

	std::vector<hep::mc_result<T>> zero_results;
	T result = hep::chi_square_dof0(zero_results.begin(), zero_results.end());

	EXPECT_NEAR( T() , result , T(1e-10) );
}

TYPED_TEST(McHelper, ChiSquareDof1WithZero)
{
	typedef TypeParam T;

	std::vector<hep::mc_result<T>> zero_results;
	T result = hep::chi_square_dof1(zero_results.begin(), zero_results.end());

	EXPECT_NEAR( T() , result , T(1e-10) );
}

TYPED_TEST(McHelper, ChiSquareDof0WithOne)
{
	typedef TypeParam T;

	std::vector<hep::mc_result<T>> one_result = {
		hep::mc_result<T>(1000, T(10.0), T(20.0))
	};
	T result = hep::chi_square_dof0(one_result.begin(), one_result.end());

	EXPECT_EQ( std::numeric_limits<T>::infinity() , result );
}

TYPED_TEST(McHelper, ChiSquareDof1WithOne)
{
	typedef TypeParam T;

	std::vector<hep::mc_result<T>> one_result = {
		hep::mc_result<T>(1000, T(10.0), T(20.0))
	};
	T result = hep::chi_square_dof1(one_result.begin(), one_result.end());

	EXPECT_EQ( std::numeric_limits<T>::infinity() , result );
}

TYPED_TEST(McHelper, ChiSquareDof0WithTwo)
{
	typedef TypeParam T;

	// two times the same result
	std::vector<hep::mc_result<T>> two_results = {
		hep::mc_result<T>(100, T(100.0), T(200.0)),
		hep::mc_result<T>(100, T(100.0), T(200.0))
	};
	T result = hep::chi_square_dof0(two_results.begin(), two_results.end());

	EXPECT_NEAR( T() , result , T(1e-10) );
}

TYPED_TEST(McHelper, ChiSquareDof1WithTwo)
{
	typedef TypeParam T;

	// two times the same result
	std::vector<hep::mc_result<T>> two_results = {
		hep::mc_result<T>(100, T(100.0), T(200.0)),
		hep::mc_result<T>(100, T(100.0), T(200.0))
	};
	T result = hep::chi_square_dof1(two_results.begin(), two_results.end());

	EXPECT_NEAR( T() , result , T(1e-10) );
}
