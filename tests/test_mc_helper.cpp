#include <gtest/gtest.h>

#include <hep/mc/mc_helper.hpp>

#include <limits>

typedef testing::Types<float, double, long double> MyT;
template <typename T> class McHelper : public testing::Test { };
TYPED_TEST_CASE(McHelper, MyT);

TYPED_TEST(McHelper, CumulativeResultWithZero)
{
	typedef TypeParam T;

	std::vector<hep::mc_result<T>> zero_results;
	auto result = hep::cumulative_result<T>(zero_results.begin(), zero_results.end());

	// assert that results are reasonable even if there is no real input
	ASSERT_EQ( 0U  , result.calls );
	ASSERT_EQ( T() , result.value );
	ASSERT_EQ( T() , result.error );
}

TYPED_TEST(McHelper, CumulativeResultWithOne)
{
	typedef TypeParam T;

	std::vector<hep::mc_result<T>> one_result = {
		hep::mc_result<T>(1000, T(10.0), T(20.0))
	};
	auto result = hep::cumulative_result<T>(one_result.begin(), one_result.end());

	ASSERT_EQ( one_result.front().calls , result.calls );
	ASSERT_NEAR( one_result.front().value , result.value , T(1e-10) );
	ASSERT_NEAR( one_result.front().error , result.error , T(1e-10) );
}

TYPED_TEST(McHelper, CumulativeResultWithTwo)
{
	typedef TypeParam T;

	// two times the same result
	std::vector<hep::mc_result<T>> two_results = {
		hep::mc_result<T>(100, T(100.0), T(200.0)),
		hep::mc_result<T>(100, T(100.0), T(200.0))
	};
	auto result = hep::cumulative_result<T>(two_results.begin(), two_results.end());

	// should give two times the calls
	ASSERT_EQ( 2 * two_results.front().calls , result.calls );
	// the same result and
	ASSERT_NEAR( two_results.front().value , result.value , T(1e-10) );
	// an error reduces by 1/sqrt(2)
	ASSERT_NEAR( two_results.front().error / std::sqrt(T(2.0)) , result.error , T(1e-10) );
}

TYPED_TEST(McHelper, ChiSquareDofWithZero)
{
	typedef TypeParam T;

	std::vector<hep::mc_result<T>> zero_results;
	T result = hep::chi_square_dof<T>(zero_results.begin(), zero_results.end());

	ASSERT_NEAR( T() , result , T(1e-10) );
}

TYPED_TEST(McHelper, ChiSquareDofWithOne)
{
	typedef TypeParam T;

	std::vector<hep::mc_result<T>> one_result = {
		hep::mc_result<T>(1000, T(10.0), T(20.0))
	};
	T result = hep::chi_square_dof<T>(one_result.begin(), one_result.end());

	ASSERT_EQ( std::numeric_limits<T>::infinity() , result );
}

TYPED_TEST(McHelper, ChiSquareDofWithTwo)
{
	typedef TypeParam T;

	// two times the same result
	std::vector<hep::mc_result<T>> two_results = {
		hep::mc_result<T>(100, T(100.0), T(200.0)),
		hep::mc_result<T>(100, T(100.0), T(200.0))
	};
	T result = hep::chi_square_dof<T>(two_results.begin(), two_results.end());

	ASSERT_NEAR( T() , result , T(1e-10) );
}
