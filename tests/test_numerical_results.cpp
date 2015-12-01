#include "gtest/gtest.h"

#include <algorithm>
#include <cstddef>
//#include <cstdio>
#include <vector>

#ifndef HEP_USE_MPI
#include "hep/mc.hpp"
#else
#include "hep/mc-mpi.hpp"
#endif

template <typename T>
T function(hep::mc_point<T> const& x)
{
	return T(3.0) / T(2.0) *
		(x.point[0] * x.point[0] + x.point[1] * x.point[1]);
}

template <typename T>
std::vector<T> reference_results();

template <>
std::vector<float> reference_results() {
	return {
		// PLAIN
		0xf.f7d32p-4f, 0x8.2b4edp-12f,
		// VEGAS
		0xf.f7d30p-4f, 0x8.2b4efp-12f,
		0x8.02be3p-3f, 0x9.04993p-13f,
		0x8.00063p-3f, 0xc.e633ep-14f,
		0xf.fe5c6p-4f, 0xb.4668cp-14f,
		0xf.ffadfp-4f, 0xa.c45a6p-14f,
		// MULTI CHANNEL
		0xf.f7d30p-4f, 0x8.2b4efp-12f,
		0x8.03f22p-3f, 0x8.31456p-12f,
		0x8.00e3dp-3f, 0x8.33374p-12f,
		0xf.fdc35p-4f, 0x8.34ac7p-12f,
		0xf.fbccep-4f, 0x8.30541p-12f
	};
}

template <>
std::vector<double> reference_results() {
	return {
		// PLAIN
		0xf.ffc86404543c0p-4, 0x8.2b672f116b6c0p-12,
		// VEGAS
		0xf.ffc86404543c0p-4, 0x8.2b672f116b6c8p-12,
		0x8.01294905b4ce0p-3, 0x9.03ebf3d5a5248p-13,
		0x8.00331091e43e0p-3, 0xc.d9e34cb3fbcf8p-14,
		0xf.fe107b0b3d6a0p-4, 0xb.518ce42907288p-14,
		0xf.ffd9e3eac8cf8p-4, 0xa.d9002081b1648p-14,
		// MULTI CHANNEL
		0xf.ffc86404543c0p-4, 0x8.2b672f116b6c8p-12,
		0x8.028e6efb5b418p-3, 0x8.32bfd3f0c9be8p-12,
		0xf.f928706e8fae8p-4, 0x8.2cded6f995ad8p-12,
		0x8.006378779d598p-3, 0x8.3454b2eb745a8p-12,
		0x8.010ec36645680p-3, 0x8.2e66988e3c488p-12
	};
}

template <>
std::vector<long double> reference_results() {
	return {
		// PLAIN
		0xf.ffc86404543bcdbp-4l, 0x8.2b672f116b780aep-12l,
		// VEGAS
		0xf.ffc86404543bcdbp-4l, 0x8.2b672f116b780adp-12l,
		0x8.01294905b4cd55ap-3l, 0x9.03ebf3d5a522088p-13l,
		0x8.00331091e43d804p-3l, 0xc.d9e34cb3fb9ebafp-14l,
		0xf.fe107b0b3d69c68p-4l, 0xb.518ce429071b2b8p-14l,
		0xf.ffd9e3eac8cf7eap-4l, 0xa.d9002081b209a36p-14l,
		// MULTI CHANNEL
		0xf.ffc86404543bcdbp-4l, 0x8.2b672f116b780adp-12l,
		0x8.028e6efb5b41702p-3l, 0x8.32bfd3f0c9aec74p-12l,
		0xf.f928706e8fae36ep-4l, 0x8.2cded6f995b63acp-12l,
		0x8.006378779d59461p-3l, 0x8.3454b2eb74a4954p-12l,
		0x8.010ec3664567e43p-3l, 0x8.2e66988e3c3f0f0p-12l
	};
}

typedef testing::Types<float, double, long double> MyT;
template <typename T> class NumericalResults : public testing::Test { };
TYPED_TEST_CASE(NumericalResults, MyT);

TYPED_TEST(NumericalResults, CheckPlainIntegration)
{
	using T = TypeParam;

	std::size_t const calls = 100000;

#ifndef HEP_USE_MPI
	auto const result = hep::plain<T>(
#else
	auto const result = hep::mpi_plain<T>(
		MPI_COMM_WORLD,
#endif
		2,
		calls,
		function<T>
	);
	auto const reference = reference_results<T>();

	EXPECT_EQ( result.value() , reference[0] );
	EXPECT_EQ( result.error() , reference[1] );

	EXPECT_EQ( result.calls() , calls );
}

TYPED_TEST(NumericalResults, CheckVegasIntegration)
{
	using T = TypeParam;

	std::size_t const calls = 100000;
	std::size_t const iterations = 5;

#ifndef HEP_USE_MPI
	auto const results = hep::vegas<T>(
#else
	auto const results = hep::mpi_vegas<T>(
		MPI_COMM_WORLD,
#endif
		2,
		std::vector<std::size_t>(iterations, calls),
		function<T>
	);
	auto const reference = reference_results<T>();

	ASSERT_EQ( results.size() , iterations );

	EXPECT_EQ( results[0].value() , reference[2] );
	EXPECT_EQ( results[0].error() , reference[3] );
	EXPECT_EQ( results[1].value() , reference[4] );
	EXPECT_EQ( results[1].error() , reference[5] );
	EXPECT_EQ( results[2].value() , reference[6] );
	EXPECT_EQ( results[2].error() , reference[7] );
	EXPECT_EQ( results[3].value() , reference[8] );
	EXPECT_EQ( results[3].error() , reference[9] );
	EXPECT_EQ( results[4].value() , reference[10] );
	EXPECT_EQ( results[4].error() , reference[11] );

	EXPECT_EQ( results[0].calls() , calls );
	EXPECT_EQ( results[1].calls() , calls );
	EXPECT_EQ( results[2].calls() , calls );
	EXPECT_EQ( results[3].calls() , calls );
	EXPECT_EQ( results[4].calls() , calls );

//	for (auto i : results)
//	{
//		std::printf("%La, %La,\n", static_cast <long double> (i.value()),
//			static_cast <long double> (i.error()));
//	}
}

TYPED_TEST(NumericalResults, CheckMultiChannelIntegration)
{
	using T = TypeParam;

	std::size_t const calls = 100000;
	std::size_t const iterations = 5;

	auto unit_densities = [](std::size_t, std::vector<T> const& random_numbers,
		std::vector<T>& coordinates, std::vector<T>& channel_densities) {

		std::copy(
			random_numbers.begin(),
			random_numbers.end(),
			coordinates.begin()
		);

		std::fill(channel_densities.begin(), channel_densities.end(), T(1.0));

		return T(1.0);
	};

#ifndef HEP_USE_MPI
	auto const results = hep::multi_channel<T>(
#else
	auto const results = hep::mpi_multi_channel<T>(
		MPI_COMM_WORLD,
#endif
		2,
		2,
		std::vector<std::size_t>(iterations, calls),
		function<T>,
		1,
		unit_densities
	);
	auto const reference = reference_results<T>();

	ASSERT_EQ( results.size() , iterations );

	EXPECT_EQ( results[0].value() , reference[12] );
	EXPECT_EQ( results[0].error() , reference[13] );
	EXPECT_EQ( results[1].value() , reference[14] );
	EXPECT_EQ( results[1].error() , reference[15] );
	EXPECT_EQ( results[2].value() , reference[16] );
	EXPECT_EQ( results[2].error() , reference[17] );
	EXPECT_EQ( results[3].value() , reference[18] );
	EXPECT_EQ( results[3].error() , reference[19] );
	EXPECT_EQ( results[4].value() , reference[20] );
	EXPECT_EQ( results[4].error() , reference[21] );

	EXPECT_EQ( results[0].calls() , calls );
	EXPECT_EQ( results[1].calls() , calls );
	EXPECT_EQ( results[2].calls() , calls );
	EXPECT_EQ( results[3].calls() , calls );
	EXPECT_EQ( results[4].calls() , calls );

//	for (auto i : results)
//	{
//		std::printf("%La, %La,\n", static_cast <long double> (i.value()),
//			static_cast <long double> (i.error()));
//	}
}
