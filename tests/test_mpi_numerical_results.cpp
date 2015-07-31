#include "gtest/gtest.h"

#include <vector>

#include "hep/mc-mpi.hpp"

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
		0xf.f7d32p-4F, 0x8.2b4edp-12F,
		0xf.f7d30p-4F, 0x8.2b505p-12F,
		0x8.02be3p-3F, 0x9.048eep-13F,
		0x8.00063p-3F, 0xc.e5fbep-14F,
		0xf.fe5c6p-4F, 0xb.46568p-14F,
		0xf.ffaddp-4F, 0xa.c42a2p-14F
	};
}

template <>
std::vector<double> reference_results() {
	return {
		0xf.ffc86404543c0p-4, 0x8.2b672f116b6c0p-12,
		0xf.ffc86404543c0p-4, 0x8.2b672f116b4c0p-12,
		0x8.01294905b4cd8p-3, 0x9.03ebf3d5a5800p-13,
		0x8.00331091e43d8p-3, 0xc.d9e34cb3fb030p-14,
		0xf.fe107b0b3d6a0p-4, 0xb.518ce429078a8p-14,
		0xf.ffd9e3eac8cf8p-4, 0xa.d9002081b21f0p-14
	};
}

template <>
std::vector<long double> reference_results() {
	return {
		0xf.ffc86404543bcdbp-4l, 0x8.2b672f116b780aep-12l,
		0xf.ffc86404543bcdbp-4l, 0x8.2b672f116b7817dp-12l,
		0x8.01294905b4cd55ap-3l, 0x9.03ebf3d5a522130p-13l,
		0x8.00331091e43d804p-3l, 0xc.d9e34cb3fb9e6f0p-14l,
		0xf.fe107b0b3d69c69p-4l, 0xb.518ce429071b33ap-14l,
		0xf.ffd9e3eac8cf7eap-4l, 0xa.d9002081b209eb4p-14l
	};
}

typedef testing::Types<float, double, long double> MyT;
template <typename T> class NumericalResults : public testing::Test { };
TYPED_TEST_CASE(NumericalResults, MyT);

TYPED_TEST(NumericalResults, CheckPlainIntegration)
{
	using T = TypeParam;

	auto const result = hep::mpi_plain<T>(
		MPI_COMM_WORLD,
		2,
		100000,
		function<T>
	);
	auto const reference = reference_results<T>();

	EXPECT_EQ( result.value , reference[0] );
	EXPECT_EQ( result.error , reference[1] );
}

TYPED_TEST(NumericalResults, CheckVegasIntegration)
{
	using T = TypeParam;

	auto const results = hep::mpi_vegas<T>(
		MPI_COMM_WORLD,
		2,
		std::vector<std::size_t>(5, 100000),
		function<T>
	);
	auto const reference = reference_results<T>();

	EXPECT_EQ( results[0].value , reference[2] );
	EXPECT_EQ( results[0].error , reference[3] );
	EXPECT_EQ( results[1].value , reference[4] );
	EXPECT_EQ( results[1].error , reference[5] );
	EXPECT_EQ( results[2].value , reference[6] );
	EXPECT_EQ( results[2].error , reference[7] );
	EXPECT_EQ( results[3].value , reference[8] );
	EXPECT_EQ( results[3].error , reference[9] );
	EXPECT_EQ( results[4].value , reference[10] );
	EXPECT_EQ( results[4].error , reference[11] );
}
