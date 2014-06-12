#include <gtest/gtest.h>

#include <hep/mc/vegas_pdf.hpp>

#include <cstddef>
#include <sstream>
#include <vector>

typedef testing::Types<float, double, long double> MyT;
template <typename T> class VegasPdf : public testing::Test { };
TYPED_TEST_CASE(VegasPdf, MyT);

TYPED_TEST(VegasPdf, ConstructorAndMembers)
{
	typedef TypeParam T;

	hep::vegas_pdf<T> pdf(2, 4);

	EXPECT_EQ( 4U , pdf.bins() );
	EXPECT_EQ( 2U , pdf.dimensions() );

	EXPECT_EQ( T(0.25) , pdf(0, 0) );
	EXPECT_EQ( T(0.50) , pdf(0, 1) );
	EXPECT_EQ( T(0.75) , pdf(0, 2) );
	EXPECT_EQ( T(1.00) , pdf(0, 3) );
	EXPECT_EQ( T(0.25) , pdf(1, 0) );
	EXPECT_EQ( T(0.50) , pdf(1, 1) );
	EXPECT_EQ( T(0.75) , pdf(1, 2) );
	EXPECT_EQ( T(1.00) , pdf(1, 3) );
}

TYPED_TEST(VegasPdf, SetAndGetAndCopyConstructor)
{
	typedef TypeParam T;

	hep::vegas_pdf<T> pdf1(2, 3);

	pdf1(0, 0) = T(0.0);
	pdf1(0, 1) = T(0.1);
	pdf1(0, 2) = T(0.2);
	pdf1(1, 0) = T(0.3);
	pdf1(1, 1) = T(0.4);
	pdf1(1, 2) = T(0.5);

	EXPECT_EQ( T(0.0) , pdf1(0, 0) );
	EXPECT_EQ( T(0.1) , pdf1(0, 1) );
	EXPECT_EQ( T(0.2) , pdf1(0, 2) );
	EXPECT_EQ( T(0.3) , pdf1(1, 0) );
	EXPECT_EQ( T(0.4) , pdf1(1, 1) );
	EXPECT_EQ( T(0.5) , pdf1(1, 2) );

	hep::vegas_pdf<T> pdf2(pdf1);

	EXPECT_EQ( T(0.0) , pdf2(0, 0) );
	EXPECT_EQ( T(0.1) , pdf2(0, 1) );
	EXPECT_EQ( T(0.2) , pdf2(0, 2) );
	EXPECT_EQ( T(0.3) , pdf2(1, 0) );
	EXPECT_EQ( T(0.4) , pdf2(1, 1) );
	EXPECT_EQ( T(0.5) , pdf2(1, 2) );
}

TYPED_TEST(VegasPdf, InverseCumulativeDistributionFunction)
{
	typedef TypeParam T;

	std::vector<T> random_numbers(1);
	std::vector<std::size_t> bin(1);

	// check a uniform pdf
	hep::vegas_pdf<T> pdf(1, 5);

	for (std::size_t i = 0; i != 99; ++i)
	{
		SCOPED_TRACE(testing::Message("i=") << i);

		T const number = T(1+i) / T(100.0);
		std::size_t const bin_index = number * pdf.bins();
		random_numbers[0] = number;

		// weight for a uniform pdf is 1
		EXPECT_NEAR( T(1.0) , pdf.icdf(random_numbers, bin) , T(1e-6) );

		// the bins are equal-sized
		EXPECT_EQ( bin_index , bin[0] );

		// "random" number must be unmodified
		EXPECT_NEAR( number , random_numbers[0] , T(1e-7) );
	}

	// check a pdf in which bin 0 is twice as large as bin 1 ...
	pdf(0, 0) = T( 2.0) / T( 3.0);
	pdf(0, 1) = T( 8.0) / T( 9.0);
	pdf(0, 2) = T(26.0) / T(27.0);
	pdf(0, 3) = T(80.0) / T(81.0);
	pdf(0, 4) = T( 1.0);

	T const weight[] = {
		T(10.0) / T(  3.0),
		T(10.0) / T(  9.0),
		T(10.0) / T( 27.0),
		T(10.0) / T( 81.0),
		T(15.0) / T(243.0)
	};

	T const offset[] = {
		T(),
		T( 2.0) / T(5.0),
		T(10.0) / T(5.0),
		T(36.0) / T(5.0),
		T(76.0) / T(5.0)
	};

	for (std::size_t i = 0; i != 99; ++i)
	{
		SCOPED_TRACE(testing::Message("i=") << i);

		T const number = T(1+i) / T(100.0);
		std::size_t const bin_index = number * pdf.bins();
		random_numbers[0] = number;

		EXPECT_NEAR( weight[bin_index] , pdf.icdf(random_numbers, bin) , T(1e-6) );
		EXPECT_EQ( bin_index , bin[0] );
		EXPECT_FLOAT_EQ( (number + offset[bin_index]) * weight[bin_index], random_numbers[0] );
	}
}

TYPED_TEST(VegasPdf, StreamOperators)
{
	typedef TypeParam T;

	hep::vegas_pdf<T> pdf1(2, 3);

	pdf1(0, 0) = T(0.0);
	pdf1(0, 1) = T(0.1);
	pdf1(0, 2) = T(0.2);
	pdf1(1, 0) = T(1.0);
	pdf1(1, 1) = T(1.1);
	pdf1(1, 2) = T(1.2);

	std::stringstream iostream;

	iostream << pdf1;
	hep::vegas_pdf<T> pdf2(2, 3);
	iostream >> pdf2;

	EXPECT_NEAR( T(0.0) , pdf2(0, 0) , T(1e-10) );
	EXPECT_NEAR( T(0.1) , pdf2(0, 1) , T(1e-10) );
	EXPECT_NEAR( T(0.2) , pdf2(0, 2) , T(1e-10) );
	EXPECT_NEAR( T(1.0) , pdf2(1, 0) , T(1e-10) );
	EXPECT_NEAR( T(1.1) , pdf2(1, 1) , T(1e-10) );
	EXPECT_NEAR( T(1.2) , pdf2(1, 2) , T(1e-10) );
}
