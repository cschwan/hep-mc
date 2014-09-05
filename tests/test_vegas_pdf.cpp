#include "gtest/gtest.h"

#include "hep/mc/vegas_pdf.hpp"

#include <cstddef>
#include <sstream>
#include <vector>

typedef testing::Types<float, double, long double> MyT;
template <typename T> class VegasPdf : public testing::Test { };
TYPED_TEST_CASE(VegasPdf, MyT);

TYPED_TEST(VegasPdf, DefaultConstructorAndBinLeft)
{
	typedef TypeParam T;

	hep::vegas_pdf<T> pdf(2, 4);

	EXPECT_EQ( 4U , pdf.bins() );
	EXPECT_EQ( 2U , pdf.dimensions() );

	EXPECT_EQ( T(0.00) , pdf.bin_left(0, 0) );
	EXPECT_EQ( T(0.25) , pdf.bin_left(0, 1) );
	EXPECT_EQ( T(0.50) , pdf.bin_left(0, 2) );
	EXPECT_EQ( T(0.75) , pdf.bin_left(0, 3) );
	EXPECT_EQ( T(1.00) , pdf.bin_left(0, 4) );
	EXPECT_EQ( T(0.00) , pdf.bin_left(1, 0) );
	EXPECT_EQ( T(0.25) , pdf.bin_left(1, 1) );
	EXPECT_EQ( T(0.50) , pdf.bin_left(1, 2) );
	EXPECT_EQ( T(0.75) , pdf.bin_left(1, 3) );
	EXPECT_EQ( T(1.00) , pdf.bin_left(1, 4) );
}

TYPED_TEST(VegasPdf, SetBinLeft)
{
	typedef TypeParam T;

	hep::vegas_pdf<T> pdf(2, 3);

	pdf.set_bin_left(0, 0, T(0.0));
	pdf.set_bin_left(0, 1, T(0.1));
	pdf.set_bin_left(0, 2, T(0.2));
	pdf.set_bin_left(0, 3, T(1.0));
	pdf.set_bin_left(1, 0, T(0.0));
	pdf.set_bin_left(1, 1, T(0.4));
	pdf.set_bin_left(1, 2, T(0.5));
	pdf.set_bin_left(1, 3, T(1.0));

	EXPECT_EQ( T(0.0) , pdf.bin_left(0, 0) );
	EXPECT_EQ( T(0.1) , pdf.bin_left(0, 1) );
	EXPECT_EQ( T(0.2) , pdf.bin_left(0, 2) );
	EXPECT_EQ( T(1.0) , pdf.bin_left(0, 3) );
	EXPECT_EQ( T(0.0) , pdf.bin_left(1, 0) );
	EXPECT_EQ( T(0.4) , pdf.bin_left(1, 1) );
	EXPECT_EQ( T(0.5) , pdf.bin_left(1, 2) );
	EXPECT_EQ( T(1.0) , pdf.bin_left(1, 3) );
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
		EXPECT_NEAR( T(1.0) , vegas_icdf(pdf, random_numbers, bin) , T(1e-6) );

		// the bins are equal-sized
		EXPECT_EQ( bin_index , bin[0] );

		// "random" number must be unmodified
		EXPECT_NEAR( number , random_numbers[0] , T(1e-7) );
	}

	// check a pdf in which bin 0 is twice as large as bin 1 ...
	pdf.set_bin_left(0, 0, T(    )          );
	pdf.set_bin_left(0, 1, T( 2.0) / T( 3.0));
	pdf.set_bin_left(0, 2, T( 8.0) / T( 9.0));
	pdf.set_bin_left(0, 3, T(26.0) / T(27.0));
	pdf.set_bin_left(0, 4, T(80.0) / T(81.0));
	pdf.set_bin_left(0, 5, T( 1.0)          );

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

		EXPECT_NEAR( weight[bin_index] , hep::vegas_icdf(pdf, random_numbers, bin) , T(1e-6) );
		EXPECT_EQ( bin_index , bin[0] );
		EXPECT_FLOAT_EQ( (number + offset[bin_index]) * weight[bin_index], random_numbers[0] );
	}
}

TYPED_TEST(VegasPdf, StreamOperators)
{
	typedef TypeParam T;

	hep::vegas_pdf<T> pdf1(2, 3);

	pdf1.set_bin_left(0, 0, T(0.0));
	pdf1.set_bin_left(0, 1, T(0.1));
	pdf1.set_bin_left(0, 2, T(0.2));
	pdf1.set_bin_left(0, 3, T(1.0));
	pdf1.set_bin_left(1, 0, T(0.0));
	pdf1.set_bin_left(1, 1, T(0.4));
	pdf1.set_bin_left(1, 2, T(0.5));
	pdf1.set_bin_left(1, 3, T(1.0));

	std::stringstream iostream;

	iostream << pdf1;
	hep::vegas_pdf<T> pdf2(2, 3);
	iostream >> pdf2;

	EXPECT_FLOAT_EQ( T(0.0) , pdf2.bin_left(0, 0) );
	EXPECT_FLOAT_EQ( T(0.1) , pdf2.bin_left(0, 1) );
	EXPECT_FLOAT_EQ( T(0.2) , pdf2.bin_left(0, 2) );
	EXPECT_FLOAT_EQ( T(1.0) , pdf2.bin_left(0, 3) );
	EXPECT_FLOAT_EQ( T(0.0) , pdf2.bin_left(1, 0) );
	EXPECT_FLOAT_EQ( T(0.4) , pdf2.bin_left(1, 1) );
	EXPECT_FLOAT_EQ( T(0.5) , pdf2.bin_left(1, 2) );
	EXPECT_FLOAT_EQ( T(1.0) , pdf2.bin_left(1, 3) );
}

TYPED_TEST(VegasPdf, CompareWithCubaBins)
{
	typedef TypeParam T;

	// uniform distribution
	hep::vegas_pdf<T> old_pdf(1, 128);

	std::vector<T> adjustment_data(128 + 2);
	for (std::size_t i = 0; i != 128; ++i)
	{
		// gaussian data
		T const power = (T(0.5) - i / T(128.0)) / T(0.1);
		adjustment_data[i] = T(1000.0) * std::exp(-power * power);
	}
	auto new_pdf = hep::vegas_refine_pdf(old_pdf, T(1.5), adjustment_data);

	EXPECT_FLOAT_EQ( T(0.000000000000000e-00) , new_pdf.bin_left(0,   0) );
	EXPECT_FLOAT_EQ( T(5.856790233298160e-02) , new_pdf.bin_left(0,   1) );
	EXPECT_FLOAT_EQ( T(1.023518853498444e-01) , new_pdf.bin_left(0,   2) );
	EXPECT_FLOAT_EQ( T(1.366676582172576e-01) , new_pdf.bin_left(0,   3) );
	EXPECT_FLOAT_EQ( T(1.646794900969992e-01) , new_pdf.bin_left(0,   4) );
	EXPECT_FLOAT_EQ( T(1.881023011333995e-01) , new_pdf.bin_left(0,   5) );
	EXPECT_FLOAT_EQ( T(2.080991704454243e-01) , new_pdf.bin_left(0,   6) );
	EXPECT_FLOAT_EQ( T(2.255537973056561e-01) , new_pdf.bin_left(0,   7) );
	EXPECT_FLOAT_EQ( T(2.409221800698078e-01) , new_pdf.bin_left(0,   8) );
	EXPECT_FLOAT_EQ( T(2.546166761713083e-01) , new_pdf.bin_left(0,   9) );
	EXPECT_FLOAT_EQ( T(2.670054526284730e-01) , new_pdf.bin_left(0,  10) );
	EXPECT_FLOAT_EQ( T(2.782234792224454e-01) , new_pdf.bin_left(0,  11) );
	EXPECT_FLOAT_EQ( T(2.885598777394001e-01) , new_pdf.bin_left(0,  12) );
	EXPECT_FLOAT_EQ( T(2.980273187882674e-01) , new_pdf.bin_left(0,  13) );
	EXPECT_FLOAT_EQ( T(3.068001682751398e-01) , new_pdf.bin_left(0,  14) );
	EXPECT_FLOAT_EQ( T(3.149850363756436e-01) , new_pdf.bin_left(0,  15) );
	EXPECT_FLOAT_EQ( T(3.226516721177158e-01) , new_pdf.bin_left(0,  16) );
	EXPECT_FLOAT_EQ( T(3.298643556560267e-01) , new_pdf.bin_left(0,  17) );
	EXPECT_FLOAT_EQ( T(3.366818450396390e-01) , new_pdf.bin_left(0,  18) );
	EXPECT_FLOAT_EQ( T(3.431179235844910e-01) , new_pdf.bin_left(0,  19) );
	EXPECT_FLOAT_EQ( T(3.491921594737729e-01) , new_pdf.bin_left(0,  20) );
	EXPECT_FLOAT_EQ( T(3.550009104495114e-01) , new_pdf.bin_left(0,  21) );
	EXPECT_FLOAT_EQ( T(3.605850397505242e-01) , new_pdf.bin_left(0,  22) );
	EXPECT_FLOAT_EQ( T(3.659030600258112e-01) , new_pdf.bin_left(0,  23) );
	EXPECT_FLOAT_EQ( T(3.709800511760474e-01) , new_pdf.bin_left(0,  24) );
	EXPECT_FLOAT_EQ( T(3.759230087274760e-01) , new_pdf.bin_left(0,  25) );
	EXPECT_FLOAT_EQ( T(3.806310694015796e-01) , new_pdf.bin_left(0,  26) );
	EXPECT_FLOAT_EQ( T(3.851954784636450e-01) , new_pdf.bin_left(0,  27) );
	EXPECT_FLOAT_EQ( T(3.896358622380245e-01) , new_pdf.bin_left(0,  28) );
	EXPECT_FLOAT_EQ( T(3.938863912817975e-01) , new_pdf.bin_left(0,  29) );
	EXPECT_FLOAT_EQ( T(3.980825073736518e-01) , new_pdf.bin_left(0,  30) );
	EXPECT_FLOAT_EQ( T(4.020754659556730e-01) , new_pdf.bin_left(0,  31) );
	EXPECT_FLOAT_EQ( T(4.060496489324251e-01) , new_pdf.bin_left(0,  32) );
	EXPECT_FLOAT_EQ( T(4.098332987385077e-01) , new_pdf.bin_left(0,  33) );
	EXPECT_FLOAT_EQ( T(4.136068332227977e-01) , new_pdf.bin_left(0,  34) );
	EXPECT_FLOAT_EQ( T(4.172217685058796e-01) , new_pdf.bin_left(0,  35) );
	EXPECT_FLOAT_EQ( T(4.208149222131841e-01) , new_pdf.bin_left(0,  36) );
	EXPECT_FLOAT_EQ( T(4.242945131735323e-01) , new_pdf.bin_left(0,  37) );
	EXPECT_FLOAT_EQ( T(4.277265787737108e-01) , new_pdf.bin_left(0,  38) );
	EXPECT_FLOAT_EQ( T(4.310974680364662e-01) , new_pdf.bin_left(0,  39) );
	EXPECT_FLOAT_EQ( T(4.343868139913527e-01) , new_pdf.bin_left(0,  40) );
	EXPECT_FLOAT_EQ( T(4.376694540540855e-01) , new_pdf.bin_left(0,  41) );
	EXPECT_FLOAT_EQ( T(4.408335842621441e-01) , new_pdf.bin_left(0,  42) );
	EXPECT_FLOAT_EQ( T(4.439977144702027e-01) , new_pdf.bin_left(0,  43) );
	EXPECT_FLOAT_EQ( T(4.470984245801854e-01) , new_pdf.bin_left(0,  44) );
	EXPECT_FLOAT_EQ( T(4.501540463795743e-01) , new_pdf.bin_left(0,  45) );
	EXPECT_FLOAT_EQ( T(4.532071044898807e-01) , new_pdf.bin_left(0,  46) );
	EXPECT_FLOAT_EQ( T(4.561702043484013e-01) , new_pdf.bin_left(0,  47) );
	EXPECT_FLOAT_EQ( T(4.591333042069219e-01) , new_pdf.bin_left(0,  48) );
	EXPECT_FLOAT_EQ( T(4.620662204631659e-01) , new_pdf.bin_left(0,  49) );
	EXPECT_FLOAT_EQ( T(4.649521465242148e-01) , new_pdf.bin_left(0,  50) );
	EXPECT_FLOAT_EQ( T(4.678380725852636e-01) , new_pdf.bin_left(0,  51) );
	EXPECT_FLOAT_EQ( T(4.706813332560154e-01) , new_pdf.bin_left(0,  52) );
	EXPECT_FLOAT_EQ( T(4.735048838124036e-01) , new_pdf.bin_left(0,  53) );
	EXPECT_FLOAT_EQ( T(4.763284343687919e-01) , new_pdf.bin_left(0,  54) );
	EXPECT_FLOAT_EQ( T(4.791079331671701e-01) , new_pdf.bin_left(0,  55) );
	EXPECT_FLOAT_EQ( T(4.818834500916072e-01) , new_pdf.bin_left(0,  56) );
	EXPECT_FLOAT_EQ( T(4.846554832347171e-01) , new_pdf.bin_left(0,  57) );
	EXPECT_FLOAT_EQ( T(4.873969493944923e-01) , new_pdf.bin_left(0,  58) );
	EXPECT_FLOAT_EQ( T(4.901384155542676e-01) , new_pdf.bin_left(0,  59) );
	EXPECT_FLOAT_EQ( T(4.928747480878182e-01) , new_pdf.bin_left(0,  60) );
	EXPECT_FLOAT_EQ( T(4.955958877971886e-01) , new_pdf.bin_left(0,  61) );
	EXPECT_FLOAT_EQ( T(4.983170275065589e-01) , new_pdf.bin_left(0,  62) );
	EXPECT_FLOAT_EQ( T(5.010355888679004e-01) , new_pdf.bin_left(0,  63) );
	EXPECT_FLOAT_EQ( T(5.037499704701970e-01) , new_pdf.bin_left(0,  64) );
	EXPECT_FLOAT_EQ( T(5.064643520724935e-01) , new_pdf.bin_left(0,  65) );
	EXPECT_FLOAT_EQ( T(5.091821352427399e-01) , new_pdf.bin_left(0,  66) );
	EXPECT_FLOAT_EQ( T(5.119032749521102e-01) , new_pdf.bin_left(0,  67) );
	EXPECT_FLOAT_EQ( T(5.146244146614806e-01) , new_pdf.bin_left(0,  68) );
	EXPECT_FLOAT_EQ( T(5.173584066191084e-01) , new_pdf.bin_left(0,  69) );
	EXPECT_FLOAT_EQ( T(5.200998727788836e-01) , new_pdf.bin_left(0,  70) );
	EXPECT_FLOAT_EQ( T(5.228413389386589e-01) , new_pdf.bin_left(0,  71) );
	EXPECT_FLOAT_EQ( T(5.256094511610799e-01) , new_pdf.bin_left(0,  72) );
	EXPECT_FLOAT_EQ( T(5.283849680855170e-01) , new_pdf.bin_left(0,  73) );
	EXPECT_FLOAT_EQ( T(5.311604850099543e-01) , new_pdf.bin_left(0,  74) );
	EXPECT_FLOAT_EQ( T(5.339824864026360e-01) , new_pdf.bin_left(0,  75) );
	EXPECT_FLOAT_EQ( T(5.368060369590243e-01) , new_pdf.bin_left(0,  76) );
	EXPECT_FLOAT_EQ( T(5.396421151359576e-01) , new_pdf.bin_left(0,  77) );
	EXPECT_FLOAT_EQ( T(5.425280411970065e-01) , new_pdf.bin_left(0,  78) );
	EXPECT_FLOAT_EQ( T(5.454139672580554e-01) , new_pdf.bin_left(0,  79) );
	EXPECT_FLOAT_EQ( T(5.483379970078687e-01) , new_pdf.bin_left(0,  80) );
	EXPECT_FLOAT_EQ( T(5.513010968663893e-01) , new_pdf.bin_left(0,  81) );
	EXPECT_FLOAT_EQ( T(5.542641967249100e-01) , new_pdf.bin_left(0,  82) );
	EXPECT_FLOAT_EQ( T(5.573066010012800e-01) , new_pdf.bin_left(0,  83) );
	EXPECT_FLOAT_EQ( T(5.603622228006688e-01) , new_pdf.bin_left(0,  84) );
	EXPECT_FLOAT_EQ( T(5.634504382466203e-01) , new_pdf.bin_left(0,  85) );
	EXPECT_FLOAT_EQ( T(5.666145684546789e-01) , new_pdf.bin_left(0,  86) );
	EXPECT_FLOAT_EQ( T(5.697786986627376e-01) , new_pdf.bin_left(0,  87) );
	EXPECT_FLOAT_EQ( T(5.730469202236491e-01) , new_pdf.bin_left(0,  88) );
	EXPECT_FLOAT_EQ( T(5.763362661785356e-01) , new_pdf.bin_left(0,  89) );
	EXPECT_FLOAT_EQ( T(5.796907213783420e-01) , new_pdf.bin_left(0,  90) );
	EXPECT_FLOAT_EQ( T(5.831227869785205e-01) , new_pdf.bin_left(0,  91) );
	EXPECT_FLOAT_EQ( T(5.865838287609512e-01) , new_pdf.bin_left(0,  92) );
	EXPECT_FLOAT_EQ( T(5.901769824682557e-01) , new_pdf.bin_left(0,  93) );
	EXPECT_FLOAT_EQ( T(5.937711470365723e-01) , new_pdf.bin_left(0,  94) );
	EXPECT_FLOAT_EQ( T(5.975446815208623e-01) , new_pdf.bin_left(0,  95) );
	EXPECT_FLOAT_EQ( T(6.013182160051523e-01) , new_pdf.bin_left(0,  96) );
	EXPECT_FLOAT_EQ( T(6.052794097762929e-01) , new_pdf.bin_left(0,  97) );
	EXPECT_FLOAT_EQ( T(6.092535927530451e-01) , new_pdf.bin_left(0,  98) );
	EXPECT_FLOAT_EQ( T(6.134429290140082e-01) , new_pdf.bin_left(0,  99) );
	EXPECT_FLOAT_EQ( T(6.176653308124899e-01) , new_pdf.bin_left(0, 100) );
	EXPECT_FLOAT_EQ( T(6.221057145868693e-01) , new_pdf.bin_left(0, 101) );
	EXPECT_FLOAT_EQ( T(6.266393008493750e-01) , new_pdf.bin_left(0, 102) );
	EXPECT_FLOAT_EQ( T(6.313473615234786e-01) , new_pdf.bin_left(0, 103) );
	EXPECT_FLOAT_EQ( T(6.362566749963183e-01) , new_pdf.bin_left(0, 104) );
	EXPECT_FLOAT_EQ( T(6.412970738224275e-01) , new_pdf.bin_left(0, 105) );
	EXPECT_FLOAT_EQ( T(6.466150940977146e-01) , new_pdf.bin_left(0, 106) );
	EXPECT_FLOAT_EQ( T(6.521595552067653e-01) , new_pdf.bin_left(0, 107) );
	EXPECT_FLOAT_EQ( T(6.579254342648549e-01) , new_pdf.bin_left(0, 108) );
	EXPECT_FLOAT_EQ( T(6.639602642056871e-01) , new_pdf.bin_left(0, 109) );
	EXPECT_FLOAT_EQ( T(6.703895452152064e-01) , new_pdf.bin_left(0, 110) );
	EXPECT_FLOAT_EQ( T(6.771573718222229e-01) , new_pdf.bin_left(0, 111) );
	EXPECT_FLOAT_EQ( T(6.843168056393930e-01) , new_pdf.bin_left(0, 112) );
	EXPECT_FLOAT_EQ( T(6.919264769944591e-01) , new_pdf.bin_left(0, 113) );
	EXPECT_FLOAT_EQ( T(7.000505380349344e-01) , new_pdf.bin_left(0, 114) );
	EXPECT_FLOAT_EQ( T(7.087586092426237e-01) , new_pdf.bin_left(0, 115) );
	EXPECT_FLOAT_EQ( T(7.181257254107776e-01) , new_pdf.bin_left(0, 116) );
	EXPECT_FLOAT_EQ( T(7.283429878696618e-01) , new_pdf.bin_left(0, 117) );
	EXPECT_FLOAT_EQ( T(7.394790466368668e-01) , new_pdf.bin_left(0, 118) );
	EXPECT_FLOAT_EQ( T(7.516898617754770e-01) , new_pdf.bin_left(0, 119) );
	EXPECT_FLOAT_EQ( T(7.652128431565691e-01) , new_pdf.bin_left(0, 120) );
	EXPECT_FLOAT_EQ( T(7.803899881190566e-01) , new_pdf.bin_left(0, 121) );
	EXPECT_FLOAT_EQ( T(7.975513221703457e-01) , new_pdf.bin_left(0, 122) );
	EXPECT_NEAR(     T(8.172850153120408e-01) , new_pdf.bin_left(0, 123) , T(1e-6) );
	EXPECT_NEAR(     T(8.402673799385958e-01) , new_pdf.bin_left(0, 124) , T(1e-6) );
	EXPECT_NEAR(     T(8.676005324026509e-01) , new_pdf.bin_left(0, 125) , T(1e-6) );
	EXPECT_NEAR(     T(9.011245100404384e-01) , new_pdf.bin_left(0, 126) , T(1e-6) );
	EXPECT_NEAR(     T(9.435695207842165e-01) , new_pdf.bin_left(0, 127) , T(1e-6) );
	EXPECT_FLOAT_EQ( T(1.000000000000000e+00) , new_pdf.bin_left(0, 128) );
}

TYPED_TEST(VegasPdf, RefineWithAlphaZero)
{
	typedef TypeParam T;

	// uniform distribution
	hep::vegas_pdf<T> old_pdf(1, 128);

	std::vector<T> adjustment_data(128 + 2);
	for (std::size_t i = 0; i != 128; ++i)
	{
		// gaussian data
		T const power = (T(0.5) - i / T(128.0)) / T(0.1);
		adjustment_data[i] = T(1000.0) * std::exp(-power * power);
	}
	auto new_pdf = hep::vegas_refine_pdf(old_pdf, T(), adjustment_data);

	for (std::size_t i = 0; i != old_pdf.bins() + 1; ++i)
	{
		EXPECT_EQ( old_pdf.bin_left(0, i) , new_pdf.bin_left(0, i) );
	}
}
