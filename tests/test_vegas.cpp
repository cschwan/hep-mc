#include "gtest/gtest.h"

#include "hep/mc/vegas.hpp"

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

	T operator()(hep::vegas_point<T> const& sample) const
	{
		std::vector<T> const& x = sample.point;

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

		case 6:
			return T(3.0) / T(2.0) * (T(1.0) -
				(T(2.0) * x[0] - T(1.0)) * (T(2.0) * x[0] - T(1.0)));

		case 7:
			return T(2.0) * x[0];

		default:
			return T(0.0);
		}
	}

private:
	std::size_t type;
};

typedef testing::Types<float, double, long double> MyT;
template <typename T> class Vegas : public testing::Test { };
TYPED_TEST_CASE(Vegas, MyT);

TYPED_TEST(Vegas, CompareWithCubaData)
{
	typedef TypeParam T;

	// uniform distribution
	hep::vegas_pdf<T> old_grid(1, 128);

	std::vector<T> adjustment_data(128 + 2);
	for (std::size_t i = 0; i != 128; ++i)
	{
		// gaussian data
		T const power = (T(0.5) - i / T(128.0)) / T(0.1);
		adjustment_data[i] = T(1000.0) * std::exp(-power * power);
	}
	auto new_grid = hep::vegas_adjust_grid(T(1.5), old_grid, adjustment_data);

	EXPECT_FLOAT_EQ( T(5.856790233298160e-02) , new_grid(0,   0) );
	EXPECT_FLOAT_EQ( T(1.023518853498444e-01) , new_grid(0,   1) );
	EXPECT_FLOAT_EQ( T(1.366676582172576e-01) , new_grid(0,   2) );
	EXPECT_FLOAT_EQ( T(1.646794900969992e-01) , new_grid(0,   3) );
	EXPECT_FLOAT_EQ( T(1.881023011333995e-01) , new_grid(0,   4) );
	EXPECT_FLOAT_EQ( T(2.080991704454243e-01) , new_grid(0,   5) );
	EXPECT_FLOAT_EQ( T(2.255537973056561e-01) , new_grid(0,   6) );
	EXPECT_FLOAT_EQ( T(2.409221800698078e-01) , new_grid(0,   7) );
	EXPECT_FLOAT_EQ( T(2.546166761713083e-01) , new_grid(0,   8) );
	EXPECT_FLOAT_EQ( T(2.670054526284730e-01) , new_grid(0,   9) );
	EXPECT_FLOAT_EQ( T(2.782234792224454e-01) , new_grid(0,  10) );
	EXPECT_FLOAT_EQ( T(2.885598777394001e-01) , new_grid(0,  11) );
	EXPECT_FLOAT_EQ( T(2.980273187882674e-01) , new_grid(0,  12) );
	EXPECT_FLOAT_EQ( T(3.068001682751398e-01) , new_grid(0,  13) );
	EXPECT_FLOAT_EQ( T(3.149850363756436e-01) , new_grid(0,  14) );
	EXPECT_FLOAT_EQ( T(3.226516721177158e-01) , new_grid(0,  15) );
	EXPECT_FLOAT_EQ( T(3.298643556560267e-01) , new_grid(0,  16) );
	EXPECT_FLOAT_EQ( T(3.366818450396390e-01) , new_grid(0,  17) );
	EXPECT_FLOAT_EQ( T(3.431179235844910e-01) , new_grid(0,  18) );
	EXPECT_FLOAT_EQ( T(3.491921594737729e-01) , new_grid(0,  19) );
	EXPECT_FLOAT_EQ( T(3.550009104495114e-01) , new_grid(0,  20) );
	EXPECT_FLOAT_EQ( T(3.605850397505242e-01) , new_grid(0,  21) );
	EXPECT_FLOAT_EQ( T(3.659030600258112e-01) , new_grid(0,  22) );
	EXPECT_FLOAT_EQ( T(3.709800511760474e-01) , new_grid(0,  23) );
	EXPECT_FLOAT_EQ( T(3.759230087274760e-01) , new_grid(0,  24) );
	EXPECT_FLOAT_EQ( T(3.806310694015796e-01) , new_grid(0,  25) );
	EXPECT_FLOAT_EQ( T(3.851954784636450e-01) , new_grid(0,  26) );
	EXPECT_FLOAT_EQ( T(3.896358622380245e-01) , new_grid(0,  27) );
	EXPECT_FLOAT_EQ( T(3.938863912817975e-01) , new_grid(0,  28) );
	EXPECT_FLOAT_EQ( T(3.980825073736518e-01) , new_grid(0,  29) );
	EXPECT_FLOAT_EQ( T(4.020754659556730e-01) , new_grid(0,  30) );
	EXPECT_FLOAT_EQ( T(4.060496489324251e-01) , new_grid(0,  31) );
	EXPECT_FLOAT_EQ( T(4.098332987385077e-01) , new_grid(0,  32) );
	EXPECT_FLOAT_EQ( T(4.136068332227977e-01) , new_grid(0,  33) );
	EXPECT_FLOAT_EQ( T(4.172217685058796e-01) , new_grid(0,  34) );
	EXPECT_FLOAT_EQ( T(4.208149222131841e-01) , new_grid(0,  35) );
	EXPECT_FLOAT_EQ( T(4.242945131735323e-01) , new_grid(0,  36) );
	EXPECT_FLOAT_EQ( T(4.277265787737108e-01) , new_grid(0,  37) );
	EXPECT_FLOAT_EQ( T(4.310974680364662e-01) , new_grid(0,  38) );
	EXPECT_FLOAT_EQ( T(4.343868139913527e-01) , new_grid(0,  39) );
	EXPECT_FLOAT_EQ( T(4.376694540540855e-01) , new_grid(0,  40) );
	EXPECT_FLOAT_EQ( T(4.408335842621441e-01) , new_grid(0,  41) );
	EXPECT_FLOAT_EQ( T(4.439977144702027e-01) , new_grid(0,  42) );
	EXPECT_FLOAT_EQ( T(4.470984245801854e-01) , new_grid(0,  43) );
	EXPECT_FLOAT_EQ( T(4.501540463795743e-01) , new_grid(0,  44) );
	EXPECT_FLOAT_EQ( T(4.532071044898807e-01) , new_grid(0,  45) );
	EXPECT_FLOAT_EQ( T(4.561702043484013e-01) , new_grid(0,  46) );
	EXPECT_FLOAT_EQ( T(4.591333042069219e-01) , new_grid(0,  47) );
	EXPECT_FLOAT_EQ( T(4.620662204631659e-01) , new_grid(0,  48) );
	EXPECT_FLOAT_EQ( T(4.649521465242148e-01) , new_grid(0,  49) );
	EXPECT_FLOAT_EQ( T(4.678380725852636e-01) , new_grid(0,  50) );
	EXPECT_FLOAT_EQ( T(4.706813332560154e-01) , new_grid(0,  51) );
	EXPECT_FLOAT_EQ( T(4.735048838124036e-01) , new_grid(0,  52) );
	EXPECT_FLOAT_EQ( T(4.763284343687919e-01) , new_grid(0,  53) );
	EXPECT_FLOAT_EQ( T(4.791079331671701e-01) , new_grid(0,  54) );
	EXPECT_FLOAT_EQ( T(4.818834500916072e-01) , new_grid(0,  55) );
	EXPECT_FLOAT_EQ( T(4.846554832347171e-01) , new_grid(0,  56) );
	EXPECT_FLOAT_EQ( T(4.873969493944923e-01) , new_grid(0,  57) );
	EXPECT_FLOAT_EQ( T(4.901384155542676e-01) , new_grid(0,  58) );
	EXPECT_FLOAT_EQ( T(4.928747480878182e-01) , new_grid(0,  59) );
	EXPECT_FLOAT_EQ( T(4.955958877971886e-01) , new_grid(0,  60) );
	EXPECT_FLOAT_EQ( T(4.983170275065589e-01) , new_grid(0,  61) );
	EXPECT_FLOAT_EQ( T(5.010355888679004e-01) , new_grid(0,  62) );
	EXPECT_FLOAT_EQ( T(5.037499704701970e-01) , new_grid(0,  63) );
	EXPECT_FLOAT_EQ( T(5.064643520724935e-01) , new_grid(0,  64) );
	EXPECT_FLOAT_EQ( T(5.091821352427399e-01) , new_grid(0,  65) );
	EXPECT_FLOAT_EQ( T(5.119032749521102e-01) , new_grid(0,  66) );
	EXPECT_FLOAT_EQ( T(5.146244146614806e-01) , new_grid(0,  67) );
	EXPECT_FLOAT_EQ( T(5.173584066191084e-01) , new_grid(0,  68) );
	EXPECT_FLOAT_EQ( T(5.200998727788836e-01) , new_grid(0,  69) );
	EXPECT_FLOAT_EQ( T(5.228413389386589e-01) , new_grid(0,  70) );
	EXPECT_FLOAT_EQ( T(5.256094511610799e-01) , new_grid(0,  71) );
	EXPECT_FLOAT_EQ( T(5.283849680855170e-01) , new_grid(0,  72) );
	EXPECT_FLOAT_EQ( T(5.311604850099543e-01) , new_grid(0,  73) );
	EXPECT_FLOAT_EQ( T(5.339824864026360e-01) , new_grid(0,  74) );
	EXPECT_FLOAT_EQ( T(5.368060369590243e-01) , new_grid(0,  75) );
	EXPECT_FLOAT_EQ( T(5.396421151359576e-01) , new_grid(0,  76) );
	EXPECT_FLOAT_EQ( T(5.425280411970065e-01) , new_grid(0,  77) );
	EXPECT_FLOAT_EQ( T(5.454139672580554e-01) , new_grid(0,  78) );
	EXPECT_FLOAT_EQ( T(5.483379970078687e-01) , new_grid(0,  79) );
	EXPECT_FLOAT_EQ( T(5.513010968663893e-01) , new_grid(0,  80) );
	EXPECT_FLOAT_EQ( T(5.542641967249100e-01) , new_grid(0,  81) );
	EXPECT_FLOAT_EQ( T(5.573066010012800e-01) , new_grid(0,  82) );
	EXPECT_FLOAT_EQ( T(5.603622228006688e-01) , new_grid(0,  83) );
	EXPECT_FLOAT_EQ( T(5.634504382466203e-01) , new_grid(0,  84) );
	EXPECT_FLOAT_EQ( T(5.666145684546789e-01) , new_grid(0,  85) );
	EXPECT_FLOAT_EQ( T(5.697786986627376e-01) , new_grid(0,  86) );
	EXPECT_FLOAT_EQ( T(5.730469202236491e-01) , new_grid(0,  87) );
	EXPECT_FLOAT_EQ( T(5.763362661785356e-01) , new_grid(0,  88) );
	EXPECT_FLOAT_EQ( T(5.796907213783420e-01) , new_grid(0,  89) );
	EXPECT_FLOAT_EQ( T(5.831227869785205e-01) , new_grid(0,  90) );
	EXPECT_FLOAT_EQ( T(5.865838287609512e-01) , new_grid(0,  91) );
	EXPECT_FLOAT_EQ( T(5.901769824682557e-01) , new_grid(0,  92) );
	EXPECT_FLOAT_EQ( T(5.937711470365723e-01) , new_grid(0,  93) );
	EXPECT_FLOAT_EQ( T(5.975446815208623e-01) , new_grid(0,  94) );
	EXPECT_FLOAT_EQ( T(6.013182160051523e-01) , new_grid(0,  95) );
	EXPECT_FLOAT_EQ( T(6.052794097762929e-01) , new_grid(0,  96) );
	EXPECT_FLOAT_EQ( T(6.092535927530451e-01) , new_grid(0,  97) );
	EXPECT_FLOAT_EQ( T(6.134429290140082e-01) , new_grid(0,  98) );
	EXPECT_FLOAT_EQ( T(6.176653308124899e-01) , new_grid(0,  99) );
	EXPECT_FLOAT_EQ( T(6.221057145868693e-01) , new_grid(0, 100) );
	EXPECT_FLOAT_EQ( T(6.266393008493750e-01) , new_grid(0, 101) );
	EXPECT_FLOAT_EQ( T(6.313473615234786e-01) , new_grid(0, 102) );
	EXPECT_FLOAT_EQ( T(6.362566749963183e-01) , new_grid(0, 103) );
	EXPECT_FLOAT_EQ( T(6.412970738224275e-01) , new_grid(0, 104) );
	EXPECT_FLOAT_EQ( T(6.466150940977146e-01) , new_grid(0, 105) );
	EXPECT_FLOAT_EQ( T(6.521595552067653e-01) , new_grid(0, 106) );
	EXPECT_FLOAT_EQ( T(6.579254342648549e-01) , new_grid(0, 107) );
	EXPECT_FLOAT_EQ( T(6.639602642056871e-01) , new_grid(0, 108) );
	EXPECT_FLOAT_EQ( T(6.703895452152064e-01) , new_grid(0, 109) );
	EXPECT_FLOAT_EQ( T(6.771573718222229e-01) , new_grid(0, 110) );
	EXPECT_FLOAT_EQ( T(6.843168056393930e-01) , new_grid(0, 111) );
	EXPECT_FLOAT_EQ( T(6.919264769944591e-01) , new_grid(0, 112) );
	EXPECT_FLOAT_EQ( T(7.000505380349344e-01) , new_grid(0, 113) );
	EXPECT_FLOAT_EQ( T(7.087586092426237e-01) , new_grid(0, 114) );
	EXPECT_FLOAT_EQ( T(7.181257254107776e-01) , new_grid(0, 115) );
	EXPECT_FLOAT_EQ( T(7.283429878696618e-01) , new_grid(0, 116) );
	EXPECT_FLOAT_EQ( T(7.394790466368668e-01) , new_grid(0, 117) );
	EXPECT_FLOAT_EQ( T(7.516898617754770e-01) , new_grid(0, 118) );
	EXPECT_FLOAT_EQ( T(7.652128431565691e-01) , new_grid(0, 119) );
	EXPECT_FLOAT_EQ( T(7.803899881190566e-01) , new_grid(0, 120) );
	EXPECT_FLOAT_EQ( T(7.975513221703457e-01) , new_grid(0, 121) );
	EXPECT_NEAR( T(8.172850153120408e-01) , new_grid(0, 122) , T(1e-6) );
	EXPECT_NEAR( T(8.402673799385958e-01) , new_grid(0, 123) , T(1e-6) );
	EXPECT_NEAR( T(8.676005324026509e-01) , new_grid(0, 124) , T(1e-6) );
	EXPECT_NEAR( T(9.011245100404384e-01) , new_grid(0, 125) , T(1e-6) );
	EXPECT_NEAR( T(9.435695207842165e-01) , new_grid(0, 126) , T(1e-6) );
	EXPECT_FLOAT_EQ( T(1.000000000000000e+00) , new_grid(0, 127) );
}

TYPED_TEST(Vegas, IntegrateSquareFunction)
{
	typedef TypeParam T;

	std::size_t const iterations = 10;
	std::size_t const steps_per_iteration = 10000;
	std::size_t const bins = 100;

	auto result = hep::vegas<T>(
		1,
		std::vector<std::size_t>(iterations, steps_per_iteration),
		integrand<T>(6),
		bins
	);

	// check if the result is within 3-sigma range
	EXPECT_NEAR( T(1.0) , result.back().value , T(3.0) * result.back().error );
}
