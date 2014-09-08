#include "gtest/gtest.h"

#include "hep/mc/vegas_pdf.hpp"

#include <cstddef>
#include <random>
#include <vector>

TEST(VegasIcdfBug, CheckRandomFloatOnes)
{
	typedef float T;

	std::size_t const dimensions = 1;
	std::size_t const bins = 128;

	std::vector<T> random_numbers(dimensions, T(1.0));
	std::vector<std::size_t> bin_indices(dimensions);

	hep::vegas_pdf<T> pdf(dimensions, bins);

	T const weight = hep::vegas_icdf(pdf, random_numbers, bin_indices);

	EXPECT_LE( T(), weight );

	for (std::size_t i = 0; i != dimensions; ++i)
	{
		EXPECT_EQ( bins - 1, bin_indices[i] );
		EXPECT_GT( T(1.0), random_numbers[i] );
	}
}
