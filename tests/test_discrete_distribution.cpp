#include "gtest/gtest.h"

#include "hep/mc/internal/discrete_distribution.hpp"

#include <cstddef>
#include <random>
#include <vector>

typedef testing::Types<float, double, long double> MyT;
template <typename T> class DiscreteDistribution : public testing::Test { };
TYPED_TEST_CASE(DiscreteDistribution, MyT);

TYPED_TEST(DiscreteDistribution, DistributionWithOneBin)
{
	typedef TypeParam T;

	std::mt19937 rng;
	std::vector<T> weights = { T(1.0) };
	discrete_distribution<std::size_t, T> distribution(weights.begin(),
		weights.end());

	// distribution with one bin should always give that one bin
	EXPECT_EQ( distribution(rng) , 0u );
	EXPECT_EQ( distribution(rng) , 0u );
	EXPECT_EQ( distribution(rng) , 0u );
	EXPECT_EQ( distribution(rng) , 0u );
	EXPECT_EQ( distribution(rng) , 0u );
}

TYPED_TEST(DiscreteDistribution, DistributionWithTwoEqualBins)
{
	typedef TypeParam T;

	std::mt19937 rng;
	std::vector<T> weights = { T(1.0), T(1.0) };
	discrete_distribution<std::size_t, T> distribution(weights.begin(),
		weights.end());

	// distribution with one bin should always give that one bin
	std::size_t frequencies[2] = { 0, 0 };

	std::size_t const count = 10000;
	for (std::size_t i = 0; i != count; ++i)
	{
		std::size_t index = distribution(rng);

		if (index >= 2)
		{
			FAIL() << "index out of bounds";
		}

		++frequencies[index];
	}

	T const relative_frequencies[] = {
		T(frequencies[0]) / count,
		T(frequencies[1]) / count
	};

	EXPECT_NEAR(relative_frequencies[0], T(0.5), T(1e-1));
	EXPECT_NEAR(relative_frequencies[1], T(0.5), T(1e-1));
}

TYPED_TEST(DiscreteDistribution, DistributionWithTwoUnequalBins)
{
	typedef TypeParam T;

	std::mt19937 rng;
	std::vector<T> weights = { T(3.0), T(1.0) };
	discrete_distribution<std::size_t, T> distribution(weights.begin(),
		weights.end());

	// distribution with one bin should always give that one bin
	std::size_t frequencies[] = { 0, 0 };

	std::size_t const count = 10000;
	for (std::size_t i = 0; i != count; ++i)
	{
		std::size_t index = distribution(rng);

		if (index >= 2)
		{
			FAIL() << "index out of bounds";
		}

		++frequencies[index];
	}

	T const relative_frequencies[] = {
		T(frequencies[0]) / count,
		T(frequencies[1]) / count
	};

	EXPECT_NEAR(relative_frequencies[0], T(0.75), T(1e-1));
	EXPECT_NEAR(relative_frequencies[1], T(0.25), T(1e-1));
}

TYPED_TEST(DiscreteDistribution, DistributionWithOneHundredBins)
{
	typedef TypeParam T;

	std::mt19937 rng;
	std::vector<T> weights(100, T(1.0));
	discrete_distribution<std::size_t, T> distribution(weights.begin(),
		weights.end());

	// distribution with one bin should always give that one bin
	std::vector<std::size_t> frequencies(100, 0);

	std::size_t const count = 100000;
	for (std::size_t i = 0; i != count; ++i)
	{
		std::size_t index = distribution(rng);

		if (index >= frequencies.size())
		{
			FAIL() << "index out of bounds";
		}

		++frequencies[index];
	}

	for (std::size_t i = 0; i != frequencies.size(); ++i)
	{
		T const relative_frequency = T(frequencies[i]) / count;
		EXPECT_NEAR(relative_frequency, T(1.0) / count, T(1e-1));
	}
}
