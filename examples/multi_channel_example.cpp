#include "hep/mc.hpp"

#include <cmath>
#include <cstddef>
#include <iostream>
#include <vector>

int main()
{
	constexpr double s0 = -10.0;
	constexpr double s1 = +10.0;

	auto const function = [](hep::multi_channel_point<double> const& x) {
		double s = x.coordinates[0];

		double const sms0 = s - s0;
		double const sms1 = s - s1;

		return (std::exp(-sms0 * sms0) + std::exp(-sms1 * sms1)) 
			/ std::sqrt(std::acos(-1.0));
	};

	auto const densities = [](
		std::size_t channel,
		std::vector<double> const& random_numbers,
		std::vector<double>& coordinates,
		std::vector<double>& channel_densities
	) {
		double s = std::tan(std::acos(-1.0) * (random_numbers[0] - 0.5));

		if (channel == 0)
		{
			s += s0;
		}
		else if (channel == 1)
		{
			s += s1;
		}

		coordinates[0] = s;

		double const sms0 = s - s0;
		double const sms1 = s - s1;

		double const g0 = 1.0 / (1.0 + sms0 * sms0) / std::acos(-1.0);
		double const g1 = 1.0 / (1.0 + sms1 * sms1) / std::acos(-1.0);

		// check the number of channels
		switch (channel_densities.size())
		{
		case 2:
			channel_densities[1] = g1;
		case 1:
			channel_densities[0] = g0;
		}

		return false;
	};

	auto const results = hep::multi_channel<double>(
		1,
		std::vector<std::size_t>(10, 1000),
		function,
		2,
		densities
	);

	for (auto const result : results)
	{
		std::cout << result.value() << " +- " << result.error() << "\t"
			<< result.channel_weights[0];

		if (result.channel_weights.size() == 2)
		{
			std::cout << " " << result.channel_weights[1];
		}

		std::cout << "\n";
	}

	return 0;
}
