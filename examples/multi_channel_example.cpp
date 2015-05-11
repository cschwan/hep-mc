#include "hep/mc.hpp"

#include <cmath>
#include <cstddef>
#include <iostream>
#include <vector>

int main()
{
	auto const function = [=](hep::multi_channel_point<double> const& x) {
		double const r = x.point[0];
		double const s = std::tan(std::acos(-1.0) * (r - 0.5));

		return std::exp(-s * s) / std::sqrt(std::acos(-1.0));
	};

	auto const densities = [=](
		std::size_t,
		std::vector<double> const& random_numbers,
		std::vector<double>& channel_densities
	) {
		double const r = random_numbers[0];
		double const s = std::tan(std::acos(-1.0) * (r - 0.5));
		double const g = 1.0 / (1.0 + s * s) / std::acos(-1.0);

		channel_densities[0] = g;
	};

	auto const results = hep::multi_channel<double>(
		1,
		std::vector<std::size_t>(10, 100000),
		function,
		1,
		densities
	);

	for (auto const result : results)
	{
		std::cout << result.value() << " +- " << result.error() << "\n";
	}

	return 0;
}
