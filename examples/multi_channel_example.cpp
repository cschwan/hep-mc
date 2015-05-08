#include "hep/mc.hpp"

#include <cmath>
#include <cstddef>
#include <iostream>
#include <vector>

constexpr double mu = 0.0;
constexpr double s_min = -1.0;
constexpr double s_max = +1.0;

int main()
{
	std::vector<std::size_t> events = {
		1000, 2000, 2000, 2000, 2000,
		2000, 2000, 2000, 2000, 2000,
		2000, 2000, 2000, 2000, 2000
	};

	double const G_min = -1.0 / 3.0 * std::pow(s_min - mu, 3.0);
	double const G_max = -1.0 / 3.0 * std::pow(s_max - mu, 3.0);
	double const A = 1.0 / (G_max - G_min);

	auto const function = [=](hep::multi_channel_point<double> const& x) {
		double const r = x.point[0];
		double const s = mu + std::cbrt(-3.0 * (r / A + G_min));

		return std::exp(-s * s);
	};

	auto const densities = [=](
		std::size_t,
		std::vector<double> const& random_numbers,
		std::vector<double>& channel_densities
	) {
		double const r = random_numbers[0];
		double const s = mu + std::cbrt(-3.0 * (r / A + G_min));
		double const g = -A * (s - mu) * (s - mu);

		channel_densities[0] = g;
	};

	auto const results = hep::multi_channel<double>(
		1,
		events,
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
