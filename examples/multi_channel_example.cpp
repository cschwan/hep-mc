#include "hep/mc.hpp"

#include <cmath>
#include <cstddef>
#include <iostream>
#include <vector>

int main()
{
	using T = double;

	std::vector<std::size_t> events = {
		1000, 2000, 2000, 2000, 2000,
		2000, 2000, 2000, 2000, 2000
	};

	auto const function = [](hep::mc_point<T> const& x) {
		return std::exp(T(-1.0) * x.point[0]);
	};

	auto const densities = [](
		std::vector<T> const& random_numbers,
		std::vector<T>& channel_densities
	) {
		T const x = T(3.0) * random_numbers[0];

		// three channels, each constant
		channel_densities[0] =                x < T(1.0) ? T(3.0) : T();
		channel_densities[1] = x >= T(1.0) && x < T(2.0) ? T(3.0) : T();
		channel_densities[2] = x >= T(2.0)               ? T(3.0) : T();
	};

	auto const results = hep::multi_channel<T>(
		1,
		events,
		function,
		3,
		densities
	);

	for (auto const result : results)
	{
		std::cout << result.value() << " +- " << result.error() << " weights: "
			<< result.channel_weights[0] << " " << result.channel_weights[1] 
			<< " " << result.channel_weights[2] << "\n";
	}

	return 0;
}
