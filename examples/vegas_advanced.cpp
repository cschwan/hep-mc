#include <hep/mc.hpp>

#include <cstddef>
#include <iostream>
#include <vector>

double square(hep::mc_point<double> const& x)
{
	return x.point[0] * x.point[0];
}

int main()
{
	// perform 5 iteration with 1000 calls each
	auto results = hep::vegas<double>(
		1,
		std::vector<std::size_t>(5, 1000),
		square,
		10                                 // grid with 10 bins only
	);

	for (std::size_t i = 0; i != results.size(); ++i)
	{
		// combine results of iteration 0 to i
		auto result = hep::cumulative_result<double>(results.begin(),
			results.begin() + (i+1));

		// print combined result
		std::cout << i << " (N=" << result.calls << ") : I=";
		std::cout << result.value << " +- " << result.error << "\n";
	}

	return 0;
}
