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
	// perform integration and gather result from each iteration
	auto results = hep::vegas<double>(
		1,                                 // dimension
		std::vector<std::size_t>(5, 1000), // 5 iterations each 1000 calls
		square                             // integrate a parabola
	);

	// this is what the algorithm should give
	double reference_result = 1.0 / 3.0;

	// print reference result
	std::cout << "integral of x^2 from 0 to 1 is " << reference_result;
	std::cout << "\n";

	// print result of each iteration
	std::cout << "results of each iteration:\n";
	for (std::size_t i = 0; i != results.size(); ++i)
	{
		auto& result = results[i];

		std::cout << i << " (N=" << result.calls << ") : I=" << result.value;
		std::cout << " +- " << result.error << "\n";
	}

	// combine iterations into single result
	auto result =
		hep::cumulative_result<double>(results.begin(), results.end());

	std::cout << "cumulative estimate (N=" << result.calls << ") : I=";
	std::cout << result.value << " +- " << result.error << "\n";

	// compute the approximated chi-square per degree of freedom
	double chi = hep::chi_square_dof<double>(results.begin(), results.end());

	std::cout << "chi-square/dof=" << chi << "\n";

	return 0;
}
