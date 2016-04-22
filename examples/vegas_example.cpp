#include "hep/mc.hpp"

#include <cstddef>
#include <iostream>
#include <vector>

// the function that shall be integrated
double square(hep::mc_point<double> const& x)
{
	return 3.0 * x.point()[0] * x.point()[0];
}

int main()
{
	// print reference result
	std::cout << ">> computing integral of 3*x^2 from 0 to 1 which is 1.0\n\n";

	// set the verbose vegas callback function
	hep::vegas_callback<double>(hep::vegas_verbose_callback<double>);

	// perform 5 iteration with 1000 calls each; this function will also call
	// vegas_verbose_callback after each iteration which in turn prints the
	// individual iterations
	auto results = hep::vegas(
		hep::make_integrand<double>(square, 1),
		std::vector<std::size_t>(5, 1000)
	);

	// results contains the estimations for each iteration. We could take the
	// result from last iteration, but here we instead choose to combine the
	// results of all iterations but the first one in a cumulative result
	auto result = hep::cumulative_result0(results.begin() + 1, results.end());
	double chi_square_dof = hep::chi_square_dof0(results.begin() + 1,
		results.end());

	// print the cumulative result
	std::cout << ">> cumulative result (excluding first iteration):\n>> N="
		<< result.calls() << " I=" << result.value() << " +- " << result.error()
		<< " chi^2/dof=" << chi_square_dof << "\n";

	return 0;
}
