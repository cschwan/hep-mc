#include <hep/mc.hpp>

#include <cstddef>
#include <iostream>
#include <vector>

// the function that will be integrated
double square(hep::mc_point<double> const& x)
{
	return x.point[0] * x.point[0];
}

int main()
{
	// this is what the integrator should give
	double reference_result = 1.0 / 3.0;

	// print reference result
	std::cout << "computing integral of x^2 from 0 to 1\n";
	std::cout << "reference result is " << reference_result << "\n\n";

	// set the verbose vegas callback function
	hep::vegas_callback<double>(hep::vegas_verbose_callback<double>);

	// perform 5 iteration with 1000 calls each; this function will also call
	// vegas_verbose_callback after each iteration which in turn prints the
	// individual iterations
	auto results = hep::vegas<double>(
		1,
		std::vector<std::size_t>(50, 1000),
		square
	);

	return 0;
}
