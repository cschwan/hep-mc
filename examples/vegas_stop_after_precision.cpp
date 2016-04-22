#include "hep/mc.hpp"

#include <cstddef>
#include <iostream>
#include <vector>

double square(hep::mc_point<double> const& x)
{
	return 3.0 * x.point()[0] * x.point()[0];
}

struct stop_after_precision
{
	stop_after_precision(double abs_error, double rel_error)
		: abs_error(abs_error)
		, rel_error(rel_error)
	{
	}

	bool operator()(std::vector<hep::vegas_result<double>> const& r)
	{
		// print the results obtained so far
		hep::vegas_verbose_callback<double>(r);

		// compute cumulative result ...
		auto const result = hep::cumulative_result0(r.begin(), r.end());

		// ... and check for the absolute error ...
		if (result.error() < abs_error)
		{
			std::cout << ">> absolute error " << result.error()
				<< " is smaller than the limit " << abs_error << "\n";

			// returning false stops all remaining iterations
			return false;
		}

		// ... and the relative error
		if (result.error() < rel_error * result.value())
		{
			std::cout << ">> relative error "
				<< (result.error() / result.value())
				<< " is smaller than the limit " << rel_error << "\n";

			// returning false stops all remaining iterations
			return false;
		}

		return true;
	}

	double abs_error;
	double rel_error;
};

int main()
{
	// stop if error is better than lower than 0.0001 or better than 1% (=0.01)
	hep::vegas_callback<double>(stop_after_precision(0.0001, 0.01));

	// print only 3 digits
	std::cout.precision(3);

	// perform 100 iterations with 1000 calls each _at maximum_
	auto results = hep::vegas(
		hep::make_integrand<double>(square, 1),
		std::vector<std::size_t>(100, 1000)
	);

	return 0;
}
