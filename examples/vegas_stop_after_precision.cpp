#include <hep/mc.hpp>

#include <cstddef>
#include <iostream>
#include <vector>

double square(hep::mc_point<double> const& x)
{
	return x.point[0] * x.point[0];
}

struct stop_after_precision
{
	stop_after_precision(double abs_error, double rel_error)
		: abs_error(abs_error)
		, rel_error(rel_error)
	{
	}

	bool operator()(std::vector<hep::vegas_iteration_result<double>> const& r)
	{
		hep::vegas_verbose_callback<double>(r);

		auto const result = hep::cumulative_result<double>(r.begin(), r.end());

		if (result.error < abs_error)
		{
			std::cout << "Absolute error " << result.error << " smaller than the limit ";
			std::cout << abs_error << "\n";
			return false;
		}

		if (result.error < rel_error * result.value)
		{
			std::cout << "Relative error smaller than the limit\n";
			return false;
		}

		return true;
	}

	double abs_error;
	double rel_error;
};

int main()
{
	// stop if error is better than 1% (=0.01) or lower than 0.0001
	hep::vegas_callback<double>(stop_after_precision(0.0001, 0.01));

	// print only 3 digits
	std::cout.precision(3);

	// perform 100 iterations with 1000 calls each _at maximum_
	auto results = hep::vegas<double>(
		1,
		std::vector<std::size_t>(100, 1000),
		square
	);

	return 0;
}
