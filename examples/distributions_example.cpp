#include "hep/mc.hpp"

double arctan(hep::mc_point<double> const& point)
{
	double const projection = 10.0 * point.point()[0] - 5.0;
	return std::atan(projection);
}

void bin_projector_function(
	hep::mc_point<double> const& point,
	hep::bin_projector<double>& projector,
	hep::function_value<double> value
) {
	double const projection = 10.0 * point.point()[0] - 5.0;
	projector.add(0, projection, value.value());
}

int main()
{
	// a differential distribution with 100 bins, covering the interval [5,+5]
	hep::distribution_parameters<double> parameters(100, -5.0, 5.0);

	// create the projector (the function that actually projects the mc_point
	// onto the value of the x-axis, the boundary values, and the number of
	// bins)
	auto const projector = hep::make_distribution_projector<double>(
		bin_projector_function,
		parameters
	);

	// now integrate and record the differential distributions
	auto const result = hep::plain_distributions<double>(
		1,
		1000000,
		arctan,
		projector
	);

	// integral is zero
	std::cout << "integral is I = " << result.value() << " +- "
		<< result.error() << "\n\n";

	auto const& distribution = result.distributions()[0];

	// print the distribution - compare this with the plot of the function
	// 'atan(x)/10.0'
	for (std::size_t i = 0; i != distribution.mid_points().size(); ++i)
	{
		std::cout << distribution.mid_points()[i] << "\t"
			<< distribution.results()[i].value() << "\t"
			<< distribution.results()[i].error() << "\n";
	}
}
