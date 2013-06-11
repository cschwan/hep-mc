#include <hep/mc.hpp>

#include <cmath>
#include <iostream>
#include <vector>

struct monomial
{
	double exponent;

	// constructor. sets the exponent
	monomial(double exponent)
		: exponent(exponent)
	{
	}

	// integrand function. raises the argument to the power 'exponent'
	double operator()(hep::mc_point<double> const& x)
	{
		return std::pow(x.point[0], exponent);
	}
};

int main()
{
	double exponent = 2.0;

	std::vector<hep::vegas_iteration_result<double>> results =
		hep::vegas<double>(
			1,                                 // integration in 1 dimension
			std::vector<std::size_t>(5, 1000), // 5 iterations, 1000 calls each
			monomial(exponent)                 // integrate a parabola
		);

	double reference_result = (1.0 / (exponent + 1.0));

	// print numbers in scientific format
	std::cout.setf(std::ios::scientific);

	// print reference result
	std::cout << "integral of x^" << exponent << " is " << reference_result;
	std::cout << "\n";

	// print result of each iteration
	std::cout << "results of each iteration:\n";
	for (std::size_t i = 0; i != results.size(); ++i)
	{
		hep::vegas_iteration_result<double>& result = results[i];

		std::cout << i << " (N=" << result.calls << ") : I=" << result.value;
		std::cout << " +- " << result.error << "\n";
	}

	return 0;
}
