#include "hep/mc.hpp"

#include <cstddef>
#include <iostream>
#include <vector>

double square(hep::mc_point<double> const& x)
{
	return x.point[0] * x.point[0];
}

bool print_grid(std::vector<hep::vegas_iteration_result<double>> const& results)
{
	std::cout << "The grid for iteration " << (results.size()-1) << ":\n";

	// grid for iteration i
	auto& grid = results.back().pdf;

	// print the grid
	double previous = 0.0;
	for (std::size_t j = 1; j != grid.bins(); ++j)
	{
		// width of this bin
		double width = (grid.bin_left(0, j) - previous);
		// middle-point of the bin
		double x = previous + width / 2.0;
		// function value of the bin
		double height = 1.0 / (10.0 * width);

		std::cout << x << "\t" << height << "\n";
		previous = grid.bin_left(0, j);
	}

	return true;
}

int main()
{
	// set the callback function
	hep::vegas_callback<double>(print_grid);

	// print only 3 digits
	std::cout.precision(3);

	// perform 5 iteration with 1000 calls each
	auto results = hep::vegas<double>(
		1,
		std::vector<std::size_t>(5, 1000),
		square,
		3                                  // 3 bins for better illustration
	);

	std::cout << "\nStarting new try using adapted grid ...\n\n";

	// perform another 5 iterations with the grid from the last iteration
	auto new_results = hep::vegas<double>(
		std::vector<std::size_t>(5, 1000),
		square,
		results.back().pdf                 // grid from 'results' last iteration
	);

	return 0;
}
