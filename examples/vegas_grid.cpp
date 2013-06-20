#include <hep/mc.hpp>

#include <cstddef>
#include <iostream>
#include <vector>

double square(hep::mc_point<double> const& x)
{
	return x.point[0] * x.point[0];
}

void print_grid(std::vector<hep::vegas_iteration_result<double>> const& results)
{
	std::cout << "The grid for iteration " << (results.size()-1) << ":\n";

	// grid for iteration i and dimension 0 (only have one dimension here)
	auto& grid = results.back().grid[0];

	// print the grid
	double previous = 0.0;
	for (std::size_t j = 0; j != grid.size(); ++j)
	{
		// width of this bin
		double width = (grid[j] - previous);
		// middle-point of the bin
		double x = previous + width / 2.0;
		// function value of the bin
		double height = 1.0 / (10.0 * width);

		std::cout << x << "\t" << height << "\n";
		previous = grid[j];
	}
}

int main()
{
	// set the callback function
	hep::vegas_callback<double>(print_grid);

	// perform 5 iteration with 1000 calls each
	auto results = hep::vegas<double>(
		1,
		std::vector<std::size_t>(5, 1000),
		square,
		10                                 // grid with 10 bins only
	);

	return 0;
}
