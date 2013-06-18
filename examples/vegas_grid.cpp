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
		std::cout << "The grid for iteration " << i << ":\n";

		// grid for iteration i and dimension 0 (only have one dimension here)
		auto& grid = results[i].grid[0];

		// print the grid
		for (std::size_t j = 0; j != grid.size(); ++j)
		{
			std::cout << grid[j] << "\n";
		}
	}

	return 0;
}
