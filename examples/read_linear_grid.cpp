#include "hep/mc.hpp"

#include <fstream>
#include <iostream>

int main(int argc, char* argv[])
{
	if (argc > 1)
	{
		// open file
		std::ifstream input(argv[1]);

		// we need to know the dimensions and bins beforehand!
		hep::vegas_pdf<double> grid(2, 5);

		// read grid itself
		input >> grid;

		// print the grid in scientific format
		std::cout << std::scientific << grid << "\n";
	}
	else
	{
		std::cout << "usage: " << argv[0] << " [grid file]\n";
	}

	return 0;
}
