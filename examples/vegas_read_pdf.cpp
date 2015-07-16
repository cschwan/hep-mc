#include "hep/mc.hpp"

#include <fstream>
#include <ios>
#include <iostream>

int main(int argc, char* argv[])
{
	if (argc > 1)
	{
		// open file
		std::ifstream input(argv[1]);

		// we need to know the dimensions and bins beforehand!
		hep::vegas_pdf<double> grid(1, 5);

		// read the pdf from the file
		input >> grid;

		// print the (left bin) boundaries of the pdf in scientific format
		std::cout << std::scientific << grid << "\n";
	}
	else
	{
		std::cout << "usage: " << argv[0] << " [p.d.f. file]\n\n";
		std::cout << "run `vegas_write_pdf > [p.d.f. file]` to generate a valid"
			" file\n";
	}

	return 0;
}
