#include <hep/mc.hpp>
#include <hep/mc-mpi.hpp>

#include <cmath>
#include <cstddef>
#include <iostream>

#include <mpi.h>

double gauss(hep::mc_point<double> const& sample)
{
	double result = 1.0;

	// acos(-1.0) = pi - there is no pi constant in C++!!
	double factor = 2.0 / std::sqrt(std::acos(-1.0));

	// multiply a (half-)gaussian for every dimension
	for (std::size_t i = 0; i != sample.point.size(); ++i)
	{
		double argument = sample.point[i];
		result *= std::exp(-argument * argument);
		result *= factor;
	}

	return result;
}

int main(int argc, char* argv[])
{
	// initialize MPI
	MPI_Init(&argc, &argv);

	// get id for this process
	int rank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// set the verbose vegas callback function
	hep::mpi_vegas_callback<double>(hep::mpi_vegas_verbose_callback<double>);

	// dimension of the gauss
	std::size_t const dimension = 10;

	// five iterations with 10^7 calls each
	std::vector<std::size_t> iterations(5, 10000000);

	// perform the integration
	auto results = hep::mpi_vegas<double>(
		MPI_COMM_WORLD,
		dimension,
		iterations,
		gauss
	);

	// clean up
	MPI_Finalize();

	return 0;
}
