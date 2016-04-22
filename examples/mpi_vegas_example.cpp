#include "hep/mc-mpi.hpp"

#include <mpi.h>

#include <cstddef>
#include <iostream>
#include <vector>

// the function that shall be integrated
double square(hep::mc_point<double> const& x)
{
	return 3.0 * x.point()[0] * x.point()[0];
}

// run this program in parallel by starting it with `mpiexec`, e.g.
//
//     mpiexec -n 8 mpi_vegas_example
//
// for the integration to be performed with eight processors. Notice that the
// approximation the integrator returns is independent of this number.
int main(int argc, char* argv[])
{
	// Initialize MPI
	MPI_Init(&argc, &argv);

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// print reference result, but only for processor with rak zero (to avoid
	// printing it more than once)
	if (rank == 0)
	{
		std::cout << ">> computing integral of 3*x^2 from 0 to 1 which is 1.0"
			"\n\n";
	}

	// set the verbose vegas callback function
	hep::mpi_vegas_callback<double>(hep::mpi_vegas_verbose_callback<double>);

	// perform 5 iteration with 1000 calls each; this function will also call
	// vegas_verbose_callback after each iteration which in turn prints the
	// individual iterations
	auto results = hep::mpi_vegas(
		MPI_COMM_WORLD,
		hep::make_integrand<double>(square, 1),
		std::vector<std::size_t>(5, 10000000)
	);

	// results contains the estimations for each iteration. We could take the
	// result from last iteration, but here we instead choose to combine the
	// results of all iterations but the first one in a cumulative result
	auto result = hep::cumulative_result0(results.begin() + 1, results.end());
	double chi_square_dof = hep::chi_square_dof0(results.begin() + 1,
		results.end());

	if (rank == 0)
	{
		// print the cumulative result
		std::cout << ">> cumulative result (excluding first iteration):\n>> N="
			<< result.calls() << " I=" << result.value() << " +- "
			<< result.error() << " chi^2/dof=" << chi_square_dof << "\n";
	}

	// clean up
	MPI_Finalize();

	return 0;
}
