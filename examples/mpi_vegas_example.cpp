#include "hep/mc.hpp"
#include "hep/mc-mpi.hpp"

#include <cmath>
#include <cstddef>
#include <iostream>
#include <vector>

#include <mpi.h>

std::vector<double> binned_gaussian(100, 0.0);

double gaussian(hep::mc_point<double> const& sample)
{
	double result = 1.0;

	// acos(-1.0) = pi - note that there is no pi constant in C++
	double factor = 2.0 / std::sqrt(std::acos(-1.0));

	// multiply a (half-)gaussian for every dimension
	for (std::size_t i = 0; i != sample.point.size(); ++i)
	{
		double argument = sample.point[i];
		result *= std::exp(-argument * argument);
		result *= factor;

		// bin the zeroth gauss
		if (i == 0)
		{
			std::size_t index = sample.point[0] * binned_gaussian.size();
			binned_gaussian[index] += result * sample.weight;
		}
	}

	return result;
}

/*
 * Start this program on your local box with:
 *
 *     mpirun -np 8 ./mpi_vegas_example
 *
 * Change '8' to the number of processes this program should use.
 */
int main(int argc, char* argv[])
{
	// initialize MPI
	MPI_Init(&argc, &argv);

	// which index has this process?
	int rank = -1;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// only print in process 0 to prevent cluttering
	if (rank == 0)
	{
		// how many processes are started?
		int size = -1;
		MPI_Comm_size(MPI_COMM_WORLD, &size);

		std::cout << "running hep::mpi_vegas with " << size << " processes\n\n";
	}

	// set the verbose vegas callback function
	hep::mpi_vegas_callback<double>(hep::mpi_vegas_verbose_callback<double>);

	// perform 5 iteration with 10^7 calls each; the integrand is a ten
	// dimensional gaussian
	std::size_t iterations = 5;
	std::size_t calls = 10000000;
	std::size_t dimensions = 10;

	// perform the integration; this function will also call the callback
	// function (see above) after each iteration which in turn prints the
	// individual iteration results
	auto results = hep::mpi_vegas<double>(
		MPI_COMM_WORLD,
		dimensions,
		std::vector<std::size_t>(iterations, calls),
		gaussian
	);

	// take binned_gaussian of every process, add entries separately and write
	// back into it
	MPI_Allreduce(
		MPI_IN_PLACE,
		&(binned_gaussian[0]),
		binned_gaussian.size(),
		MPI_DOUBLE,
		MPI_SUM,
		MPI_COMM_WORLD
	);

	auto result = hep::cumulative_result0(results.begin(), results.end());
	double chi_square_dof = hep::chi_square_dof0(results.begin(),
		results.end());

	if (rank == 0)
	{
		std::cout << "cumulative result : " << result.value << " +- ";
		std::cout << result.error << "\n";
		std::cout << "chi^2/dof : " << chi_square_dof << "\n\n";

		std::cout << "# binned gaussian:\n";
		for (auto i = binned_gaussian.begin(); i != binned_gaussian.end(); ++i)
		{
			double average_function_value = *i;

			// average over the iterations
			average_function_value /= iterations;

			// the binned gaussian is the integral over the bin - the average is
			// obtained by dividing through the bin length or muplying by the
			// number of bins
			average_function_value *= binned_gaussian.size();

			// print result
			std::cout << average_function_value << "\n";
		}

		/* - run the program
		 * - save the program's output (after `# binned gaussian`) into a file
		 *   called 'data'
		 * - start gnuplot and type
		 *
		 *     plot 'data' using ($0/100.0):1, 2.0*exp(-x*x)/sqrt(pi)
		 *
		 *   to verify the result.
		 */
	}

	// clean up
	MPI_Finalize();

	return 0;
}
