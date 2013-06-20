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

	// set the verbose vegas callback function
	hep::mpi_vegas_callback<double>(hep::mpi_vegas_verbose_callback<double>);

	// perform 5 iteration with 10^7 calls each; this function will also call
	// mpi_vegas_verbose_callback after each iteration which in turn prints the
	// individual iterations
	auto results = hep::mpi_vegas<double>(
		MPI_COMM_WORLD,
		10,                                    // 10 dimensions
		std::vector<std::size_t>(5, 10000000),
		gauss
	);

	auto result = hep::cumulative_result<double>(results.begin(),
		results.end());
	double chi_square_dof = hep::chi_square_dof<double>(results.begin(),
		results.end());

	std::cout << "cumulative result : " << result.value << " +- " <<
		result.error << "\n";
	std::cout << "chi^2/dof : " << chi_square_dof << "\n";

	// clean up
	MPI_Finalize();

	return 0;
}
