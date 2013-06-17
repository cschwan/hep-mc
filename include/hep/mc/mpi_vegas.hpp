#ifndef HEP_MC_MPI_VEGAS_HPP
#define HEP_MC_MPI_VEGAS_HPP

#include <hep/mc/mpi_datatype.hpp>
#include <hep/mc/vegas.hpp>

#include <cstddef>
#include <random>
#include <vector>

#include <mpi.h>

namespace hep
{

/**
 * \ingroup algorithms
 *
 * Implements the MPI-parallelized VEGAS algorithm. See \ref vegas() for a more
 * detailed description on the VEGAS algorithm. In contrast to the
 * single-threaded versions this function makes sure that every random number
 * generator is seeded differently so every MPI process yields an independent
 * result.
 *
 * \see The file mpi_vegas.cpp provides an example
 *
 * \param communicator The MPI communicator that is used to communicate between
 *        the different MPI processes.
 * \param dimensions The number of parameters `function` accepts.
 * \param iteration_calls The number of function calls that are used to obtain
 *        a result for each iteration. `iteration_calls.size()` determines the
 *        number of iterations.
 * \param function The function that will be integrated over the hypercube. See
 *        \ref integrands for further explanation.
 * \param bins The number of bins that the grid will contain for each dimension.
 * \param alpha The \f$ \alpha \f$ parameter of VEGAS. This parameter should be
 *        between `1` and `2`.
 * \param generator The random number generator that will be used to generate
 *        random points from the hypercube. This generator is not seeded.
 */
template <typename T, typename F, typename R = std::mt19937>
std::vector<vegas_iteration_result<T>> mpi_vegas(
	MPI_Comm communicator,
	std::size_t dimensions,
	std::vector<std::size_t> const& iteration_calls,
	F&& function,
	std::size_t bins = 128,
	T alpha = T(1.5),
	R&& generator = std::mt19937()
) {
	int rank = 0;
	MPI_Comm_rank(communicator, &rank);
	int world = 0;
	MPI_Comm_size(communicator, &world);

	// seed every random number generator differently
	std::size_t r = rank * 10;
	std::seed_seq sequence{r+0, r+1, r+2, r+3, r+4, r+5, r+6, r+7, r+8, r+9};
	generator.seed(sequence);

	// create a fresh grid
	auto grid = vegas_grid<T>(dimensions, bins);

	// vector holding all iteration results
	std::vector<vegas_iteration_result<T>> results;
	results.reserve(iteration_calls.size());

	// perform iterations
	for (auto i = iteration_calls.begin(); i != iteration_calls.end(); ++i)
	{
		std::size_t const calls = (*i / world) +
			(static_cast <std::size_t> (rank) < (*i % world) ? 1 : 0);
		auto result = vegas_iteration(calls, *i, grid, function, generator);

		// add up results
		MPI_Allreduce(
			MPI_IN_PLACE,
			&(result.adjustment_data[0]),
			result.adjustment_data.size(),
			mpi_datatype<T>(),
			MPI_SUM,
			communicator
		);

		// calculate accumulated results
		results.push_back(vegas_iteration_result<T>(*i, grid,
			result.adjustment_data));

		grid = vegas_adjust_grid(alpha, grid, result.adjustment_data);
	}

	return results;
}

}

#endif
