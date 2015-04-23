#ifndef HEP_MC_MPI_VEGAS_HPP
#define HEP_MC_MPI_VEGAS_HPP

/*
 * hep-mc - A Template Library for Monte Carlo Integration
 * Copyright (C) 2013-2015  Christopher Schwan
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "hep/mc/generator_helper.hpp"
#include "hep/mc/global_configuration.hpp"
#include "hep/mc/mpi_helper.hpp"
#include "hep/mc/mpi_vegas_callback.hpp"
#include "hep/mc/vegas.hpp"
#include "hep/mc/vegas_iteration_result.hpp"
#include "hep/mc/vegas_pdf.hpp"

#include <cstddef>
#include <random>
#include <utility>
#include <vector>

#include <mpi.h>

namespace hep
{

/// \addtogroup vegas_group
/// @{

/**
 * Implements the MPI-parallelized VEGAS algorithm. This function can be used to start from an
 * already adapted grid, e.g. one by \ref vegas_iteration_result.pdf obtained by a previous \ref
 * vegas call.
 */
template <typename T, typename F, typename R = std::mt19937>
inline std::vector<vegas_iteration_result<T>> mpi_vegas(
	MPI_Comm communicator,
	std::vector<std::size_t> const& iteration_calls,
	F&& function,
	vegas_pdf<T> const& start_pdf,
	T alpha = T(1.5),
	R&& generator = std::mt19937()
) {
	int rank = 0;
	MPI_Comm_rank(communicator, &rank);
	int world = 0;
	MPI_Comm_size(communicator, &world);

	// create a fresh grid
	auto pdf = start_pdf;

	// vector holding all iteration results
	std::vector<vegas_iteration_result<T>> results;
	results.reserve(iteration_calls.size());

	std::size_t const usage = pdf.dimensions() * random_number_usage<T, R>();

	// perform iterations
	for (auto i = iteration_calls.begin(); i != iteration_calls.end(); ++i)
	{
		generator.discard(usage * discard_before(*i, rank, world));

		std::size_t const calls = (*i / world) +
			(static_cast <std::size_t> (rank) < (*i % world) ? 1 : 0);
		auto result = vegas_iteration(calls, *i, pdf, function, generator);

		generator.discard(usage * discard_after(*i, calls, rank, world));

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
		results.push_back(vegas_iteration_result<T>(*i, pdf, result.adjustment_data));

		if (!mpi_vegas_callback<T>()(communicator, results))
		{
			break;
		}

		pdf = vegas_refine_pdf(pdf, alpha, result.adjustment_data);
	}

	return results;
}

/**
 * Implements the MPI-parallelized VEGAS algorithm. See \ref vegas for a more detailed description
 * on the VEGAS algorithm. In contrast to the single-threaded versions this function makes sure that
 * every random number generator is seeded differently so every MPI process yields an independent
 * result. After each iteration the intermediate results are passed to the function set by \ref
 * mpi_vegas_callback which can e.g. be used to print them out. The callback function is able to
 * stop the integration if it returns `false`. In this case less iterations are performed than
 * requested.
 *
 * \param communicator The MPI communicator that is used to communicate between the different MPI
 *        processes.
 * \param dimensions The number of parameters `function` accepts.
 * \param iteration_calls The number of function calls that are used to obtain a result for each
 *        iteration. `iteration_calls.size()` determines the number of iterations.
 * \param function The function that will be integrated over the hypercube. See \ref integrands for
 *        further explanation.
 * \param bins The number of bins that the grid will contain for each dimension.
 * \param alpha The \f$ \alpha \f$ parameter of VEGAS. This parameter is usually set between `1` and
 *        `2`.
 * \param generator The random number generator that will be used to generate random points from the
 * hypercube. This generator is properly seeded.
 */
template <typename T, typename F, typename R = std::mt19937>
inline std::vector<vegas_iteration_result<T>> mpi_vegas(
	MPI_Comm communicator,
	std::size_t dimensions,
	std::vector<std::size_t> const& iteration_calls,
	F&& function,
	std::size_t bins = 128,
	T alpha = T(1.5),
	R&& generator = std::mt19937()
) {
	return mpi_vegas(
		communicator,
		iteration_calls,
		function,
		vegas_pdf<T>(dimensions, bins),
		alpha,
		std::forward<R>(generator)
	);
}

/// @}

}

#endif
