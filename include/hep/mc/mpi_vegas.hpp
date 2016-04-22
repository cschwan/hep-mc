#ifndef HEP_MC_MPI_VEGAS_HPP
#define HEP_MC_MPI_VEGAS_HPP

/*
 * hep-mc - A Template Library for Monte Carlo Integration
 * Copyright (C) 2013-2016  Christopher Schwan
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

#include "hep/mc/internal/generator_helper.hpp"
#include "hep/mc/internal/mpi_helper.hpp"
#include "hep/mc/global_configuration.hpp"
#include "hep/mc/integrand.hpp"
#include "hep/mc/mpi_vegas_callback.hpp"
#include "hep/mc/vegas.hpp"
#include "hep/mc/vegas_pdf.hpp"
#include "hep/mc/vegas_result.hpp"

#include <cstddef>
#include <random>
#include <utility>
#include <vector>

#include <mpi.h>

namespace hep
{

/// \addtogroup vegas_group
/// @{

/// Implements the MPI-parallelized VEGAS algorithm. This function can be used
/// to start from an already adapted grid, e.g. one by \ref vegas_result.pdf
/// obtained by a previous \ref vegas call.
template <typename I, typename R = std::mt19937>
inline std::vector<vegas_result<numeric_type_of<I>>> mpi_vegas(
	MPI_Comm communicator,
	I&& integrand,
	std::vector<std::size_t> const& iteration_calls,
	vegas_pdf<numeric_type_of<I>> const& start_pdf,
	numeric_type_of<I> alpha = numeric_type_of<I>(1.5),
	R&& generator = std::mt19937()
) {
	using T = numeric_type_of<I>;

	int rank = 0;
	MPI_Comm_rank(communicator, &rank);
	int world = 0;
	MPI_Comm_size(communicator, &world);

	// create a fresh grid
	auto pdf = start_pdf;

	// vector holding all iteration results
	std::vector<vegas_result<T>> results;
	results.reserve(iteration_calls.size());

	// reserve a buffer for the MPI call to sum `adjustment_data`, `sum`, and
	// `sum_of_squares`
	std::vector<T> buffer(pdf.dimensions() * pdf.bins() + 2);

	std::size_t const usage = pdf.dimensions() * random_number_usage<T, R>();

	// perform iterations
	for (auto i = iteration_calls.begin(); i != iteration_calls.end(); ++i)
	{
		generator.discard(usage * discard_before(*i, rank, world));

		std::size_t const calls = (*i / world) +
			(static_cast <std::size_t> (rank) < (*i % world) ? 1 : 0);
		auto const result = vegas_iteration(integrand, calls, pdf, generator);

		generator.discard(usage * discard_after(*i, calls, rank, world));

		auto const& new_result = allreduce_result(
			communicator,
			result,
			buffer,
			result.adjustment_data(),
			*i
		);

		results.emplace_back(
			new_result.distributions(),
			new_result.calls(),
			new_result.sum(),
			new_result.sum_of_squares(),
			pdf,
			buffer
		);

		if (!mpi_vegas_callback<T>()(communicator, results))
		{
			break;
		}

		pdf = vegas_refine_pdf(pdf, alpha, buffer);
	}

	return results;
}

/// Implements the MPI-parallelized VEGAS algorithm. See \ref vegas for a more
/// detailed description on the VEGAS algorithm. In contrast to the
/// single-threaded versions this function makes sure that every random number
/// generator is seeded differently so every MPI process yields an independent
/// result. After each iteration the intermediate results are passed to the
/// function set by \ref mpi_vegas_callback which can e.g. be used to print them
/// out. The callback function is able to stop the integration if it returns
/// `false`. In this case less iterations are performed than requested.
template <typename I, typename R = std::mt19937>
inline std::vector<vegas_result<numeric_type_of<I>>> mpi_vegas(
	MPI_Comm communicator,
	I&& integrand,
	std::vector<std::size_t> const& iteration_calls,
	std::size_t bins = 128,
	numeric_type_of<I> alpha = numeric_type_of<I>(1.5),
	R&& generator = std::mt19937()
) {
	using T = numeric_type_of<I>;

	return mpi_vegas(
		communicator,
		std::forward<I>(integrand),
		iteration_calls,
		vegas_pdf<T>(integrand.dimensions(), bins),
		alpha,
		std::forward<R>(generator)
	);
}

/// @}

}

#endif
