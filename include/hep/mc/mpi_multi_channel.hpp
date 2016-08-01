#ifndef HEP_MC_MPI_MULTI_CHANNEL_HPP
#define HEP_MC_MPI_MULTI_CHANNEL_HPP

/*
 * hep-mc - A Template Library for Monte Carlo Integration
 * Copyright (C) 2015-2016  Christopher Schwan
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
#include "hep/mc/integrand.hpp"
#include "hep/mc/mpi_multi_channel_callback.hpp"
#include "hep/mc/multi_channel.hpp"
#include "hep/mc/multi_channel_result.hpp"

#include <cstddef>
#include <random>
#include <utility>
#include <vector>

#include <mpi.h>

namespace hep
{

/// \addtogroup multi_channel_group
/// @{

/// Implements the MPI-parallelized adaptive multi channel algorithm. See
/// \ref multi_channel_group for a detailed description of the parameters.
template <typename I, typename R = std::mt19937>
inline std::vector<multi_channel_result<numeric_type_of<I>>> mpi_multi_channel(
	MPI_Comm communicator,
	I&& integrand,
	std::vector<std::size_t> const& iteration_calls,
	std::vector<numeric_type_of<I>> const& channel_weights,
	std::size_t min_calls_per_channel = 0,
	R&& generator = std::mt19937()
) {
	using T = numeric_type_of<I>;

	int rank = 0;
	MPI_Comm_rank(communicator, &rank);
	int world = 0;
	MPI_Comm_size(communicator, &world);

	auto weights = channel_weights;

	std::vector<multi_channel_result<T>> results;
	results.reserve(iteration_calls.size());

	std::vector<T> buffer;

	// hep::discrete_distribution consumes as many random numbers as an
	// additional dimension
	std::size_t const usage = (1 + integrand.dimensions()) *
		random_number_usage<T, R>();

	for (auto i = iteration_calls.begin(); i != iteration_calls.end(); ++i)
	{
		generator.discard(usage * discard_before(*i, rank, world));

		std::size_t const calls = (*i / world) +
			(static_cast <std::size_t> (rank) < (*i % world) ? 1 : 0);

		auto const result = multi_channel_iteration(
			integrand,
			calls,
			weights,
			generator
		);

		generator.discard(usage * discard_after(*i, calls, rank, world));

		auto const new_result = allreduce_result(
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
			buffer,
			weights
		);

		if (!mpi_multi_channel_callback<T>()(communicator, results))
		{
			break;
		}

		T const minimum_weight = T(min_calls_per_channel) / *i;

		weights = multi_channel_refine_weights(weights, buffer, minimum_weight);
	}

	return results;
}

/// Implements the MPI-parallelized adaptive multi channel algorithm. See
/// \ref multi_channel_group for a detailed description of the parameters.
template <typename I, typename R = std::mt19937>
inline std::vector<multi_channel_result<numeric_type_of<I>>> mpi_multi_channel(
	MPI_Comm communicator,
	I&& integrand,
	std::vector<std::size_t> const& iteration_calls,
	std::size_t min_calls_per_channel = 0,
	R&& generator = std::mt19937()
) {
	using T = numeric_type_of<I>;

	return mpi_multi_channel(
		communicator,
		std::forward<I>(integrand),
		iteration_calls,
		std::vector<T>(integrand.channels(), T(1.0) / T(integrand.channels())),
		min_calls_per_channel,
		std::forward<R>(generator)
	);
}

/// @}

}

#endif
