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

#include "hep/mc/buffer_helper.hpp"
#include "hep/mc/generator_helper.hpp"
#include "hep/mc/mpi_helper.hpp"
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
template <class T, class F, class D, class P, class R = std::mt19937>
inline std::vector<multi_channel_result<T>> mpi_multi_channel(
	MPI_Comm communicator,
	std::size_t dimensions,
	std::size_t map_dimensions,
	std::vector<std::size_t> const& iteration_calls,
	F&& function,
	std::vector<T> const& channel_weights,
	D&& densities,
	P&& projector,
	R&& generator = std::mt19937()
) {
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
	std::size_t const usage = (1 + dimensions) * random_number_usage<T, R>();

	for (auto i = iteration_calls.begin(); i != iteration_calls.end(); ++i)
	{
		generator.discard(usage * discard_before(*i, rank, world));

		std::size_t const calls = (*i / world) +
			(static_cast <std::size_t> (rank) < (*i % world) ? 1 : 0);

		auto const result = multi_channel_iteration(
			dimensions,
			map_dimensions,
			calls,
			function,
			weights,
			densities,
			projector,
			generator
		);

		generator.discard(usage * discard_after(*i, calls, rank, world));

		auto const new_result = allreduce_result(
			communicator,
			result,
			buffer,
			result.adjustment_data()
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

		weights = multi_channel_refine_weights(weights, buffer);
	}

	return results;
}

/// Implements the MPI-parallelized adaptive multi channel algorithm. See
/// \ref multi_channel_group for a detailed description of the parameters.
template <typename T, typename F, typename D, typename R = std::mt19937>
inline std::vector<multi_channel_result<T>> mpi_multi_channel(
	MPI_Comm communicator,
	std::size_t dimensions,
	std::size_t map_dimensions,
	std::vector<std::size_t> const& iteration_calls,
	F&& function,
	std::size_t channels,
	D&& densities,
	R&& generator = std::mt19937()
) {
	return mpi_multi_channel(
		communicator,
		dimensions,
		map_dimensions,
		iteration_calls,
		std::forward<F>(function),
		std::vector<T>(channels, T(1.0) / T(channels)),
		std::forward<D>(densities),
		default_projector<T>(),
		std::forward<R>(generator)
	);
}

///
template <class T, class F, class D, class P, class R = std::mt19937>
inline std::vector<multi_channel_result<T>> mpi_multi_channel_distributions(
	MPI_Comm communicator,
	std::size_t dimensions,
	std::size_t map_dimensions,
	std::vector<std::size_t> const& iteration_calls,
	F&& function,
	std::size_t channels,
	D&& densities,
	P&& projector,
	R&& generator = std::mt19937()
) {
	return mpi_multi_channel(
		communicator,
		dimensions,
		map_dimensions,
		iteration_calls,
		std::forward<F>(function),
		std::vector<T>(channels, T(1.0) / T(channels)),
		std::forward<D>(densities),
		std::forward<P>(projector),
		std::forward<R>(generator)
	);
}

/// @}

}

#endif
