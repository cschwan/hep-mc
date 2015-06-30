#ifndef HEP_MC_MPI_MULTI_CHANNEL_HPP
#define HEP_MC_MPI_MULTI_CHANNEL_HPP

/*
 * hep-mc - A Template Library for Monte Carlo Integration
 * Copyright (C) 2015  Christopher Schwan
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
#include "hep/mc/mpi_helper.hpp"
#include "hep/mc/mpi_multi_channel_callback.hpp"
#include "hep/mc/multi_channel.hpp"
#include "hep/mc/multi_channel_result.hpp"

#include <cstddef>
#include <random>
#include <vector>

#include <mpi.h>

namespace hep
{

template <typename T, typename F, typename D, typename R = std::mt19937>
inline std::vector<multi_channel_result<T>> mpi_multi_channel(
	MPI_Comm communicator,
	std::size_t dimensions,
	std::size_t map_dimensions,
	std::vector<std::size_t> iteration_calls,
	F&& function,
	std::size_t channels,
	D&& densities,
	R&& generator = std::mt19937()
) {
	int rank = 0;
	MPI_Comm_rank(communicator, &rank);
	int world = 0;
	MPI_Comm_size(communicator, &world);

	std::vector<T> weights(channels, T(1.0) / T(channels));

	std::vector<multi_channel_result<T>> results;
	results.reserve(iteration_calls.size());

	// FIXME: this assumes std::discrete_distribution takes two random numbers
	// for the selection of a random channel; this is probably not true for all
	// implementations
	std::size_t const usage = 2 + dimensions * random_number_usage<T, R>();

	for (auto i = iteration_calls.begin(); i != iteration_calls.end(); ++i)
	{
		generator.discard(usage * discard_before(*i, rank, world));

		std::size_t const calls = (*i / world) +
			(static_cast <std::size_t> (rank) < (*i % world) ? 1 : 0);

		auto result = multi_channel_iteration(
			dimensions,
			map_dimensions,
			calls,
			*i,
			function,
			weights,
			densities,
			generator
		);

		generator.discard(usage * discard_after(*i, calls, rank, world));

		MPI_Allreduce(
			MPI_IN_PLACE,
			&(result.adjustment_data[0]),
			result.adjustment_data.size(),
			mpi_datatype<T>(),
			MPI_SUM,
			communicator
		);

		results.emplace_back(*i, result.adjustment_data, weights);

		if (!mpi_multi_channel_callback<T>()(communicator, results))
		{
			break;
		}

		weights = multi_channel_refine_weights(weights, result.adjustment_data);
	}

	return results;
}

}

#endif
