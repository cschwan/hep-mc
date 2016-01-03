#ifndef HEP_MC_BUFFER_HELPER_HPP
#define HEP_MC_BUFFER_HELPER_HPP

/*
 * hep-mc - A Template Library for Monte Carlo Integration
 * Copyright (C) 2016  Christopher Schwan
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

#include "hep/mc/distribution_result.hpp"
#include "hep/mc/mc_result.hpp"
#include "hep/mc/plain_result.hpp"

#include <cstddef>
#include <vector>

#include <mpi.h>

namespace
{

template <typename T>
hep::plain_result<T> allreduce_result(
	MPI_Comm communicator,
	hep::plain_result<T> const& result,
	std::vector<T>& buffer,
	std::vector<T> const& in_buffer
) {
	// pack everything into `buffer` ...
	buffer = in_buffer;
	buffer.push_back(result.sum());
	buffer.push_back(result.sum_of_squares());
	
	for (auto const distribution : result.distributions())
	{
		for (auto const bin : distribution.results())
		{
			buffer.push_back(bin.sum());
			buffer.push_back(bin.sum_of_squares());
		}
	}

	// merge `buffer` of all processes into `buffer` again ...
	MPI_Allreduce(
		MPI_IN_PLACE,
		&buffer[0],
		buffer.size(),
		mpi_datatype<T>(),
		MPI_SUM,
		communicator
	);

	// ... and extract the results
	std::size_t index = in_buffer.size();
	T sum = buffer[index++];
	T sum_of_squares = buffer[index++];

	std::vector<hep::distribution_result<T>> distributions;
	std::vector<hep::mc_result<T>> bins;
	for (std::size_t j = 0; j != result.distributions().size(); ++j)
	{
		bins.clear();
		auto const& distribution = result.distributions()[j];
		std::size_t const size = result.distributions().size();

		for (std::size_t k = 0; k != size; ++k)
		{
			bins.emplace_back(result.calls(), buffer[index], buffer[index + 1]);
			index += 2;
		}

		distributions.emplace_back(distribution.mid_points(), bins);
	}
	
	// resize `buffer` - contains the merged `additional_data`
	buffer.resize(in_buffer.size());

	return hep::plain_result<T>(
		distributions,
		result.calls(),
		sum,
		sum_of_squares
	);
}

}

#endif
