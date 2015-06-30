#ifndef HEP_MC_MPI_MULTI_CHANNEL_CALLBACK_HPP
#define HEP_MC_MPI_MULTI_CHANNEL_CALLBACK_HPP

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

#include "hep/mc/multi_channel_callback.hpp"
#include "hep/mc/multi_channel_result.hpp"

#include <functional>
#include <vector>

#include <mpi.h>

namespace hep
{

template <typename T>
inline bool mpi_multi_channel_default_callback(
	MPI_Comm,
	std::vector<multi_channel_result<T>> const&
) {
	return true;
}

template <typename T>
inline bool mpi_multi_channel_verbose_callback(
	MPI_Comm communicator,
	std::vector<multi_channel_result<T>> const& results
) {
	int rank = -1;
	MPI_Comm_rank(communicator, &rank);

	if (rank == 0)
	{
		multi_channel_verbose_callback<T>(results);
	}

	return true;
}

template <typename T>
using mpi_multi_channel_callback_type =
	std::function<bool(MPI_Comm, std::vector<multi_channel_result<T>>)>;

template <typename T>
inline mpi_multi_channel_callback_type<T> mpi_multi_channel_callback(
	mpi_multi_channel_callback_type<T> callback = nullptr
) {
	static mpi_multi_channel_callback_type<T> object =
		mpi_multi_channel_default_callback<T>;

	if (callback != nullptr)
	{
		object = callback;
	}

	return object;
}

}

#endif

