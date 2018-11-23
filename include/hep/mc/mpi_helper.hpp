#ifndef HEP_MC_MPI_HELPER_HPP
#define HEP_MC_MPI_HELPER_HPP

/*
 * hep-mc - A Template Library for Monte Carlo Integration
 * Copyright (C) 2016-2018  Christopher Schwan
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

namespace hep
{

/// \cond INTERNAL

template <typename T>
MPI_Datatype mpi_datatype();

template <>
inline MPI_Datatype mpi_datatype<unsigned int>()
{
    return MPI_UNSIGNED;
}

template <>
inline MPI_Datatype mpi_datatype<unsigned long int>()
{
    return MPI_UNSIGNED_LONG;
}

template <>
inline MPI_Datatype mpi_datatype<unsigned long long int>()
{
    return MPI_UNSIGNED_LONG_LONG;
}

template <>
inline MPI_Datatype mpi_datatype<float>()
{
    return MPI_FLOAT;
}

template <>
inline MPI_Datatype mpi_datatype<double>()
{
    return MPI_DOUBLE;
}

template <>
inline MPI_Datatype mpi_datatype<long double>()
{
    return MPI_LONG_DOUBLE;
}

template <typename T>
hep::plain_result<T> allreduce_result(
    MPI_Comm communicator,
    hep::plain_result<T> const& result,
    std::vector<T>& buffer,
    std::vector<T> const& in_buffer,
    std::size_t total_calls
) {
    // pack everything into `buffer` ...
    buffer = in_buffer;
    buffer.push_back(result.sum());
    buffer.push_back(result.sum_of_squares());

    std::vector<std::size_t> size_t_buffer;
    size_t_buffer.push_back(result.non_zero_calls());
    size_t_buffer.push_back(result.finite_calls());

    for (auto const& distribution : result.distributions())
    {
        for (auto const& bin : distribution.results())
        {
            buffer.push_back(bin.sum());
            buffer.push_back(bin.sum_of_squares());
            size_t_buffer.push_back(bin.non_zero_calls());
            size_t_buffer.push_back(bin.finite_calls());
        }
    }

    // ... merge `buffer` of all processes into `buffer` again ...
    MPI_Allreduce(
        MPI_IN_PLACE,
        &buffer[0],
        buffer.size(),
        mpi_datatype<T>(),
        MPI_SUM,
        communicator
    );

    MPI_Allreduce(
        MPI_IN_PLACE,
        &size_t_buffer[0],
        size_t_buffer.size(),
        mpi_datatype<std::size_t>(),
        MPI_SUM,
        communicator
    );

    // ... and extract the results
    std::size_t index = in_buffer.size();
    T sum = buffer[index++];
    T sum_of_squares = buffer[index++];

    std::vector<hep::distribution_result<T>> distributions;
    for (auto const& distribution : result.distributions())
    {
        std::vector<hep::mc_result<T>> bins;
        bins.reserve(distribution.results().size());

        for (std::size_t i = 0; i != distribution.results().size(); ++i)
        {
            bins.emplace_back(
                total_calls,
                size_t_buffer[index     - in_buffer.size()],
                size_t_buffer[index + 1 - in_buffer.size()],
                buffer[index],
                buffer[index + 1]
            );
            index += 2;
        }

        distributions.emplace_back(distribution.parameters(), bins);
    }

    // resize `buffer` - contains the merged `additional_data`
    buffer.resize(in_buffer.size());

    return hep::plain_result<T>(
        distributions,
        total_calls,
        size_t_buffer[0],
        size_t_buffer[1],
        sum,
        sum_of_squares
    );
}

/// \endcond

}

#endif
