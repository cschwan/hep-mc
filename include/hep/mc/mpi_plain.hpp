#ifndef HEP_MC_MPI_PLAIN_HPP
#define HEP_MC_MPI_PLAIN_HPP

/*
 * hep-mc - A Template Library for Monte Carlo Integration
 * Copyright (C) 2013-2019  Christopher Schwan
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
#include "hep/mc/integrand.hpp"
#include "hep/mc/mpi_callback.hpp"
#include "hep/mc/mpi_helper.hpp"
#include "hep/mc/plain.hpp"
#include "hep/mc/plain_chkpt.hpp"
#include "hep/mc/plain_result.hpp"

#include <mpi.h>

#include <cstddef>
#include <random>
#include <vector>

namespace hep
{

/// \addtogroup plain_group
/// @{

/// MPI version of the plain Monte Carlo integrator.
template <typename I, typename Checkpoint = default_plain_chkpt<numeric_type_of<I>>,
    typename Callback = mpi_callback<Checkpoint>>
inline Checkpoint mpi_plain(
    MPI_Comm communicator,
    I&& integrand,
    std::vector<std::size_t> const& iteration_calls,
    Checkpoint chkpt = make_plain_chkpt<numeric_type_of<I>>(),
    Callback callback = mpi_callback<Checkpoint>()
) {
    using T = numeric_type_of<I>;

    int rank = 0;
    MPI_Comm_rank(communicator, &rank);
    int world = 0;
    MPI_Comm_size(communicator, &world);

    auto generator = chkpt.generator();

    std::vector<T> buffer;

    std::size_t const usage = integrand.dimensions() *
        random_number_usage<T, decltype (generator)>();

    // perform iterations
    for (auto calls : iteration_calls)
    {
        generator.discard(usage * discard_before(calls, rank, world));

        // the number of function calls for each MPI process
        std::size_t const sub_calls = (calls / world) +
            (static_cast <std::size_t> (rank) < (calls % world) ? 1 : 0);
        auto const result = plain_iteration(integrand, sub_calls, generator);

        generator.discard(usage * discard_after(calls, sub_calls, rank, world));

        auto const& new_result = allreduce_result(
            communicator,
            result,
            buffer,
            std::vector<T>(),
            calls
        );

        chkpt.add(new_result, generator);

        if (!callback(communicator, chkpt))
        {
            break;
        }
    }

    return chkpt;
}

/// @}

}

#endif
