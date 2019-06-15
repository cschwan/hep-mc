#ifndef HEP_MC_MPI_MULTI_CHANNEL_HPP
#define HEP_MC_MPI_MULTI_CHANNEL_HPP

/*
 * hep-mc - A Template Library for Monte Carlo Integration
 * Copyright (C) 2015-2019  Christopher Schwan
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

/// MPI version of \ref multi_channel.
template <typename I, typename Checkpoint = default_multi_channel_chkpt<numeric_type_of<I>>,
    typename Callback = mpi_callback<Checkpoint>>
inline Checkpoint mpi_multi_channel(
    MPI_Comm communicator,
    I&& integrand,
    std::vector<std::size_t> const& iteration_calls,
    Checkpoint chkpt = make_multi_channel_chkpt<numeric_type_of<I>>(),
    Callback callback = mpi_callback<Checkpoint>()
) {
    using T = numeric_type_of<I>;

    chkpt.channels(integrand.channels());

    int rank = 0;
    MPI_Comm_rank(communicator, &rank);
    int world = 0;
    MPI_Comm_size(communicator, &world);

    auto generator = chkpt.generator();
    auto weights = chkpt.channel_weights();

    std::vector<T> buffer;

    // hep::discrete_distribution consumes as many random numbers as an
    // additional dimension
    std::size_t const usage = (1 + integrand.dimensions()) *
        random_number_usage<T, decltype (generator)>();

    for (auto const calls : iteration_calls)
    {
        generator.discard(usage * discard_before(calls, rank, world));

        std::size_t const sub_calls = (calls / world) +
            (static_cast <std::size_t> (rank) < (calls % world) ? 1 : 0);

        auto const sub_result = multi_channel_iteration(integrand, sub_calls, weights, generator);

        generator.discard(usage * discard_after(calls, sub_calls, rank, world));

        auto const result = multi_channel_result<T>{allreduce_result(communicator, sub_result,
            buffer, sub_result.adjustment_data(), calls), buffer, weights};

        chkpt.add(result, generator);

        if (!callback(communicator, chkpt))
        {
            break;
        }

        weights = multi_channel_refine_weights(weights, result.adjustment_data(),
            chkpt.min_weight(), chkpt.beta());
    }

    return chkpt;
}

/// @}

}

#endif
