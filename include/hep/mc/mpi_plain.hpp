#ifndef HEP_MC_MPI_PLAIN_HPP
#define HEP_MC_MPI_PLAIN_HPP

/*
 * hep-mc - A Template Library for Monte Carlo Integration
 * Copyright (C) 2013-2018  Christopher Schwan
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
#include "hep/mc/integrand.hpp"
#include "hep/mc/mpi_helper.hpp"
#include "hep/mc/mpi_plain_callback.hpp"
#include "hep/mc/plain.hpp"
#include "hep/mc/plain_result.hpp"

#include <cstddef>
#include <random>
#include <vector>

#include <mpi.h>

namespace hep
{

/// \addtogroup plain_group
/// @{

/// MPI-parallelized PLAIN Monte Carlo integrator. This function integrates `function` over the
/// unit-hypercube with the specified `dimensions` using `calls` function evaluations at randomly
/// chosen points determined by `generator`. The generator is not seeded.
///
/// \param communicator The MPI communicator that is used to communicate between the different MPI
/// processes.
/// \param integrand The function that will be integrated over the hypercube. See \ref integrands
/// for further explanation.
/// \param iteration_calls The number of function calls that are used to obtain the result.
/// \param generator The random number generator that will be used to generate random points from
/// the hypercube. This generator is properly seeded.
template <typename I, typename R = std::mt19937>
inline std::vector<plain_result<numeric_type_of<I>>> mpi_plain(
    MPI_Comm communicator,
    I&& integrand,
    std::vector<std::size_t> const& iteration_calls,
    R&& generator = std::mt19937()
) {
    using T = numeric_type_of<I>;

    int rank = 0;
    MPI_Comm_rank(communicator, &rank);
    int world = 0;
    MPI_Comm_size(communicator, &world);

    // vector holding all iteration results
    std::vector<plain_result<T>> results;
    results.reserve(iteration_calls.size());

    std::vector<T> buffer;

    std::size_t const usage = integrand.dimensions() *
        random_number_usage<T, decltype (generator)>();

    // perform iterations
    for (auto i = iteration_calls.begin(); i != iteration_calls.end(); ++i)
    {
        generator.discard(usage * discard_before(*i, rank, world));

        // the number of function calls for each MPI process
        std::size_t const calls = (*i / world) +
            (static_cast <std::size_t> (rank) < (*i % world) ? 1 : 0);
        auto const result = plain_iteration(integrand, calls, generator);

        generator.discard(usage * discard_after(*i, calls, rank, world));

        auto const& new_result = allreduce_result(
            communicator,
            result,
            buffer,
            std::vector<T>(),
            *i
        );

        results.push_back(new_result);

        if (!mpi_plain_callback<T>()(communicator, results))
        {
            break;
        }
    }

    return results;
}

/// @}

}

#endif
