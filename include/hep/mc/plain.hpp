#ifndef HEP_MC_PLAIN_HPP
#define HEP_MC_PLAIN_HPP

/*
 * hep-mc - A Template Library for Monte Carlo Integration
 * Copyright (C) 2012-2019  Christopher Schwan
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

#include "hep/mc/accumulator.hpp"
#include "hep/mc/callback.hpp"
#include "hep/mc/integrand.hpp"
#include "hep/mc/mc_point.hpp"
#include "hep/mc/plain_chkpt.hpp"
#include "hep/mc/plain_result.hpp"

#include <cstddef>
#include <limits>
#include <random>
#include <vector>

namespace hep
{

/// \addtogroup plain_group
/// @{

/// Performs exactly one iteration using the PLAIN Monte Carlo integration algorithm.
template <typename I, typename R>
inline plain_result<numeric_type_of<I>> plain_iteration(
    I&& integrand,
    std::size_t calls,
    R& generator
) {
    using T = numeric_type_of<I>;

    // the accumulator takes care of the actual evaluation of the integrand and the generation of
    // possible distribution(s)
    auto accumulator = make_accumulator(integrand);

    // storage for random numbers
    std::vector<T> random_numbers(integrand.dimensions());

    // perform as many calls as requested
    for (std::size_t i = 0; i != calls; ++i)
    {
        // fill container with random numbers
        for (std::size_t j = 0; j != integrand.dimensions(); ++j)
        {
            random_numbers[j] = std::generate_canonical<T,
                std::numeric_limits<T>::digits>(generator);
        }

        mc_point<T> const point(random_numbers);

        // evaluate the integrand with the specified point and, if there are any distributions
        // requested, take care of them as well
        accumulator.invoke(integrand, point);
    }

    return accumulator.result(calls);
}

/// PLAIN Monte Carlo integrator. This function integrates `integrand` over the unit-hypercube using
/// `calls` function evaluations with randomly chosen points determined by `generator`.
template <typename I, typename Checkpoint = default_plain_chkpt<numeric_type_of<I>>,
    typename Callback = callback<Checkpoint>>
inline Checkpoint plain(
    I&& integrand,
    std::vector<std::size_t> const& iteration_calls,
    Checkpoint chkpt = make_plain_chkpt<numeric_type_of<I>>(),
    Callback callback = hep::callback<Checkpoint>()
) {
    auto generator = chkpt.generator();

    // perform iterations
    for (auto const calls : iteration_calls)
    {
        auto const result = plain_iteration(integrand, calls, generator);

        chkpt.add(result, generator);

        if (!callback(chkpt))
        {
            break;
        }
    }

    return chkpt;
}

/// @}

}

#endif
