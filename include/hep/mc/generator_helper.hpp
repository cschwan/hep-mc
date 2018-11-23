#ifndef HEP_MC_GENERATOR_HELPER_HPP
#define HEP_MC_GENERATOR_HELPER_HPP

/*
 * hep-mc - A Template Library for Monte Carlo Integration
 * Copyright (C) 2014-2018  Christopher Schwan
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

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <type_traits>

namespace hep
{

/// \cond INTERNAL

template <typename T, typename R>
inline std::size_t random_number_usage()
{
    using S = typename std::remove_reference<R>::type;

    // the number of random bits
    std::size_t const b = std::numeric_limits<T>::digits;

    // the number of different numbers the generator can generate
    long double const r = static_cast <long double> (S::max())
        - static_cast <long double> (S::min()) + 1.0L;

    // the number of bits needed to hold the value of 'r'
    std::size_t const log2r = std::log2(r);

    std::size_t const k = std::max<std::size_t>(1, (b + log2r - 1UL) / log2r);

    return k;
}

inline std::size_t discard_before(std::size_t total_calls, std::size_t rank, std::size_t world)
{
    std::size_t const before = (total_calls / world) * rank +
        (((total_calls % world) < rank) ? total_calls % world : rank);

    return before;
}

inline std::size_t discard_after(
    std::size_t total_calls,
    std::size_t calls,
    std::size_t rank,
    std::size_t world
) {
    std::size_t const before = discard_before(total_calls, rank, world);
    std::size_t const after = ((before + calls) < total_calls) ?
        (total_calls - before - calls) : 0;

    return after;
}

/// \endcond

}

#endif
