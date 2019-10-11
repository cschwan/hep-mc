#ifndef HEP_PS_PLAIN_CHKPT_HPP
#define HEP_PS_PLAIN_CHKPT_HPP

/*
 * hep-mc - A Template Library for Monte Carlo Integration
 * Copyright (C) 2019  Christopher Schwan
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

#include "hep/mc/chkpt.hpp"
#include "hep/mc/plain_result.hpp"

#include <iosfwd>
#include <random>

namespace hep
{

/// \addtogroup checkpoints
/// @{

/// Checkpoint created and used by the \ref plain_group.
template <typename T>
using plain_chkpt = chkpt<plain_result<T>>;

/// Checkpoint with random number generators created by using the \ref plain_group.
template <typename RandomNumberEngine, typename T>
using plain_chkpt_with_rng = chkpt_with_rng<RandomNumberEngine, plain_chkpt<T>>;

/// Helper function to create an initial checkpoint to start the \ref plain_group.
template <typename T, typename RandomNumberEngine = std::mt19937>
plain_chkpt_with_rng<RandomNumberEngine, T> make_plain_chkpt(
    RandomNumberEngine const& generator = RandomNumberEngine()
) {
    return plain_chkpt_with_rng<RandomNumberEngine, T>{generator};
}

/// Helper function to create a checkpoint by reading from the stream `in`.
template <typename T, typename RandomNumberEngine>
plain_chkpt_with_rng<RandomNumberEngine, T> make_plain_chkpt(std::istream& in)
{
    if (in.peek() == std::istream::traits_type::eof())
    {
        return make_plain_chkpt<T, RandomNumberEngine>();
    }

    return plain_chkpt_with_rng<RandomNumberEngine, T>{in};
}

/// Return type of \ref make_plain_chkpt with default arguments.
template <typename T>
using default_plain_chkpt = decltype (make_plain_chkpt<T>());

/// @}

}

#endif
