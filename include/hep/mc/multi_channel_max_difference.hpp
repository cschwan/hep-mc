#ifndef HEP_MC_MULTI_CHANNEL_MAX_DIFFERENCE_HPP
#define HEP_MC_MULTI_CHANNEL_MAX_DIFFERENCE_HPP

/*
 * hep-mc - A Template Library for Monte Carlo Integration
 * Copyright (C) 2015-2018  Christopher Schwan
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

#include "hep/mc/multi_channel_result.hpp"

#include <cmath>
#include <cstddef>

namespace hep
{

/// \addtogroup results
/// @{

/// Returns the maximum difference \f$ D \f$ defined in Ref. \cite WeightOptimization as
/// \f[
///     D = \max_{i,j} | W_i ( \alpha ) - W_j ( \alpha ) |
/// \f]
/// with \f$ W_i ( \alpha ) \f$ being stored in `adjustment_data` of of the given `result`.
template <typename T>
inline T multi_channel_max_difference(multi_channel_result<T> const& result)
{
    using std::fabs;
    using std::fmax;

    T max = T();

    for (std::size_t i = 0; i != result.adjustment_data().size() - 1; ++i)
    {
        for (std::size_t j = i + 1; j != result.adjustment_data().size(); ++j)
        {
            T const wi = result.adjustment_data()[i];
            T const wj = result.adjustment_data()[j];

            max = fmax(max, fabs(wi - wj));
        }
    }

    return max;
}

/// @}

}

#endif
