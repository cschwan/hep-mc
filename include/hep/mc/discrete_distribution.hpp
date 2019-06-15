#ifndef HEP_MC_DISCRETE_DISTRIBUTION_HPP
#define HEP_MC_DISCRETE_DISTRIBUTION_HPP

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

#include <algorithm>
#include <iterator>
#include <limits>
#include <numeric>
#include <random>
#include <vector>

namespace hep
{

/// \cond INTERNAL

// Implements a subset of the functionality of `std::discrete_distribution`, but it uses
// `std::generate_canonical<T, ...>()` exactly once per random integer.
template <typename I = int, typename T = double>
class discrete_distribution
{
public:
    /// Constructor. Creates a new object using the weights pointed to by the range given with
    /// `begin` and `end`.
    template <typename Iterator>
    discrete_distribution(Iterator begin, Iterator end)
        : weight_sums(std::distance(begin, end))
    {
        std::partial_sum(begin, end, weight_sums.begin());

        // normalize the sums
        for (auto& k : weight_sums)
        {
            k /= weight_sums.back();
        }
    }

    /// Creates a new random integer using the specified random number generator. The integer is
    /// generated using exactly one call to `std::generate_canonical`.
    template <typename R>
    I operator()(R& generator) const
    {
        T const value = std::generate_canonical<T, std::numeric_limits<T>::digits>(generator);

        auto const iterator = std::lower_bound(weight_sums.begin(), weight_sums.end(), value);

        I const result = std::distance(weight_sums.begin(), iterator);

        return result;
    }

private:
    std::vector<T> weight_sums;
};

/// \endcond

}

#endif
