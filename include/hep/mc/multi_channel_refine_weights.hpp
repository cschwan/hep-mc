#ifndef HEP_MC_MULTI_CHANNEL_REFINE_WEIGHTS_HPP
#define HEP_MC_MULTI_CHANNEL_REFINE_WEIGHTS_HPP

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

#include <cmath>
#include <cstddef>
#include <vector>

namespace hep
{

/// Uses `adjustment_data` from a previous call of \ref multi_channel_iteration to refine `weights`.
/// The procedure is the one suggested in Ref. \cite WeightOptimization with the following
/// modifications:
/// - the weights are check if they are smaller then the given `minimum_weight`. If this is the case
///   they are set to the value of `minimum_weight` which, after the normalization of all weights,
///   make them a little smaller then the given minimum weight,
/// - `adjustment_data` is raised to the power given by `beta`. The reference given above suggest
///   `beta = 0.5`, but the default value is smaller which sometimes gives a more stable
///   convergence.
template <typename T>
inline std::vector<T> multi_channel_refine_weights(
    std::vector<T> const& weights,
    std::vector<T> const& adjustment_data,
    T minimum_weight,
    T beta
) {
    using std::fmax;
    using std::pow;

    std::vector<T> new_weights(weights.size());

    T sum_of_new_weights = T();

    for (std::size_t i = 0; i != new_weights.size(); ++i)
    {
        new_weights[i] = weights[i] * pow(adjustment_data[i], beta);
        sum_of_new_weights += new_weights[i];
    }

    T new_sum = T();

    for (T& weight : new_weights)
    {
        if (weight == T())
        {
            // do not enable disabled channels (with weight zero) by setting them to the minimum
            // weight
            continue;
        }

        weight /= sum_of_new_weights;
        weight = fmax(weight, minimum_weight);
        new_sum += weight;
    }

    for (T& weight : new_weights)
    {
        weight /= new_sum;
    }

    return new_weights;
}

}

#endif
