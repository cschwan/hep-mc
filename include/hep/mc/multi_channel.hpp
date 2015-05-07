#ifndef HEP_MC_MULTI_CHANNEL_HPP
#define HEP_MC_MULTI_CHANNEL_HPP

/*
 * hep-mc - A Template Library for Monte Carlo Integration
 * Copyright (C) 2015  Christopher Schwan
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

#include "hep/mc/multi_channel_point.hpp"
#include "hep/mc/multi_channel_result.hpp"

#include <cmath>
#include <cstddef>
#include <limits>
#include <random>
#include <vector>

namespace hep
{

/// \addtogroup multi_channel_group
/// @{

/// Multi-channel integrator. Integrates `function` over the unit-hypercube with
/// the specified number of dimensions using `calls` number of function calls.
/// The integration is performed using random numbers from the specified
/// generator and the given `densities`. This must be a function with the
/// following signature:
/// \code
/// void my_density_functions(
///     std::size_t channel,
///     std::vector<T> const& point,
///     std::vector<T>& channel_densities
/// );
/// \endcode
/// The variable `point` gives the uniformly generated point in the
/// unit-hypercube, `channel` is the number of the channel that was reandomly
/// selected and  for the selected `channel`, and `channel_densities` is the
/// vector where the result of the densities for this combination of `point` and
/// `channel` must be stored.
template <typename T, typename F, typename D, typename R>
inline multi_channel_result<T> multi_channel_iteration(
	std::size_t dimensions,
	std::size_t calls,
	F&& function,
	std::vector<T> const& channel_weights,
	D&& densities,
	R&& generator
) {
	T sum = T();
	T sum_of_squares = T();

	// for kahan summation
	T compensation = T();

	std::size_t const channels = channel_weights.size();

	std::vector<T> random_numbers(dimensions);
	std::vector<T> channel_densities(channels);
	std::vector<T> adjustment_data(channels);

	// distribution that randomly selects a channel
	std::discrete_distribution<std::size_t> channel_selector(
		channel_weights.begin(), channel_weights.end());

	for (std::size_t i = 0; i != calls; ++i)
	{
		// generate as many random numbers as we need
		for (std::size_t j = 0; j != dimensions; ++j)
		{
			random_numbers[j] = std::generate_canonical<T,
				std::numeric_limits<T>::digits>(generator);
		}

		// randomly select a channel
		std::size_t const channel = channel_selector(generator);

		// compute the densities for `random_numbers` for every channel
		densities(channel, static_cast <std::vector<T> const> (random_numbers),
			channel_densities);

		T total_density = T();
		for (std::size_t j = 0; j != channels; ++j)
		{
			total_density += channel_weights[j] * channel_densities[j];
		}

		multi_channel_point<T> const point(calls, random_numbers, channel,
			total_density);

		T const value = function(point) * point.weight;

		// perform kahan summation
		T const y = value - compensation;
		T const t = sum + y;
		compensation = (t - sum) - y;
		sum = t;

		T const square = value * value;

		sum_of_squares += square;

		// these are the values W that are used to update the alphas
		for (std::size_t j = 0; j != channels; ++j)
		{
			adjustment_data[j] += channel_densities[j] * square /
				total_density;
		}
	}

	sum *= T(calls);
	sum_of_squares *= T(calls) * T(calls);

	return multi_channel_result<T>(calls, sum, sum_of_squares, adjustment_data,
		channel_weights);
}

/// Uses `adjustment_data` to from a previous call of
/// \ref multi_channel_iteration to refine `weights`.
template <typename T>
inline std::vector<T> multi_channel_refine_weights(
	std::vector<T> const& weights,
	std::vector<T> const& adjustment_data
) {
	std::vector<T> new_weights(weights.size());

	T sum_of_new_weights = T();

	for (std::size_t i = 0; i != new_weights.size(); ++i)
	{
		new_weights[i] = weights[i] * std::sqrt(adjustment_data[i]);
		sum_of_new_weights += new_weights[i];
	}

	for (std::size_t i = 0; i != new_weights.size(); ++i)
	{
		new_weights[i] /= sum_of_new_weights;
	}

	return new_weights;
}

/// Performs `iteration_calls.size()` multi channel iterations by calling
/// \ref multi_channel_iteration with the specified parameters and refining
/// the weights after each iteration with \ref multi_channel_refine_weights.
template <typename T, typename F, typename D, typename R = std::mt19937>
inline std::vector<multi_channel_result<T>> multi_channel(
	std::size_t dimensions,
	std::vector<std::size_t> iteration_calls,
	F&& function,
	std::size_t channels,
	D&& densities,
	R&& generator = std::mt19937()
) {
	// start with an equal probability for every channel
	std::vector<T> weights(channels, T(1.0) / T(channels));

	std::vector<multi_channel_result<T>> results;
	results.reserve(iteration_calls.size());

	for (std::size_t const calls : iteration_calls)
	{
		auto const result = multi_channel_iteration(dimensions, calls, function,
			weights, densities, generator);
		results.push_back(result);

		weights = multi_channel_refine_weights(weights, result.adjustment_data);
	}

	return results;
}

/// @}

}

#endif
