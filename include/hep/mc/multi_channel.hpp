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

#include "hep/mc/kahan_accumulator.hpp"
#include "hep/mc/multi_channel_callback.hpp"
#include "hep/mc/multi_channel_point.hpp"
#include "hep/mc/multi_channel_result.hpp"

#include <cmath>
#include <cstddef>
#include <limits>
#include <random>
#include <type_traits>
#include <utility>
#include <vector>

namespace hep
{

/// \addtogroup multi_channel_group
/// @{

/// Uses `adjustment_data` to from a previous call of
/// \ref multi_channel_iteration to refine `weights`. The procedure is the one
/// suggested in Ref. \cite WeightOptimization .
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

/// Performs exactly one iteration of the multi channel integration. The
/// parameter `channel_weights` lets the user specify the weights of each
/// channel. Note that they must add up to one. See \ref multi_channel_group for
/// a description of the remaining parameters.
template <typename T, typename F, typename D, typename R>
inline multi_channel_result<T> multi_channel_iteration(
	std::size_t dimensions,
	std::size_t map_dimensions,
	std::size_t calls,
	std::size_t total_calls,
	F&& function,
	std::vector<T> const& channel_weights,
	D&& densities,
	R&& generator
) {
	kahan_accumulator<T> accumulator;

	std::size_t const channels = channel_weights.size();

	std::vector<T> random_numbers(dimensions);
	std::vector<T> coordinates(map_dimensions);
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
		T const weight = densities(channel, static_cast <std::vector<T> const&>
			(random_numbers), coordinates, channel_densities);

		if (weight == T())
		{
			// function does not contribute, skip evaluation
			continue;
		}

		T total_density = T();
		for (std::size_t j = 0; j != channels; ++j)
		{
			channel_densities[j] /= weight;
			total_density += channel_weights[j] * channel_densities[j];
		}

		multi_channel_point2<T, typename std::remove_reference<D>::type> const
			point(total_calls, random_numbers, coordinates, channel,
			total_density, densities);

		T const value = function(point) * point.weight * T(total_calls);

		accumulator.add(value);

		T const square = value * value;

		// these are the values W that are used to update the alphas
		for (std::size_t j = 0; j != channels; ++j)
		{
			adjustment_data[j] += channel_densities[j] * square / total_density;
		}
	}

	return multi_channel_result<T>(calls, accumulator.sum(),
		accumulator.sum_of_squares(), adjustment_data, channel_weights);
}

/// Performs `iteration_calls.size()` multi channel iterations by calling
/// \ref multi_channel_iteration with the specified parameters and refining
/// the weights after each iteration with \ref multi_channel_refine_weights. The
/// weights that are used for the first iteration must be given by the parameter
/// `channel_weights`. See \ref multi_channel_group for a description of the
/// remaining parameters.
template <typename T, typename F, typename D, typename R = std::mt19937>
inline std::vector<multi_channel_result<T>> multi_channel(
	std::size_t dimensions,
	std::size_t map_dimensions,
	std::vector<std::size_t> const& iteration_calls,
	F&& function,
	std::vector<T> const& channel_weights,
	D&& densities,
	R&& generator = std::mt19937()
) {
	auto weights = channel_weights;

	std::vector<multi_channel_result<T>> results;
	results.reserve(iteration_calls.size());

	for (std::size_t const calls : iteration_calls)
	{
		auto const result = multi_channel_iteration(
			dimensions,
			map_dimensions,
			calls,
			calls,
			function,
			weights,
			densities,
			generator
		);

		results.push_back(result);

		if (!multi_channel_callback<T>()(results))
		{
			break;
		}

		weights = multi_channel_refine_weights(weights,
			result.adjustment_data());
	}

	return results;
}

/// Performs `iteration_calls.size()` multi channel iterations by calling
/// \ref multi_channel_iteration with the specified parameters and refining
/// the weights after each iteration with \ref multi_channel_refine_weights. The
/// weights used for the first iteration are \f$ \alpha = 1 / M \f$ with \f$ M
/// \f$ the number of channels. See \ref multi_channel_group for a description
/// of the remaining parameters.
template <typename T, typename F, typename D, typename R = std::mt19937>
inline std::vector<multi_channel_result<T>> multi_channel(
	std::size_t dimensions,
	std::size_t map_dimensions,
	std::vector<std::size_t> const& iteration_calls,
	F&& function,
	std::size_t channels,
	D&& densities,
	R&& generator = std::mt19937()
) {
	return multi_channel(
		dimensions,
		map_dimensions,
		iteration_calls,
		std::forward<F>(function),
		std::vector<T>(channels, T(1.0) / T(channels)),
		std::forward<D>(densities),
		std::forward<R>(generator)
	);
}

/// @}

}

#endif
