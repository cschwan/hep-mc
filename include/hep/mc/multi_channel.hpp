#ifndef HEP_MC_MULTI_CHANNEL_HPP
#define HEP_MC_MULTI_CHANNEL_HPP

/*
 * hep-mc - A Template Library for Monte Carlo Integration
 * Copyright (C) 2015-2016  Christopher Schwan
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

#include "hep/mc/internal/accumulator.hpp"
#include "hep/mc/internal/discrete_distribution.hpp"
#include "hep/mc/integrand.hpp"
#include "hep/mc/multi_channel_callback.hpp"
#include "hep/mc/multi_channel_map.hpp"
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

/// Uses `adjustment_data` from a previous call of \ref multi_channel_iteration
/// to refine `weights`. The procedure is the one suggested in Ref.
/// \cite WeightOptimization with the following modifications:
/// - the weights are check if they are smaller then the given
///   `minimum_weight`. If this is the case they are set to the value of
///   `minimum_weight` which, after the normalization of all weights, make
///   them a little smaller then the given minimum weight,
/// - `adjustment_data` is raised to the power given by `beta`. The
///   reference given above suggest `beta = 0.5`, but the default value is
///   smaller which sometimes gives a more stable convergence.
template <typename T>
inline std::vector<T> multi_channel_refine_weights(
	std::vector<T> const& weights,
	std::vector<T> const& adjustment_data,
	T minimum_weight,
	T beta = T(0.25)
) {
	std::vector<T> new_weights(weights.size());

	T sum_of_new_weights = T();

	for (std::size_t i = 0; i != new_weights.size(); ++i)
	{
		new_weights[i] = weights[i] * std::pow(adjustment_data[i], beta);
		sum_of_new_weights += new_weights[i];
	}

	T new_sum = T();

	for (T& weight : new_weights)
	{
		if (weight == T())
		{
			// do not enable disabled channels (with weight zero) by setting
			// them to the minimum weight
			continue;
		}

		weight /= sum_of_new_weights;
		weight = std::fmax(weight, minimum_weight);
		new_sum += weight;
	}

	for (T& weight : new_weights)
	{
		weight /= new_sum;
	}

	return new_weights;
}

/// Performs exactly one iteration of the multi channel integration. The
/// parameter `channel_weights` lets the user specify the weights of each
/// channel. Note that they must add up to one. See \ref multi_channel_group for
/// a description of the remaining parameters.
template <typename I, typename R>
inline multi_channel_result<numeric_type_of<I>> multi_channel_iteration(
	I&& integrand,
	std::size_t calls,
	std::vector<numeric_type_of<I>> const& channel_weights,
	R&& generator
) {
	using T = numeric_type_of<I>;

	auto accumulator = make_accumulator(integrand);

	std::size_t const channels = channel_weights.size();

	std::vector<T> random_numbers(integrand.dimensions());
	std::vector<T> coordinates(integrand.map_dimensions());
	std::vector<T> densities(channels);
	std::vector<T> adjustment_data(channels);

	// distribution that randomly selects a channel
	discrete_distribution<std::size_t, T> channel_selector(
		channel_weights.begin(), channel_weights.end());

	for (std::size_t i = 0; i != calls; ++i)
	{
		// generate as many random numbers as we need
		for (std::size_t j = 0; j != integrand.dimensions(); ++j)
		{
			random_numbers[j] = std::generate_canonical<T,
				std::numeric_limits<T>::digits>(generator);
		}

		// randomly select a channel
		std::size_t const channel = channel_selector(generator);

		using map_type = typename std::remove_reference<
			typename std::remove_reference<I>::type::map_type>::type;

		// calculate `coordinates` and possibly `densities`
		integrand.map()(
			channel,
			random_numbers,
			coordinates,
			densities,
			multi_channel_map::calculate_coordinates
		);

		multi_channel_point2<T, map_type> const point(
			random_numbers,
			coordinates,
			channel,
			densities,
			channel_weights,
			integrand.map()
		);

		T const value = accumulator.invoke(integrand, point);

		if (value == T())
		{
			continue;
		}

		T const square = value * value * point.weight();

		// these are the values W that are used to update the alphas
		for (std::size_t j = 0; j != channels; ++j)
		{
			adjustment_data[j] += densities[j] * square;
		}
	}

	return multi_channel_result<T>(
		accumulator.distributions(calls),
		calls,
		accumulator.sum(),
		accumulator.sum_of_squares(),
		adjustment_data,
		channel_weights
	);
}

/// Performs `iteration_calls.size()` multi channel iterations by calling
/// \ref multi_channel_iteration with the specified parameters and refining
/// the weights after each iteration with \ref multi_channel_refine_weights. The
/// weights that are used for the first iteration must be given by the parameter
/// `channel_weights`. See \ref multi_channel_group for a description of the
/// remaining parameters.
template <typename I, typename R = std::mt19937>
inline std::vector<multi_channel_result<numeric_type_of<I>>> multi_channel(
	I&& integrand,
	std::vector<std::size_t> const& iteration_calls,
	std::vector<numeric_type_of<I>> const& channel_weights,
	std::size_t min_calls_per_channel = 0,
	R&& generator = std::mt19937()
) {
	using T = numeric_type_of<I>;

	auto weights = channel_weights;

	std::vector<multi_channel_result<T>> results;
	results.reserve(iteration_calls.size());

	for (std::size_t const calls : iteration_calls)
	{
		auto const result = multi_channel_iteration(integrand, calls, weights,
			generator);

		results.push_back(result);

		if (!multi_channel_callback<T>()(results))
		{
			break;
		}

		T const minimum_weight = T(min_calls_per_channel) / calls;

		weights = multi_channel_refine_weights(weights,
			result.adjustment_data(), minimum_weight);
	}

	return results;
}

/// Performs `iteration_calls.size()` multi channel iterations by calling
/// \ref multi_channel_iteration with the specified parameters and refining
/// the weights after each iteration with \ref multi_channel_refine_weights. The
/// weights used for the first iteration are \f$ \alpha = 1 / M \f$ with \f$ M
/// \f$ the number of channels. See \ref multi_channel_group for a description
/// of the remaining parameters.
template <typename I, typename R = std::mt19937>
inline std::vector<multi_channel_result<numeric_type_of<I>>> multi_channel(
	I&& integrand,
	std::vector<std::size_t> const& iteration_calls,
	std::size_t min_calls_per_channel = 0,
	R&& generator = std::mt19937()
) {
	using T = numeric_type_of<I>;

	return multi_channel(
		std::forward<I>(integrand),
		iteration_calls,
		std::vector<T>(integrand.channels(), T(1.0) / T(integrand.channels())),
		min_calls_per_channel,
		std::forward<R>(generator)
	);
}

/// @}

}

#endif
