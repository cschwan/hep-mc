#ifndef HEP_MC_MULTI_CHANNEL_WEIGHT_INFO_HPP
#define HEP_MC_MULTI_CHANNEL_WEIGHT_INFO_HPP

/*
 * hep-mc - A Template Library for Monte Carlo Integration
 * Copyright (C) 2016  Christopher Schwan
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

#include <algorithm>
#include <cstddef>
#include <iterator>
#include <vector>

namespace hep
{

/// \addtogroup results
/// @{

/// Class that displays the available information about the a-posterioi weights
/// of a \ref multi_channel_result.
template <typename T>
class multi_channel_weight_info
{
public:
	/// Constructor.
	multi_channel_weight_info(multi_channel_result<T> const& result)
		: channels_(result.channel_weights().size())
		, weights_(result.channel_weights().size())
		, calls_(result.channel_weights().size())
		, minimal_weight_count_()
	{
		std::iota(channels_.begin(), channels_.end(), 0);
		std::stable_sort(channels_.begin(), channels_.end(),
			[&](std::size_t a, std::size_t b) {
				return result.channel_weights().at(a) <
					result.channel_weights().at(b);
		});

		std::transform(channels_.begin(), channels_.end(), weights_.begin(),
			[&](std::size_t index) {
				return result.channel_weights().at(index);
		});

		std::transform(weights_.begin(), weights_.end(), calls_.begin(),
			[&](T weight) {
				return static_cast <std::size_t> (result.calls() * weight);
		});

		minimal_weight_count_ = std::distance(calls_.begin(),
			std::upper_bound(calls_.begin(), calls_.end(), calls_.front()));
	}

	/// Returns the number of expected calls for each weight in the same order
	/// as the channel indices returned by \ref channels().
	std::vector<std::size_t> const& calls() const
	{
		return calls_;
	}

	/// Returns the channel indices sorted in ascending order of their
	/// corresponding weight.
	std::vector<std::size_t> const& channels() const
	{
		return channels_;
	}

	/// Returns the number of channels that have a weight that, multiplied with
	/// the number of calls, corresponds to the minimum number of calls.
	std::size_t minimal_weight_count() const
	{
		return minimal_weight_count_;
	}

	/// Returns the a-priori weights, starting with the lowest weights up to the
	/// largest one.
	std::vector<T> const& weights() const
	{
		return weights_;
	}

private:
	std::vector<std::size_t> channels_;
	std::vector<T> weights_;
	std::vector<std::size_t> calls_;
	std::size_t minimal_weight_count_;
};

/// Returns a vector containing the indices of the channels with the smallest
/// weights.
template <typename T>
inline std::vector<std::size_t> minimal_weight_channels(
	multi_channel_weight_info<T> const& info
) {
	std::vector<std::size_t> result(
		info.channels().begin(),
		info.channels().begin() + info.minimal_weight_count()
	);

	return result;
}

/// @}

}

#endif
