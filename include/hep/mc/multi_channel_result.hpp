#ifndef HEP_MC_MULTI_CHANNEL_RESULT_HPP
#define HEP_MC_MULTI_CHANNEL_RESULT_HPP

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

#include "hep/mc/plain_result.hpp"

#include <cstddef>
#include <vector>

namespace hep
{

/// \addtogroup results
/// @{

/// Result of a multi-channel integration.
template <typename T>
class multi_channel_result : public plain_result<T>
{
public:
	/// Constructor. The parameter `adjustment_data` must contain two additional
	/// values, being the sum of the sum of squares of the integration.
	multi_channel_result(
		std::vector<distribution_result<T>> const& distributions,
		std::size_t calls,
		T sum,
		T sum_of_squares,
		std::vector<T> const& adjustment_data,
		std::vector<T> const& channel_weights
	)
		: plain_result<T>(distributions, calls, sum, sum_of_squares)
		, adjustment_data_(adjustment_data)
		, channel_weights_(channel_weights)
	{
	}

	/// This is the data used by \ref multi_channel_refine_weights to refine the
	/// \ref channel_weights used in the same iteration. The refined weights are
	/// then used in a subsequent iteration.
	std::vector<T> const& adjustment_data() const
	{
		return adjustment_data_;
	}

	/// The weight for each channel.
	std::vector<T> const& channel_weights() const
	{
		return channel_weights_;
	}

private:
	std::vector<T> adjustment_data_;
	std::vector<T> channel_weights_;
};

/// @}

}

#endif
