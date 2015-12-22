#ifndef HEP_MC_MULTI_CHANNEL_DISTRIBUTION_RESULT_HPP
#define HEP_MC_MULTI_CHANNEL_DISTRIBUTION_RESULT_HPP

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

#include "hep/mc/distribution_result.hpp"
#include "hep/mc/multi_channel_result.hpp"

namespace hep
{

/// \addtogroup results
/// @{

/// Captures the result of distributions and the total Multi Channel integration
/// result.
template <typename T>
using multi_channel_distribution_result =
	distribution_result<multi_channel_result<T>>;


/// @}

}

#endif

