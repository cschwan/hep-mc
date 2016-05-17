#ifndef HEP_MC_DISTRIBUTION_RESULT_HPP
#define HEP_MC_DISTRIBUTION_RESULT_HPP

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

#include "hep/mc/distribution_parameters.hpp"
#include "hep/mc/mc_result.hpp"

#include <vector>

namespace hep
{

/// \addtogroup distributions
/// @{

/// Captures the result of the integration of a distribution.
template <typename T>
class distribution_result
{
public:
	/// Constructor.
	distribution_result(
		distribution_parameters<T> const& parameters,
		std::vector<mc_result<T>> const& results
	)
		: parameters_(parameters)
		, results_(results)
	{
	}

	/// Returns the parameters associated with this distribution.
	distribution_parameters<T> const& parameters() const
	{
		return parameters_;
	}

	/// Returns the middle point of each bin of this distribution.
	std::vector<T> mid_points() const
	{
		std::vector<T> result;
		result.reserve(parameters_.bins());

		T x = parameters_.x_min() + T(0.5) * parameters_.bin_size();

		for (std::size_t bin = 0; bin != parameters_.bins(); ++bin)
		{
			result.push_back(x);
			x += parameters_.bin_size();
		}

		return result;
	}

	/// Returns the result for each bin, corresponding to the bin positions
	/// returned by \ref mid_points.
	std::vector<mc_result<T>> const& results() const
	{
		return results_;
	}

private:
	distribution_parameters<T> parameters_;
	std::vector<mc_result<T>> results_;
};

/// @}

}

#endif
