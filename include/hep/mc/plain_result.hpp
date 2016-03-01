#ifndef HEP_MC_PLAIN_RESULT_HPP
#define HEP_MC_PLAIN_RESULT_HPP

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

#include "hep/mc/distribution_result.hpp"
#include "hep/mc/mc_result.hpp"

#include <cstddef>
#include <vector>

namespace hep
{

/// \addtogroup results
/// @{

/// Return type of the \ref plain MC integrator.
template <typename T>
class plain_result : public mc_result<T>
{
public:
	/// Constructor.
	plain_result(
		std::vector<distribution_result<T>> const& distributions,
		std::size_t calls,
		T sum,
		T sum_of_squares
	)
		: mc_result<T>(calls, sum, sum_of_squares)
		, distributions_(distributions)
	{
	}

	/// Returns the differential distributions accumulated during the
	/// integration.
	std::vector<distribution_result<T>> const& distributions() const
	{
		return distributions_;
	}

private:
	std::vector<distribution_result<T>> distributions_;
};

/// @}

}

#endif

