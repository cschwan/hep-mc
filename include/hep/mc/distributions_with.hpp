#ifndef HEP_MC_DISTRIBUTIONS_WITH_HPP
#define HEP_MC_DISTRIBUTIONS_WITH_HPP

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

#include "hep/mc/mc_result.hpp"

#include <cstddef>
#include <utility>
#include <vector>

namespace hep
{

/// \addtogroup results
/// @{

/// Holds the result of type `R` for an integration and, in addition, the
/// results for the generated distributions.
template <typename R>
class distributions_with
{
public:
	/// Shortcut to the numeric type of `R`.
	using numeric_type = typename R::numeric_type;

	/// Constructor. This stores `distributions` and forwards all remaining
	/// arguments to the constructor of `R`.
	template <typename... A>
	distributions_with(
		std::vector<distribution_result<numeric_type>> const& distributions,
		std::size_t calls,
		numeric_type sum,
		numeric_type sum_of_squares,
		A&&... args
	)
		: integral_(calls, sum, sum_of_squares, std::forward<A>(args)...)
		, distributions_(distributions)
	{
	}

	/// Returns the results for every distribution.
	std::vector<distribution_result<numeric_type>> distributions() const
	{
		return distributions_;
	}

	/// Returns the result of the total integration.
	R integral() const
	{
		return integral_;
	}

private:
	R integral_;
	std::vector<distribution_result<numeric_type>> distributions_;
};

/// @}

}

#endif
