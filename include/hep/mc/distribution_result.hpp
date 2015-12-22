#ifndef HEP_MC_DISTRIBUTION_RESULT_HPP
#define HEP_MC_DISTRIBUTION_RESULT_HPP

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
class distribution_result : public R
{
public:
	/// Shortcut to the numeric type of `R`.
	using numeric_type = typename R::numeric_type;

	/// Constructor. This stores `distribution_results` and forwards all
	/// remaining arguments to the constructor of `R`.
	template <typename... A>
	distribution_result(
		std::vector<std::vector<mc_result<numeric_type>>> const&
			distribution_results,
		std::size_t calls,
		numeric_type sum,
		numeric_type sum_of_squares,
		A&&... args
	)
		: R(calls, sum, sum_of_squares, std::forward<A>(args)...)
		, distribution_results_(distribution_results)
	{
	}

	/// Returns the results for each bin in every distribution. The bin with
	/// index `j` of distribution `i` is accessed as follows:
	/// `distribution_results().at(i).at(j)`.
	std::vector<std::vector<mc_result<numeric_type>>> distribution_results()
	const
	{
		return distribution_results_;
	}

private:
	std::vector<std::vector<mc_result<numeric_type>>> distribution_results_;
};

/// @}

}

#endif
