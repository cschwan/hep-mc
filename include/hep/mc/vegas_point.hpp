#ifndef HEP_MC_VEGAS_POINT_HPP
#define HEP_MC_VEGAS_POINT_HPP

/*
 * hep-mc - A Template Library for Monte Carlo Integration
 * Copyright (C) 2014  Christopher Schwan
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

#include "hep/mc/mc_point.hpp"
#include "hep/mc/vegas_pdf.hpp"

#include <cstddef>
#include <vector>

namespace hep
{

/// \addtogroup vegas
/// @{

/// A point within the unit-hypercube with the additional information in which bins the point lies.
template <typename T>
struct vegas_point : public mc_point<T>
{
	/**
	 * Constructor. Creates a new point the unit-hypercube using the probability distribution
	 * function `pdf` and the random numbers in `random_numbers` for a Monte Carlo iteration with
	 * sample size specified with `total_calls`. For each dimension the point falls into a bin with
	 * indices given by `bin`.
	 */
	vegas_point(
		std::size_t total_calls,
		std::vector<T>& random_numbers,
		std::vector<std::size_t>& bin,
		vegas_pdf<T> const& pdf
	)
		: mc_point<T>(total_calls, random_numbers)
		, bin(bin)
	{
		this->weight *= pdf.icdf(random_numbers, bin);
	}

	/// The indices that determine the bins of the point in the binned pdf.
	std::vector<std::size_t> const& bin;
};

/// @}

}

#endif
