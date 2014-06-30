#ifndef HEP_MC_VEGAS_ITERATION_RESULT_HPP
#define HEP_MC_VEGAS_ITERATION_RESULT_HPP

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

#include "hep/mc/mc_result.hpp"
#include "hep/mc/vegas_pdf.hpp"

#include <cstddef>
#include <vector>

namespace hep
{

/// \addtogroup vegas
/// @{

/// The result of a single \ref vegas_iteration.
template <typename T>
struct vegas_iteration_result : public mc_result<T>
{
	/// Constructor.
	vegas_iteration_result(
		std::size_t calls,
		vegas_pdf<T> const& pdf,
		std::vector<T> const& adjustment_data
	)
		: mc_result<T>(calls, *(adjustment_data.end() - 2), *(adjustment_data.end() - 1))
		, pdf(pdf)
		, adjustment_data(adjustment_data)
	{
	}

	/// The pdf used to obtain this result.
	vegas_pdf<T> pdf;

	/// The data used to adjust the `pdf` for a subsequent iteration.
	std::vector<T> adjustment_data;
};

/// @}

}

#endif
