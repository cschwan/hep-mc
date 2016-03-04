#ifndef HEP_MC_VEGAS_RESULT_HPP
#define HEP_MC_VEGAS_RESULT_HPP

/*
 * hep-mc - A Template Library for Monte Carlo Integration
 * Copyright (C) 2014-2016  Christopher Schwan
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
#include "hep/mc/plain_result.hpp"
#include "hep/mc/vegas_pdf.hpp"

#include <cstddef>
#include <vector>

namespace hep
{

/// \addtogroup results
/// @{

/// The result of a single \ref vegas_iteration.
template <typename T>
class vegas_result : public plain_result<T>
{
public:
	/// Constructor.
	vegas_result(
		std::vector<distribution_result<T>> const& distributions,
		std::size_t calls,
		T sum,
		T sum_of_squares,
		vegas_pdf<T> const& pdf,
		std::vector<T> const& adjustment_data
	)
		: plain_result<T>(distributions, calls, sum, sum_of_squares)
		, pdf_(pdf)
		, adjustment_data_(adjustment_data)
	{
	}

	/// The pdf used to obtain this result.
	vegas_pdf<T> const& pdf() const
	{
		return pdf_;
	}

	/// The data used to adjust the \ref pdf for a subsequent iteration.
	std::vector<T> const& adjustment_data() const
	{
		return adjustment_data_;
	}

private:
	vegas_pdf<T> pdf_;
	std::vector<T> adjustment_data_;
};

/// @}

}

#endif
