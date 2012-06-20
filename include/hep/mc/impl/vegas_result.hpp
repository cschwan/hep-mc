#ifndef HEP_MC_IMPL_VEGAS_RESULT_HPP
#define HEP_MC_IMPL_VEGAS_RESULT_HPP

/*
 * hep-mc - A Template Library for Monte Carlo Integration
 * Copyright (C) 2012  Christopher Schwan
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

#include <hep/mc/vegas_result.hpp>

#include <algorithm>
#include <cmath>

namespace hep
{

template <typename T>
vegas_result<T>::vegas_result(std::vector<std::size_t> const& steps)
	: values()
	, errors()
	, steps(steps)
	, sum_of_inv_variances()
	, sum_of_averages()
{
	values.reserve(steps.size());
	errors.reserve(steps.size());
}

template <typename T>
void vegas_result<T>::add_iteration(
	std::size_t samples,
	T const& sum,
	T const& sum_of_squares
) {
	T const tmp = std::sqrt(sum_of_squares * samples);

	// compute inverse variance
	T const inv_variance = T(samples - 1.0) /
		std::max((tmp + sum) * (tmp - sum), T(0x1p-104));

	sum_of_inv_variances += inv_variance;
	T const variance = T(1.0) / sum_of_inv_variances;
	sum_of_averages += inv_variance * sum;
	T const average = variance * sum_of_averages;

	values.push_back(average);
	errors.push_back(std::sqrt(variance));
}

}

#endif
