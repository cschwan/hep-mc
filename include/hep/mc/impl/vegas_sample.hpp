#ifndef HEP_MC_IMPL_VEGAS_SAMPLE_HPP
#define HEP_MC_IMPL_VEGAS_SAMPLE_HPP

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

#include <hep/mc/vegas_sample.hpp>

namespace hep
{

template <typename T>
vegas_sample<T>::vegas_sample(
	std::size_t dimensions,
	std::size_t steps,
	T const& random,
	std::vector<std::vector<T>> const& grid
)
	: point(dimensions)
	, bin(dimensions)
	, weight(1.0 / steps)
{
	std::size_t const bins = grid[0].size();

	for (std::size_t i = 0; i != dimensions; ++i)
	{
		// compute position of 'random' in bins, as a floating point number
		T const bin_position = random * bins;

		// in which bin is 'random' (integer) ?
		std::size_t const position = bin_position;

		// compute value of grid at the previous position
		T const grid_previous = (position == 0) ? T() : grid[i][position - 1];

		// compute difference of grid values at 'position'
		T const difference = grid[i][position] - grid_previous;

		// TODO: explain
		point[i] = grid_previous + (bin_position - position) * difference;

		// save the index of the bin in which point lies
		bin[i] = position;

		// multiply weight for each dimension
		weight *= difference * bins;
	}
}

}

#endif
