#ifndef HEP_MC_VEGAS_SAMPLE_HPP
#define HEP_MC_VEGAS_SAMPLE_HPP

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

#include <cstddef>
#include <vector>

namespace hep
{

/**
 * 
 */
template <typename T>
struct vegas_sample
{
	/**
	 * 
	 */
	vegas_sample(
		std::size_t dimensions,
		std::size_t steps,
		T const& random,
		std::vector<std::vector<T>> const& grid
	);

	/**
	 * 
	 */
	std::vector<T> point;

	/**
	 * 
	 */
	std::vector<std::size_t> bin;

	/**
	 * 
	 */
	T weight;
};

}

#include <hep/mc/impl/vegas_sample.hpp>

#endif
