#ifndef HEP_MC_MC_POINT_HPP
#define HEP_MC_MC_POINT_HPP

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
 * Instances of this class represent Monte Carlo points, i.e. a random point
 * \f$ \vec{x} \in [0,1]^d \f$ from the hypercube in \f$ d \f$ dimensions
 * together with a \ref weight.
 */
template <typename T>
struct mc_point
{
	/**
	 * Constructor.
	 */
	mc_point(std::size_t samples, std::vector<T>& point)
		: weight(T(1.0) / T(samples))
		, point(point)
	{
	}

	/**
	 * The weight \f$ w \f$ of this point. Depending on the integration
	 * algorithm the weight might be constant over the entire hypercube (e.g.
	 * \ref plain) or dependent of the region in which it lies (e.g.
	 * \ref vegas).
	 */
	T weight;

	/**
	 * The coordinates of this point. The size of this vector determines the
	 * dimension \f$ d \f$ of the hypercube, the entries with index \f$ i \f$
	 * the coordinates \f$ x_i \f$ so that \f$ \vec{x} = \left( x_0, x_1,
	 * \ldots, x_{d-1} \right) \f$.
	 */
	std::vector<T>& point;
};

}

#endif
