#ifndef HEP_MC_MC_POINT_HPP
#define HEP_MC_MC_POINT_HPP

/*
 * hep-mc - A Template Library for Monte Carlo Integration
 * Copyright (C) 2012-2013  Christopher Schwan
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

/// \addtogroup integrands
/// @{

/**
 * A random point \f$ \vec{x} \in [0,1]^d \f$ in the \f$ d \f$-dimensional hypercube with an
 * associated weight.
 */
template <typename T>
struct mc_point
{
	/**
	 * Constructor. The weight is computed using the inverse of `calls`, `point` determines the
	 * location of this point in the hypercube of dimension `point.size()`.
	 */
	mc_point(std::size_t calls, std::vector<T> const& point)
		: weight(T(1.0) / T(calls))
		, point(point)
	{
	}

	/**
	 * The weight \f$ w \f$ of this point. Depending on the integration algorithm the weight might
	 * be constant over the entire hypercube (e.g. \ref plain) or dependent on the region in which
	 * it lies (e.g. \ref vegas).
	 */
	T weight;

	/**
	 * The coordinates of this point of the hypercube. The dimension can be obtained using
	 * `point.size()`.
	 */
	std::vector<T> const& point;
};

/// @}

}

#endif
