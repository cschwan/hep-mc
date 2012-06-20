#ifndef HEP_MC_VEGAS_HPP
#define HEP_MC_VEGAS_HPP

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

#include <cstddef>
#include <random>
#include <vector>

namespace hep
{

/**
 * VEGAS Monte Carlo integrator. \c vegas integrates the specified \c function
 * with the specified \c dimensions using <tt>steps.size()</tt> iterations
 * each  having the number integrand evaluation specified in \c steps. The grid
 * has a resolution specified with \c bins in each dimension.
 */
template <typename T, typename F, typename A, typename R = std::mt19937>
vegas_result<T> vegas(
	std::size_t dimensions,
	std::vector<std::size_t> steps,
	std::size_t batch_size,
	std::size_t bins,
	F& function,
	A const& aux_variable,
	std::size_t seed = 0,
	R&& generator = std::mt19937()
);

}

#include <hep/mc/impl/vegas.hpp>

#endif
