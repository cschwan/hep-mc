#ifndef HEP_MC_PLAIN_HPP
#define HEP_MC_PLAIN_HPP

/*
 * hep-mc - A Template Library for Monte Carlo Integration
 * Copyright (C) 2012-2015  Christopher Schwan
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

#include "hep/mc/distribution_accumulator.hpp"
#include "hep/mc/distributions.hpp"
#include "hep/mc/mc_point.hpp"
#include "hep/mc/mc_result.hpp"

#include <cstddef>
#include <limits>
#include <random>
#include <utility>
#include <vector>

namespace
{

template <typename T, typename F, typename D, typename R>
inline hep::distribution_result<T> plain_iteration(
	std::size_t dimensions,
	std::size_t calls,
	F&& function,
	D&& distributions,
	R&& generator
) {
	auto accumulator = make_distribution_accumulator(distributions);

	// storage for random numbers
	std::vector<T> random_numbers(dimensions);

	// iterate over calls
	for (std::size_t i = 0; i != calls; ++i)
	{
		// fill container with random numbers
		for (std::size_t j = 0; j != dimensions; ++j)
		{
			random_numbers[j] = std::generate_canonical<T,
				std::numeric_limits<T>::digits>(generator);
		}

		hep::mc_point<T> const point(random_numbers);

		// evaluate function at position specified in random_numbers
		T const value = function(point);

		accumulator.add(point, value);
	}

	return accumulator.result();
}

template <typename T, typename F, typename R>
inline hep::mc_result<T> plain_iteration(
	std::size_t dimensions,
	std::size_t calls,
	F&& function,
	R&& generator
) {
	// use default_distribution, which has one distribution with exactly one
	// bin, i.e. every point lands in there -> simple MC integration
	return plain_iteration<T>(
		dimensions,
		calls,
		std::forward<F>(function),
		hep::default_distribution<T>(),
		std::forward<R>(generator)
	);
}

}

namespace hep
{

/// \addtogroup plain_group
/// @{

/// PLAIN Monte Carlo integrator. This function integrates `function` over the
/// unit-hypercube with the specified `dimensions` using `calls` function
/// evaluations with randomly chosen points determined by `generator`. The
/// generator is not seeded.

/// \param dimensions The number of parameters `function` accepts.
/// \param calls The number of function calls that are used to obtain the
///        result.
/// \param function The function that will be integrated over the hypercube. See
///        \ref integrands for further explanation.
/// \param generator The random number generator that will be used to generate
///        random points from the hypercube. This generator is not seeded.
template <typename T, typename F, typename R = std::mt19937>
inline mc_result<T> plain(
	std::size_t dimensions,
	std::size_t calls,
	F&& function,
	R&& generator = std::mt19937()
) {
	return plain_iteration<T>(
		dimensions,
		calls,
		std::forward<F>(function),
		std::forward<R>(generator)
	);
}

/// @}

}

#endif
