#ifndef HEP_MC_PLAIN_HPP
#define HEP_MC_PLAIN_HPP

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

#include <hep/mc/mc_point.hpp>
#include <hep/mc/mc_result.hpp>

#include <cstddef>
#include <random>
#include <vector>

namespace hep
{

/**
 * \defgroup plain PLAIN Integrator
 * @{
 */

/**
 * PLAIN Monte Carlo integrator. This function integrates `function` over the
 * unit-hypercube with the specified `dimensions` using `calls` function
 * evaluations with randomly chosen points determined by `generator`. The
 * generator is not seeded.
 *
 * \param dimensions The number of parameters `function` accepts.
 * \param calls The number of function calls that are used to obtain the result.
 * \param function The function that will be integrated over the hypercube. See
 *        \ref integrands for further explanation.
 * \param generator The random number generator that will be used to generate
 *        random points from the hypercube. This generator is not seeded.
 */
template <typename T, typename F, typename R = std::mt19937>
mc_result<T> plain(
	std::size_t dimensions,
	std::size_t calls,
	F&& function,
	R&& generator = std::mt19937()
) {
	// default-initialize sum and sum_of_squares
	T sum = T();
	T sum_of_squares = T();

	// generates random number in the range [0,1]
	std::uniform_real_distribution<T> distribution;

	// compensation variable for kahan summation
	T compensation = T();

	std::vector<T> random_numbers(dimensions);

	// iterate over calls
	for (std::size_t i = 0; i != calls; ++i)
	{
		mc_point<T> point(calls, random_numbers);

		// fill container with random numbers
		for (std::size_t j = 0; j != dimensions; ++j)
		{
			random_numbers[j] = distribution(generator);
		}

		// evaluate function at position specified in random_numbers
		T const value = function(static_cast <mc_point<T> const> (point));

		// perform kahan summation 'sum += value' - this improves precision if
		// T is single precision and many values are added
		T const y = value - compensation;
		T const t = sum + y;
		compensation = (t - sum) - y;
		sum = t;

		sum_of_squares += value * value;
	}

	return mc_result<T>(calls, sum, sum_of_squares);
}

/**
 * @}
 */

}

#endif
