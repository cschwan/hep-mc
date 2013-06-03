#ifndef HEP_MC_PLAIN_HPP
#define HEP_MC_PLAIN_HPP

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

#include <cmath>
#include <cstddef>
#include <functional>
#include <random>
#include <vector>

#include <hep/mc/mc_point.hpp>
#include <hep/mc/mc_result.hpp>

namespace hep
{

/**
 * PLAIN Monte Carlo integrator. \c plain integrates \c function over the
 * unit-hypercube with \c dimensions using \c steps randomly chosen points
 * determined by the number generator \c random_number_generator. The generator
 * is not seeded. \c function must have the following form:
 * \code
 * T integrand(std::vector<T> const& x)
 * {
 *     // return value of the function at x
 * }
 * \endcode
 * where \c T must be substituted with the type \c plain is called.
 *
 * The result of this function is an object of the type \ref plain_result
 * containing both the value of the integral and an error.
 */
template <typename T, typename F, typename R = std::mt19937>
mc_result<T> plain(
	std::size_t dimensions,
	std::size_t steps,
	F function,
	R generator = std::mt19937()
) {
	// default-initialize sum and sum_of_squares
	T sum = T();
	T sum_of_squares = T();

	mc_point<T> point(dimensions);
	point.weight = T(1.0) / T(steps);

	// distribution [0, 1] for the random number generator
	std::uniform_real_distribution<T> distribution;

	// compensation variable for kahan summation
	T compensation = T();

	// iterate over samples
	for (std::size_t i = 0; i != steps; ++i)
	{
		// fill container with random numbers
		for (std::size_t j = 0; j != dimensions; ++j)
		{
			point.point[j] = distribution(generator);
		}

		// evaluate function at position specified in random_numbers
		T const evaluation = function(std::cref(point));

		// do kahan summation
		T const y = evaluation - compensation;
		T const t = sum + y;
		compensation = (t - sum) - y;
		sum = t;

		sum_of_squares += evaluation * evaluation;
	}

	T const value = sum / T(steps);

	return mc_result<T>(steps, value, std::sqrt(sum_of_squares /
		(T(steps) * T(steps)) - value * value / T(steps)));
}

}

#endif
