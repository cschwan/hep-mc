#ifndef HEP_MC_IMPL_PLAIN_HPP
#define HEP_MC_IMPL_PLAIN_HPP

/*
 * hep-ga - An Efficient Numeric Template Library for Geometric Algebra
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

#include <hep/mc/plain.hpp>

#include <vector>

namespace hep
{

template <typename T, typename F, typename A, typename R>
plain_result<T> plain(
	std::size_t dimension,
	std::size_t steps,
	std::size_t seed,
	F& function,
	A const& auxilliary_variable,
	R&& random_number_generator
) {
	// default-initialize sum and sum_of_squares
	T sum = T();
	T sum_of_squares = T();

	// container holding random numbers
	std::vector<T> random_numbers(dimension);

	// distribution [0, 1] for the random number generator
	std::uniform_real_distribution<T> distribution;

	// seed random number generator
	random_number_generator.seed(seed);

	// iterate over samples
	for (std::size_t i = 0; i != steps; ++i)
	{
		// fill container with random numbers
		for (std::size_t j = 0; j != dimension; ++j)
		{
			random_numbers[j] = distribution(random_number_generator);
		}

		// evaluate function at position specified in random_numbers
		T evaluation = function(random_numbers, auxilliary_variable);

		sum += evaluation;
		sum_of_squares += evaluation * evaluation;
	}

	// average sum
	T value = sum / T(steps);

	// compute error for value
	T error = std::sqrt(sum_of_squares / (T(steps) * T(steps)) 
		- value * value / T(steps));

	return { value, error };
}

}

#endif
