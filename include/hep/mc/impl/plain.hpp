#ifndef HEP_MC_IMPL_PLAIN_HPP
#define HEP_MC_IMPL_PLAIN_HPP

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

#include <hep/mc/plain.hpp>

#include <vector>

namespace hep
{

template <typename T>
plain_result<T>::plain_result(
	std::size_t steps,
	T const& sum,
	T const& sum_of_squares
)
	: steps(steps)
	, value(sum / T(steps))
	, error(std::sqrt(sum_of_squares / (T(steps) * T(steps)) 
	        - value * value / T(steps)))
{
}

template <typename T, typename F, typename A, typename R>
plain_result<T> plain(
	std::size_t dimensions,
	std::size_t steps,
	F& function,
	A const& auxilliary_variable,
	std::size_t seed,
	R&& generator
) {
	// default-initialize sum and sum_of_squares
	T sum = T();
	T sum_of_squares = T();

	// container holding random numbers
	std::vector<T> random_numbers(dimensions);

	// distribution [0, 1] for the random number generator
	std::uniform_real_distribution<T> distribution;

	// seed random number generator
	generator.seed(seed);

	// iterate over samples
	for (std::size_t i = 0; i != steps; ++i)
	{
		// fill container with random numbers
		for (std::size_t j = 0; j != dimensions; ++j)
		{
			random_numbers[j] = distribution(generator);
		}

		// evaluate function at position specified in random_numbers
		T evaluation = function(random_numbers, auxilliary_variable);

		sum += evaluation;
		sum_of_squares += evaluation * evaluation;
	}

	return plain_result<T>(steps, sum, sum_of_squares);
}

}

#endif
