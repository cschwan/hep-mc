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
#include <random>
#include <vector>

namespace hep
{

/**
 * The result of a PLAIN Monte Carlo integration.
 */
template <typename T>
struct plain_result
{
	/**
	 *
	 */
	plain_result(std::size_t steps, T const& sum, T const& sum_of_squares)
		: steps(steps)
		, value(sum / T(steps))
		, error(std::sqrt(sum_of_squares / (T(steps) * T(steps))
		        - value * value / T(steps)))
	{
	}

	/**
	 *
	 */
	std::size_t steps;

	/**
	 * Expectation value of the computation. The expectation value \f$ E \f$ is
	 * the average of the integrand \f$ f \f$ uniformly sampled at random points
	 * \f$ \vec{x} \in [0,1]^d \f$:
	 * \f[
	 *     E = \frac{1}{N} \sum_{i = 1}^N f \left( \vec{x}_i \right)
	 * \f]
	 */
	T value;

	/**
	 * Error (standard-deviation) of the expectation value. The error
	 * \f$ \sigma \f$ is computed as:
	 * \f[
	 *     \sigma = \sqrt{ \frac{1}{N} \left( \frac{1}{N} \sum_{i = 1}^N f^2
	 *     \left( \vec{x}_i \right) - E^2 \right) }
	 * \f]
	 */
	T error;
};

/**
 * PLAIN Monte Carlo integrator. \c plain integrates \c function over the
 * unit-hypercube with \c dimension dimensions using \c steps randomly chosen
 * points determined by the number generator \c random_number_generator. The
 * generator is seeded with \c seed. \c function must have the following
 * form:
 * \code
 * T integrand(std::vector<T> const& x, A const& variable)
 * {
 *     // return value of the function at x
 * }
 * \endcode
 * where \c T and \c A are substituted with the types the \c plain function is
 * called. \c variable is set by the \c plain function to the value given by
 * \c auxilliary_variable.
 *
 * The result of this function is an object of the type \c plain_result
 * containing both the value of the integral and an error.
 */
template <typename T, typename F, typename A, typename R = std::mt19937>
plain_result<T> plain(
	std::size_t dimensions,
	std::size_t steps,
	F& function,
	A const& auxilliary_variable,
	std::size_t seed = 0,
	R&& generator = std::mt19937()
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
