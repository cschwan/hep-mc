#ifndef HEP_MC_PLAIN_HPP
#define HEP_MC_PLAIN_HPP

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

#include <cstddef>
#include <random>

namespace hep
{

/**
 * The result of a PLAIN Monte Carlo integration.
 */
template <typename T>
struct plain_result
{
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
	std::size_t dimension,
	std::size_t steps,
	std::size_t seed,
	F& function,
	A const& auxilliary_variable,
	R&& random_number_generator = std::mt19937()
);

}

#include <hep/mc/impl/plain.hpp>

#endif
