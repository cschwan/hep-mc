#ifndef HEP_MC_MC_RESULT_HPP
#define HEP_MC_MC_RESULT_HPP

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

namespace hep
{

/**
 * The result of a Monte Carlo integration.
 */
template <typename T>
struct mc_result
{
	/**
	 * Constructor. Sets \ref steps, \ref value, and \ref error.
	 */
	mc_result(std::size_t samples, T sum, T sum_of_squares)
		: samples(samples)
		, value(sum / T(samples))
		, error(std::sqrt(sum_of_squares - value * value * T(samples)) /
			T(samples))
	{
	}

	/**
	 * The number of Monte Carlo steps performed for obtaining this result.
	 */
	std::size_t samples;

	/**
	 * Expectation value of this result. The expectation value \f$ E \f$ is the
	 * average of the integrand \f$ f \f$ uniformly sampled over \f$ N \f$
	 * random points \f$ \vec{x}_i \in [0,1]^d \f$:
	 * \f[
	 *     E = \frac{1}{N} \sum_{i = 1}^N f \left( \vec{x}_i \right)
	 * \f]
	 */
	T value;

	/**
	 * Error of the expectation value. The error \f$ S \f$ is computed as:
	 * \f[
	 *     S = \sqrt{ \frac{1}{N} \left( \frac{1}{N} \sum_{i = 1}^N f^2
	 *     \left( \vec{x}_i \right) - E^2 \right) }
	 * \f]
	 */
	T error;
};

}

#endif
