#ifndef HEP_MC_MC_RESULT_HPP
#define HEP_MC_MC_RESULT_HPP

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

#include <cmath>
#include <cstddef>

namespace hep
{

/**
 * The estimation of a Monte Carlo integration. Every Monte Carlo integrator
 * returns one or more instances of this class.
 */
template <typename T>
struct mc_result
{
	/**
	 * Constructor.
	 */
	mc_result(std::size_t calls, T sum, T sum_of_squares)
		: calls(calls)
		, value(sum / T(calls))
		, error(std::sqrt(sum_of_squares - value * value * T(calls)) / T(calls))
	{
	}

	/**
	 * The number of function evaluations \f$ N \f$ performed to obtain this
	 * result.
	 */
	std::size_t calls;

	/**
	 * Expectation value \f$ E \f$ of this result.
	 */
	T value;

	/**
	 * Error \f$ S \f$ of the expectation value.
	 */
	T error;
};

/**
 * Creates a \ref mc_result using the parameters `calls` `value` and `error`.
 */
template <typename T>
inline mc_result<T> create_result(std::size_t calls, T value, T error)
{
	T sum = T(calls) * value;
	T sum_of_squares = T(calls) * (value * value + T(calls) * error * error);

	return mc_result<T>(calls, sum, sum_of_squares);
}

}

#endif
