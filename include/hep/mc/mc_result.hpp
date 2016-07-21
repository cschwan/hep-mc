#ifndef HEP_MC_MC_RESULT_HPP
#define HEP_MC_MC_RESULT_HPP

/*
 * hep-mc - A Template Library for Monte Carlo Integration
 * Copyright (C) 2012-2016  Christopher Schwan
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
#include <vector>

namespace hep
{

/// \addtogroup results
/// @{

/**
 * The estimation of a Monte Carlo integration. Every Monte Carlo integrator
 * returns one or more instances of this class. The PLAIN Monte Carlo
 * integrator, for example, calculates the parameters as follows:
 * \f{align}{
 *     E &= \frac{1}{N} \sum_{i=1}^N f ( \vec{x}_i ) \\
 *     S^2 &= \frac{1}{N-1} \left[ \frac{1}{N} \sum_{i=1}^N f^2 ( \vec{x}_i )
 *            - E^2  \right]
 * \f}
 */
template <typename T>
class mc_result
{
public:
	/// The numeric type used for member variables.
	using numeric_type = T;

	/// Constructor.
	mc_result(std::size_t calls, T sum, T sum_of_squares)
		: calls_(calls)
		, sum_(sum)
		, sum_of_squares_(sum_of_squares)
	{
	}

	/// The number of function evaluations \f$ N \f$ performed to obtain this
	/// result.
	std::size_t calls() const
	{
		return calls_;
	}

	/// Expectation value \f$ E \f$ of this result.
	T value() const
	{
		return sum_ / T(calls_);
	}

	/// Variance \f$ S^2 \f$ of the expectation value.
	T variance() const
	{
		return (sum_of_squares_ - sum_ * sum_ / T(calls_)) / T(calls_)
			/ T(calls_ - 1);
	}

	/// Standard deviation \f$ S \f$ of the expectation value.
	T error() const
	{
		return std::sqrt(variance());
	}

	/// Returns the sum, i.e. \f$ \sum_{i=1}^N f ( \vec{x}_i ) \f$.
	T sum() const
	{
		return sum_;
	}

	/// Returns the sum of squares, i.e. \f$ \sum_{i=1}^N f^2 ( \vec{x}_i ) \f$.
	T sum_of_squares() const
	{
		return sum_of_squares_;
	}

private:
	std::size_t calls_;
	T sum_;
	T sum_of_squares_;
};

/// Creates a \ref mc_result using the parameters `calls`, `value` and `error`.
template <typename T>
inline mc_result<T> create_result(std::size_t calls, T value, T error)
{
	T sum = T(calls) * value;
	T sum_of_squares = T(calls) * (value * value +
		T(calls - 1) * error * error);

	return mc_result<T>(calls, sum, sum_of_squares);
}

/// @}

}

#endif
