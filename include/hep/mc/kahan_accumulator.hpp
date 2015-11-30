#ifndef HEP_MC_KAHAN_ACCUMALATOR_HPP
#define HEP_MC_KAHAN_ACCUMALATOR_HPP

/*
 * hep-mc - A Template Library for Monte Carlo Integration
 * Copyright (C) 2015  Christopher Schwan
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

namespace hep
{

/// Implements kahan summation.
template <typename T>
class kahan_accumulator
{
public:
	kahan_accumulator()
		: count_(0)
		, compensation_()
		, sum_()
		, sum_of_squares_()
	{
	}

	/// Adds `value` to he accumulator and updates the internal parameters.
	void add(T value)
	{
		// perform kahan summation 'sum_ += value'
		T const y = value - compensation_;
		T const t = sum_ + y;
		compensation_ = (t - sum_) - y;
		sum_ = t;

		// no kahan summation for `sum_of_squares_`, should be OK without
		sum_of_squares_ += value * value;

		++count_;
	}

	/// Returns the number of times \ref add was called.
	std::size_t count() const
	{
		return count_;
	}

	/// Returns the sum of all previously added values
	T sum() const
	{
		return sum_;
	}

	/// Returns the sum of all previously added squares of values.
	T sum_of_squares() const
	{
		return sum_of_squares_;
	}

private:
	std::size_t count_;
	T compensation_;
	T sum_;
	T sum_of_squares_;
};

}

#endif
