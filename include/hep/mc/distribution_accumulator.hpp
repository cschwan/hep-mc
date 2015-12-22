#ifndef HEP_MC_DISTRIBUTION_ACCUMULATOR_HPP
#define HEP_MC_DISTRIBUTION_ACCUMULATOR_HPP

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

#include "hep/mc/distributions.hpp"
#include "hep/mc/distribution_result.hpp"
#include "hep/mc/mc_result.hpp"

#include <cstddef>
#include <type_traits>
#include <utility>
#include <vector>

namespace hep
{

/// \addtogroup internal
/// @{

template <typename T, typename D>
class distribution_accumulator;

template <typename T>
class distribution_accumulator<T, one_bin_projector>
{
public:
	distribution_accumulator(distributions<T, one_bin_projector> const&)
		: count_(0)
		, compensation_()
		, sum_()
		, sum_of_squares_()
	{
	}

	template <typename P>
	void add(P const&, T value)
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

	std::size_t count() const
	{
		return count_;
	}

	std::vector<std::vector<mc_result<T>>> distribution_results() const
	{
		return std::vector<std::vector<mc_result<T>>>();
	}

	T sum() const
	{
		return sum_;
	}

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

/// Creates an accumulator using the specified `distributions`.
template <typename D>
inline distribution_accumulator<
	typename std::remove_reference<D>::type::numeric_type,
	typename std::remove_reference<D>::type::projector
> make_distribution_accumulator(D&& distributions) {
	return std::forward<D>(distributions);
}

template <typename R, typename D, typename... A>
inline R make_result(D const& accumulator, A&&... args)
{
	return R(
		accumulator.distribution_results(),
		accumulator.count(),
		accumulator.sum(),
		accumulator.sum_of_squares(),
		std::forward<A>(args)...
	);
}

/// @}

}

#endif

