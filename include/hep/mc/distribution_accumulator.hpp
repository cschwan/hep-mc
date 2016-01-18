#ifndef HEP_MC_DISTRIBUTION_ACCUMULATOR_HPP
#define HEP_MC_DISTRIBUTION_ACCUMULATOR_HPP

/*
 * hep-mc - A Template Library for Monte Carlo Integration
 * Copyright (C) 2015-2016  Christopher Schwan
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

#include "hep/mc/bin_projector.hpp"
#include "hep/mc/distribution_projector.hpp"
#include "hep/mc/distribution_result.hpp"
#include "hep/mc/function_value.hpp"

#include <cstddef>
#include <type_traits>
#include <utility>
#include <vector>

namespace
{

template <typename T, typename P>
class distribution_accumulator
{
public:
	distribution_accumulator(hep::distribution_projector<T, P> const& projector)
		: accumulator_(hep::default_projector<T>())
		, bin_projector_(projector.parameters())
		, bin_projector_function_(projector.projector())
	{
	}

	template <typename M, typename F>
	void add(M const& point, F const& function, T value)
	{
		// `accumulator` takes care of the total integration result
		accumulator_.add(point, function, value);

		bin_projector_function_(
			point,
			bin_projector_,
			hep::function_value2<T, F>(function, value)
		);
	}

	std::vector<hep::distribution_result<T>> distributions(
		std::size_t calls
	) const {
		return bin_projector_.distributions(calls);
	}

	T sum() const
	{
		return accumulator_.sum();
	}

	T sum_of_squares() const
	{
		return accumulator_.sum_of_squares();
	}

private:
	distribution_accumulator<T, hep::one_bin_projector> accumulator_;
	hep::bin_projector<T> bin_projector_;
	P bin_projector_function_;
};

template <typename T>
class distribution_accumulator<T, hep::one_bin_projector>
{
public:
	distribution_accumulator(hep::default_projector<T> const&)
		: compensation_()
		, sum_()
		, sum_of_squares_()
	{
	}

	template <typename M, typename F>
	void add(M const&, F const&, T value)
	{
		// perform kahan summation 'sum_ += value'
		T const y = value - compensation_;
		T const t = sum_ + y;
		compensation_ = (t - sum_) - y;
		sum_ = t;

		// no kahan summation for `sum_of_squares_`, should be OK without
		sum_of_squares_ += value * value;
	}

	std::vector<hep::distribution_result<T>> distributions(std::size_t) const
	{
		return std::vector<hep::distribution_result<T>>();
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
	T compensation_;
	T sum_;
	T sum_of_squares_;
};

template <typename P>
inline distribution_accumulator<
	typename std::remove_reference<P>::type::numeric_type,
	typename std::remove_reference<P>::type::projector_type
> make_distribution_accumulator(P&& projector) {
	return std::forward<P>(projector);
}

}

#endif
