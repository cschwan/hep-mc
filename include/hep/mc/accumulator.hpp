#ifndef HEP_MC_ACCUMULATOR_HPP
#define HEP_MC_ACCUMULATOR_HPP

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
#include "hep/mc/distribution_result.hpp"

#include <cstddef>
#include <type_traits>
#include <utility>
#include <vector>

namespace
{

template <typename T, bool distributions>
class accumulator;

template <typename T>
class accumulator<T, true>
{
public:
	accumulator(std::vector<hep::distribution_parameters<T>> const& parameters)
		: compensation_()
		, sum_()
		, sum_of_squares_()
		, bin_projector_(parameters)
	{
	}

	template <typename I, typename P>
	T invoke(I& integrand, P const& point)
	{
		// call the integrand function with the supplied point. Distributions
		// are generated here
		T const value = integrand.function()(point, bin_projector_) *
			point.weight();

		// perform kahan summation 'sum_ += value'
		T const y = value - compensation_;
		T const t = sum_ + y;
		compensation_ = (t - sum_) - y;
		sum_ = t;

		// no kahan summation for `sum_of_squares_`, should be OK without
		sum_of_squares_ += value * value;

		return value;
	}

	std::vector<hep::distribution_result<T>> distributions(
		std::size_t calls
	) const {
		return bin_projector_.distributions(calls);
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
	hep::bin_projector<T> bin_projector_;
};

template <typename T>
class accumulator<T, false>
{
public:
	accumulator(
		std::vector<hep::distribution_parameters<T>> const&
	)
		: compensation_()
		, sum_()
		, sum_of_squares_()
	{
	}

	template <typename I, typename P>
	T invoke(I& integrand, P const& point)
	{
		// call the integrand function with the supplied point. No distributions
		// are generated here
		T const value = integrand.function()(point) * point.weight();

		// perform kahan summation 'sum_ += value'
		T const y = value - compensation_;
		T const t = sum_ + y;
		compensation_ = (t - sum_) - y;
		sum_ = t;

		// no kahan summation for `sum_of_squares_`, should be OK without
		sum_of_squares_ += value * value;

		return value;
	}

	std::vector<hep::distribution_result<T>> distributions(std::size_t) const
	{
		// return empty distributions
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

template <typename I>
inline accumulator<typename I::numeric_type, I::has_distributions>
make_accumulator(
	I const& integrand
) {
	using T = typename I::numeric_type;
	constexpr bool has_distributions = I::has_distributions;

	return accumulator<T, has_distributions>(integrand.parameters());
}

}

#endif
