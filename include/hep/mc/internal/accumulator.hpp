#ifndef HEP_MC_INTERNAL_ACCUMULATOR_HPP
#define HEP_MC_INTERNAL_ACCUMULATOR_HPP

/*
 * hep-mc - A Template Library for Monte Carlo Integration
 * Copyright (C) 2015-2017  Christopher Schwan
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

#include "hep/mc/internal/accumulator_fwd.hpp"
#include "hep/mc/projector.hpp"
#include "hep/mc/distribution_parameters.hpp"
#include "hep/mc/distribution_result.hpp"
#include "hep/mc/plain_result.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <vector>

namespace
{

template <typename T>
inline void accumulate(T& sum, T& sum_of_squares, T& compensation, T value)
{
	T const y = value - compensation;
	T const t = sum + y;
	compensation = (t - sum) - y;
	sum = t;

	sum_of_squares += value * value;
}

template <typename T>
class accumulator<T, true>
{
public:
	accumulator(std::vector<hep::distribution_parameters<T>> const& parameters)
		: parameters_(parameters)
		, indices_()
		, sums_()
		, compensations_()
	{
		std::size_t index = 2;
		indices_.reserve(parameters.size());

		for (auto const& params : parameters)
		{
			indices_.push_back(index);
			index += 2 * params.bins();
		}

		sums_.resize(2 + index);
		compensations_.resize(index / 2);
	}

	template <typename I, typename P>
	T invoke(I& integrand, P const& point)
	{
		using std::isfinite;

		hep::projector<T> projector(*this, point);

		// call the integrand function with the supplied point. Distributions
		// are generated here
		T value = integrand.function()(point, projector);

		if (value != T())
		{
			value *= point.weight();

			if (isfinite(value))
			{
				accumulate(sums_[0], sums_[1], compensations_[0], value);
			}
			else
			{
				value = T();
			}
		}

		return value;
	}

	void add_to_distribution(std::size_t index, T projection, T value)
	{
		using std::isfinite;

		if (!isfinite(value))
		{
			return;
		}

		// TODO: index might be larger than the than allowed; throw?

		T const x = projection - parameters_[index].x_min();

		if (x < T())
		{
			// point is left of the range that we are binning
			return;
		}

		std::size_t const bin = x / parameters_[index].bin_size();

		if (bin >= parameters_[index].bins())
		{
			// point is right of the range that we are binning
			return;
		}

		std::size_t const new_index = indices_[index] + 2 * bin;

		accumulate(
			sums_[new_index],
			sums_[new_index + 1],
			compensations_[new_index / 2],
			value
		);
	}

	std::vector<hep::distribution_parameters<T>> const& parameters() const
	{
		return parameters_;
	}

	hep::plain_result<T> result(std::size_t calls) const
	{
		std::vector<hep::distribution_result<T>> result;
		result.reserve(parameters_.size());

		std::size_t index = 2;

		// loop over all distributions
		for (auto const& params : parameters_)
		{
			std::vector<hep::mc_result<T>> bin_results;
			bin_results.reserve(params.bins());

			T const inv_bin_size = T(1.0) / params.bin_size();

			// loop over the bins of the current distribution
			for (std::size_t bin = 0; bin != params.bins(); ++bin)
			{
				bin_results.emplace_back(
					calls,
					inv_bin_size                * sums_[index],
					inv_bin_size * inv_bin_size * sums_[index + 1]
				);

				index += 2;
			}

			result.emplace_back(params, bin_results);
		}

		return hep::plain_result<T>(result, calls, sums_[0], sums_[1]);
	}

private:
	std::vector<hep::distribution_parameters<T>> parameters_;
	std::vector<std::size_t> indices_;
	std::vector<T> sums_;
	std::vector<T> compensations_;
};

template <typename T>
class accumulator<T, false>
{
public:
	accumulator(std::vector<hep::distribution_parameters<T>> const&)
		: sums_()
	{
	}

	template <typename I, typename P>
	T invoke(I& integrand, P const& point)
	{
		// call the integrand function with the supplied point. No distributions
		// are generated here
		T value = integrand.function()(point);

		if (value != T())
		{
			value *= point.weight();

			if (std::isfinite(value))
			{
				accumulate(sums_[0], sums_[1], sums_[2], value);
			}
			else
			{
				value = T();
			}
		}

		return value;
	}

	hep::plain_result<T> result(std::size_t calls) const
	{
		return hep::plain_result<T>(
			std::vector<hep::distribution_result<T>>{},
			calls,
			sums_[0],
			sums_[1]
		);
	}

private:
	std::array<T, 3> sums_;
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

namespace hep
{

template <typename T>
inline void projector<T>::add(std::size_t index, T projection, T value)
{
	// grant selective access to the following function (only)
	accumulator_.add_to_distribution(
		index,
		projection,
		value * point_.weight()
	);
}

template <typename T>
inline std::vector<distribution_parameters<T>> const&
projector<T>::parameters() const {
	return accumulator_.parameters();
}

}

#endif
