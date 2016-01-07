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

#include "hep/mc/distribution_projector.hpp"
#include "hep/mc/distribution_result.hpp"

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
		, projector_(projector)
		, compensation_(projector.parameters().size())
		, sum_(projector.parameters().size())
		, sum_of_squares_(projector.parameters().size())
		, projections_(projector.parameters().size())
	{
		for (std::size_t i = 0; i != projector.parameters().size(); ++i)
		{
			std::size_t const bins = projector.parameters()[i].bins();

			compensation_[i].resize(bins);
			sum_[i].resize(bins);
			sum_of_squares_[i].resize(bins);
		}
	}

	template <typename M>
	void add(M const& point, T value)
	{
		// `accumulator` takes care of the total integration result
		accumulator_.add(point, value);

		// project the point
		projector_.projector()(point, projections_);

		for (std::size_t i = 0; i != sum_.size(); ++i)
		{
			T const x = projections_[i] - projector_.parameters()[i].x_min();

			if (x < T())
			{
				// point lies left from our leftmost bin
				continue;
			}

			std::size_t const j = x / projector_.parameters()[i].bin_size();

			if (j >= projector_.parameters()[i].bins())
			{
				// point lies right from our rightmost bin
				continue;
			}

			// kahan summation for each bin
			T const y = value - compensation_[i][j];
			T const t = sum_[i][j] + y;
			compensation_[i][j] = (t - sum_[i][j]) - y;
			sum_[i][j] = t;

			sum_of_squares_[i][j] += value * value;
		}
	}

	std::vector<hep::distribution_result<T>> distributions(
		std::size_t calls
	) const {
		std::vector<hep::distribution_result<T>> result;
		result.reserve(sum_.size());

		// loop over all distributions
		for (std::size_t dist = 0; dist != sum_.size(); ++dist)
		{
			std::vector<T> mid_points;
			std::vector<hep::mc_result<T>> bin_results;

			mid_points.reserve(sum_[dist].size());
			bin_results.reserve(sum_[dist].size());

			auto const& parameters = projector_.parameters()[dist];
			T const inv_bin_size = T(1.0) / parameters.bin_size();

			// loop over the bins of the current distribution
			for (std::size_t bin = 0; bin != sum_[dist].size(); ++bin)
			{
				mid_points.push_back(parameters.x_min() +
					T(bin + 0.5) * parameters.bin_size());

				bin_results.emplace_back(
					calls,
					inv_bin_size * sum_[dist][bin],
					inv_bin_size * inv_bin_size * sum_of_squares_[dist][bin]
				);
			}

			result.emplace_back(mid_points, bin_results);
		}

		return result;
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
	hep::distribution_projector<T, P> projector_;
	std::vector<std::vector<T>> compensation_;
	std::vector<std::vector<T>> sum_;
	std::vector<std::vector<T>> sum_of_squares_;
	std::vector<T> projections_;
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
