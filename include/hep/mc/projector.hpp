#ifndef HEP_MC_PROJECTOR_HPP
#define HEP_MC_PROJECTOR_HPP

/*
 * hep-mc - A Template Library for Monte Carlo Integration
 * Copyright (C) 2016  Christopher Schwan
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

#include "hep/mc/distribution_parameters.hpp"
#include "hep/mc/distribution_result.hpp"
#include "hep/mc/mc_result.hpp"

#include <cstddef>
#include <vector>

namespace hep
{

/// \addtogroup distributions
/// @{

/// Class that project the multi-dimensional Monte Carlo point into a discrete
/// bin and generates the corresponding distributions.
template <typename T>
class projector
{
public:
	/// Constructor.
	projector(std::vector<distribution_parameters<T>> const& parameters)
		: parameters_(parameters)
		, compensation_(parameters.size())
		, sum_(parameters.size())
		, sum_of_squares_(parameters.size())
	{
		for (std::size_t i = 0; i != parameters.size(); ++i)
		{
			std::size_t const bins = parameters[i].bins();

			compensation_[i].resize(bins);
			sum_[i].resize(bins);
			sum_of_squares_[i].resize(bins);
		}
	}

	/// Projects a point for the distribution with the specified `index` to the
	/// x-axis at the value `projection` and sets the value to `value`. Note
	/// that you have to manually multiply with the corresponding weight.
	void add(std::size_t index, T projection, T value)
	{
		T const x = projection - parameters_[index].x_min();

		if (x < T())
		{
			return;
		}

		std::size_t const bin = x / parameters_[index].bin_size();

		if (bin >= parameters_[index].bins())
		{
			return;
		}

		// kahan summation for each bin
		T const y = value - compensation_[index][bin];
		T const t = sum_[index][bin] + y;
		compensation_[index][bin] = (t - sum_[index][bin]) - y;
		sum_[index][bin] = t;

		sum_of_squares_[index][bin] += value * value;
	}

	/// Returns the accumulated distributions for `calls` ellapsed calls.
	std::vector<distribution_result<T>> distributions(std::size_t calls) const
	{
		std::vector<distribution_result<T>> result;
		result.reserve(sum_.size());

		// loop over all distributions
		for (std::size_t dist = 0; dist != sum_.size(); ++dist)
		{
			std::vector<T> mid_points;
			std::vector<mc_result<T>> bin_results;

			mid_points.reserve(sum_[dist].size());
			bin_results.reserve(sum_[dist].size());

			auto const& params = parameters_[dist];
			T const inv_bin_size = T(1.0) / params.bin_size();

			// loop over the bins of the current distribution
			for (std::size_t bin = 0; bin != sum_[dist].size(); ++bin)
			{
				mid_points.push_back(params.x_min() +
					T(bin + 0.5) * params.bin_size());

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

private:
	std::vector<distribution_parameters<T>> parameters_;
	std::vector<std::vector<T>> compensation_;
	std::vector<std::vector<T>> sum_;
	std::vector<std::vector<T>> sum_of_squares_;
};

/// @}

}

#endif
