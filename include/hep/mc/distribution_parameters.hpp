#ifndef HEP_MC_DISTRIBUTION_PARAMETERS_HPP
#define HEP_MC_DISTRIBUTION_PARAMETERS_HPP

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
#include <string>

namespace hep
{

/// \addtogroup distributions
/// @{

/// Defines the parameters of a one-dimensional differential distribution.
template <typename T>
class distribution_parameters
{
public:
	/// Constructor.
	distribution_parameters(
		std::size_t bins,
		T x_min,
		T x_max,
		std::string const& name
	)
		: bins_(bins)
		, x_min_(x_min)
		, bin_size_((x_max - x_min) / bins)
		, name_(name)
	{
	}

	/// Number of bins for this distribution.
	std::size_t bins() const
	{
		return bins_;
	}

	/// Name of the distribution.
	std::string const& name() const
	{
		return name_;
	}

	/// Highest value of the x-axis that is still part of the differential
	/// distribution.
	T x_max() const
	{
		return x_min_ + bin_size_ * bins_;
	}

	/// Lowest value of the x-axis that is still part of the differential
	/// distribution.
	T x_min() const
	{
		return x_min_;
	}

	/// Size of each bin in every distribution.
	T bin_size() const
	{
		return bin_size_;
	}

private:
	std::size_t bins_;
	T x_min_;
	T bin_size_;
	std::string name_;
};

/// Shortcut for calling the constructor that automatically determines the
/// numeric type.
template <typename T>
distribution_parameters<T> make_dist_params(
	std::size_t bins,
	T x_min,
	T x_max,
	std::string const& name = ""
) {
	return distribution_parameters<T>(bins, x_min, x_max, name);
}

/// @}

}

#endif

