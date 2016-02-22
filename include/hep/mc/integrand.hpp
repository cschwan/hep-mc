#ifndef HEP_MC_INTEGRAND_HPP
#define HEP_MC_INTEGRAND_HPP

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

#include <cstddef>
#include <utility>
#include <vector>

namespace hep
{

/// \addtogroup integrand_group
/// @{

/// Class representing a function that can be integrated using the PLAIN-like
/// algorithms, which currently are PLAIN and VEGAS.
template <typename T, typename F, bool distributions>
class integrand
{
public:
	/// Numeric type of the integrand.
	using numeric_type = T;

	/// The type of the integrand function that is integrated.
	using function_type = F;

	/// Signals whether this integrand wants to generate distributions or not.
	static constexpr bool has_distributions = distributions;

	/// Constructor. Instead of using the constructor directly you should
	/// consider using one of the helper functions \ref make_integrand.
	template <typename G>
	integrand(
		G&& function,
		std::size_t dimensions,
		std::vector<distribution_parameters<T>> const& parameters
	)
		: function_(std::forward<G>(function))
		, parameters_(parameters)
		, dimensions_(dimensions)
	{
	}

	/// Returns the dimension of the integrand.
	std::size_t dimensions() const
	{
		return dimensions_;
	}

	/// Returns the integrand's function.
	F& function()
	{
		return function_;
	}

	/// Returns the parameters of the distribution(s).
	std::vector<distribution_parameters<T>> const& parameters() const
	{
		return parameters_;
	}

private:
	F function_;
	std::vector<distribution_parameters<T>> parameters_;
	std::size_t dimensions_;
};

/// PLAIN/VEGAS constructor. This function constructs an integrand using the
/// given `function` that must accept points from the \f$ d \d$-dimensional
/// hypercube, where \f$ d \f$ is given by the parameter `dimensions`.
template <typename T, typename F>
inline integrand<T, F, false> make_integrand(
	F&& function,
	std::size_t dimensions
) {
	return integrand<T, F, false>(
		std::forward<F>(function),
		dimensions,
		std::vector<distribution_parameters<T>>()
	);
}

/// PLAIN/VEGAS distributions constructor.
template <typename T, typename F, typename... D>
inline integrand<T, F, true> make_integrand(
	F&& function,
	std::size_t dimensions,
	D&&... parameters
) {
	return integrand<T, F, true>(
		std::forward<F>(function),
		dimensions,
		std::vector<distribution_parameters<T>>{std::forward<D>(parameters)...}
	);
}

/// @}

}

#endif
