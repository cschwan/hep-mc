#ifndef HEP_MC_MULTI_CHANNEL_INTEGRAND_HPP
#define HEP_MC_MULTI_CHANNEL_INTEGRAND_HPP

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

/// \addtogroup integrands
/// @{

/// Class representing a function that can be integrated using the multi channel
/// algorithms.
template <typename T, typename F, typename D, bool distributions>
class multi_channel_integrand : public integrand<T, F, distributions>
{
public:
	/// Type of the density function.
	using density_type = D;

	/// Constructor. Instead of using the constructor directly you should
	/// consider using one of the helper functions \ref
	/// make_multi_channel_integrand.
	template <typename G, typename E>
	multi_channel_integrand(
		G&& function,
		std::size_t dimensions,
		E&& densities,
		std::size_t map_dimensions,
		std::size_t channels,
		std::vector<distribution_parameters<T>> const& parameters
	)
		: integrand<T, G, distributions>(std::forward<G>(function), dimensions,
			parameters)
		, densities_(std::forward<E>(densities))
		, map_dimensions_(map_dimensions)
		, channels_(channels)
	{
	}

	/// Returns the density functions.
	D& densities()
	{
		return densities_;
	}

	/// Returns the size of the vector the mappings map onto.
	std::size_t map_dimensions() const
	{
		return map_dimensions_;
	}

	/// Returns the number of channels.
	std::size_t channels() const
	{
		return channels_;
	}

private:
	D densities_;
	std::size_t map_dimensions_;
	std::size_t channels_;
};

/// Multi channel integrand constructor. For a description of the parameters see
/// \ref make_integrand. In addition, Multi Channel integrators need an
/// additonal function `densities` that compute the PDFs and the CDFs for a
/// randomly selected channel. This function would look like:
/// \code
/// T densities(
///     std::size_t channel,
///     std::vector<T> const& random_numbers,
///     std::vector<T>& coordinates,
///	    std::vector<T>& channel_densities
/// ) {
///     // for the selected `channel`, which is an integer from the half-open
///     // interval [0, channels) and the given `random_numbers` with 
///     // `random_numbers.size() == dimensions` compute the CDFs and store them
///     // in `coordinates` which is large as specified with `map_dimensions`.
///     // Also compute the PDFs and store them in `channel_densities`
///
///     return /* jacobian */;
/// }
/// \endcode
template <typename T, typename F, typename D>
inline multi_channel_integrand<T, F, D, false> make_multi_channel_integrand(
	F&& function,
	std::size_t dimensions,
	D&& densities,
	std::size_t map_dimensions,
	std::size_t channels
) {
	return multi_channel_integrand<T, F, D, false>(
		std::forward<F>(function),
		dimensions,
		std::forward<D>(densities),
		map_dimensions,
		channels,
		std::vector<distribution_parameters<T>>()
	);
}

/// Multi channel integrand constructor for distributions.
template <typename T, typename F, typename D, typename... Ds>
inline multi_channel_integrand<T, F, D, true> make_multi_channel_integrand(
	F&& function,
	std::size_t dimensions,
	D&& densities,
	std::size_t map_dimensions,
	std::size_t channels,
	Ds&&... parameters
) {
	return multi_channel_integrand<T, F, D, true>(
		std::forward<F>(function),
		dimensions,
		std::forward<D>(densities),
		map_dimensions,
		channels,
		std::vector<distribution_parameters<T>>{std::forward<Ds>(parameters)...}
	);
}

/// @}

}

#endif
