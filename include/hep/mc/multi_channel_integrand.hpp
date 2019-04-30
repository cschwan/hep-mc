#ifndef HEP_MC_MULTI_CHANNEL_INTEGRAND_HPP
#define HEP_MC_MULTI_CHANNEL_INTEGRAND_HPP

/*
 * hep-mc - A Template Library for Monte Carlo Integration
 * Copyright (C) 2016-2018  Christopher Schwan
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
#include "hep/mc/integrand.hpp"

#include <cstddef>
#include <type_traits>
#include <utility>
#include <vector>

namespace hep
{

/// \addtogroup integrands
/// @{

/// Class representing a function that can be integrated using the multi channel algorithms.
template <typename T, typename F, typename M, bool distributions>
class multi_channel_integrand : public integrand<T, F, distributions>
{
public:
    /// Type of the density function.
    using map_type = M;

    /// Constructor. Instead of using the constructor directly you should consider using one of the
    /// helper functions \ref make_multi_channel_integrand.
    template <typename G, typename N>
    multi_channel_integrand(
        G&& function,
        std::size_t dimensions,
        N&& map,
        std::size_t map_dimensions,
        std::size_t channels,
        std::vector<distribution_parameters<T>> const& parameters
    )
        : integrand<T, F, distributions>(std::forward<G>(function), dimensions, parameters)
        , map_(std::forward<N>(map))
        , map_dimensions_(map_dimensions)
        , channels_(channels)
    {
    }

    /// Returns the density functions.
    map_type& map()
    {
        return map_;
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
    M map_;
    std::size_t map_dimensions_;
    std::size_t channels_;
};

/// Template alias for a \ref multi_channel_integrand with its types `F` and `M` decayed with
/// `std::decay`.
template <typename T, typename F, typename M, bool distributions>
using multi_channel_integrand_type = multi_channel_integrand<T, typename std::decay<F>::type,
    typename std::decay<M>::type, distributions>;

/// Multi channel integrand constructor. For a description of the parameters see \ref
/// make_integrand. In addition, Multi Channel integrators need an additonal function `map` that
/// computes the PDFs and the CDFs for a randomly selected channel. This function should look like:
/// \code
/// T map(
///     std::size_t channel,
///     std::vector<T> const& random_numbers,
///     std::vector<T>& coordinates,
///        std::vector<T>& densities,
///     hep::multi_channel_map action
/// );
/// \endcode
/// This function is called by the multi channel integrator, first with the parameter `action` set
/// to \ref multi_channel_map::calculate_coordinates, which signals that the vector `coordinates`
/// must be filled using the CDFs. The return value is ignored for this function call. If the
/// integrand returns a non-zero value `map` is called again with `action` set to \ref
/// multi_channel_map::calculate_densities, which means that the vector `densities` must be
/// populated with all PDFs for the given `channel` and `random_numbers`. The return value is the
/// jacobian for the given `channel`.
template <typename T, typename F, typename M>
inline multi_channel_integrand_type<T, F, M, false>
make_multi_channel_integrand(
    F&& function,
    std::size_t dimensions,
    M&& map,
    std::size_t map_dimensions,
    std::size_t channels
) {
    return multi_channel_integrand_type<T, F, M, false>(
        std::forward<F>(function),
        dimensions,
        std::forward<M>(map),
        map_dimensions,
        channels,
        std::vector<distribution_parameters<T>>()
    );
}

/// Multi channel integrand constructor for distributions.
template <typename T, typename F, typename M, typename... Ds>
inline multi_channel_integrand_type<T, F, M, true> make_multi_channel_integrand(
    F&& function,
    std::size_t dimensions,
    M&& map,
    std::size_t map_dimensions,
    std::size_t channels,
    Ds&&... parameters
) {
    return multi_channel_integrand_type<T, F, M, true>(
        std::forward<F>(function),
        dimensions,
        std::forward<M>(map),
        map_dimensions,
        channels,
        std::vector<distribution_parameters<T>>{std::forward<Ds>(parameters)...}
    );
}

/// @}

}

#endif
