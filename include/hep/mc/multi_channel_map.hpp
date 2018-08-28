#ifndef HEP_MC_MULTI_CHANNEL_MAP_HPP
#define HEP_MC_MULTI_CHANNEL_MAP_HPP

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

namespace hep
{

/// \addtogroup integrands
/// @{

/// Enumeration that specifies which parameters of the user specified map should be calculated.
enum class multi_channel_map
{
    /// Signals the map-function that the parameter `coordinates` should be calculated. The
    /// parameter `densities` can be calculated at time point.
    calculate_coordinates,

    /// Signals that the map function is called with the previously calculated `coordinates` and
    /// that the paramter `densities` now must be calculated. If the integrand does not contribute,
    /// i.e. if the user-defined function returned zero, then this step will be skipped because the
    /// densities are not needed.
    calculate_densities
};

/// @}

}

#endif
