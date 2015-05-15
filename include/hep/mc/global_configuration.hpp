#ifndef HEP_MC_GLOBAL_CONFIGURATION_HPP
#define HEP_MC_GLOBAL_CONFIGURATION_HPP

/*
 * hep-mc - A Template Library for Monte Carlo Integration
 * Copyright (C) 2014-2015  Christopher Schwan
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

/// \addtogroup global_configuration
/// @{

/// Returns `true` if the VEGAS routines should use the grid refinement method
/// of CUBA \cite Cuba for integrands with sharp edges. If the return value is
/// `false` the original method is used. Default is `false`.
inline bool& vegas_cuba_refinement()
{
	static bool use_cuba_refinement = false;

	return use_cuba_refinement;
}

/// If `enabled` is `true` the VEGAS routines use CUBA's grid refinement method.
///
/// \see \ref vegas_cuba_refinement
inline void vegas_cuba_refinement(bool enabled)
{
	vegas_cuba_refinement() = enabled;
}

/// @}

}

#endif
