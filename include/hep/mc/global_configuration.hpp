#ifndef HEP_MC_GLOBAL_CONFIGURATION
#define HEP_MC_GLOBAL_CONFIGURATION

/*
 * hep-mc - A Template Library for Monte Carlo Integration
 * Copyright (C) 2013  Christopher Schwan
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

/**
 * This function returns `true` the MPI functions use a single random number generator. This means
 * that the generator of each processor discards an appropriate amount of random numbers to
 * ascertain independent numbers are used. This make sure the result is independent of the number
 * of processesor and the same as the non-MPI routine. Since each process consumes as many random
 * numbers as calls are made in total, this may slow down the integration if the integrand evaluates
 * too fast. In this case a faster random number generator should be used.
 */
inline bool& mpi_single_generator()
{
	static bool use_single_generator = false;

	return use_single_generator;
}

/**
 * If `enabled` is `true` the MPI integration routines use a single random number generator so that
 * results are independent from the number of processors, in particular the same result than the
 * non-MPI integrators yield.
 *
 * \see \ref mpi_single_generator
 */
inline void mpi_single_generator(bool enabled)
{
	mpi_single_generator() = enabled;
}

/**
 * Returns `true` if the VEGAS routines should use the grid refinement method of CUBA \cite Cuba for
 * integrands with sharp edges. If the return value is `false` the original method is used. Default
 * is `false`.
 */
inline bool& vegas_cuba_refinement()
{
	static bool use_cuba_refinement = false;

	return use_cuba_refinement;
}

/**
 * If `enabled` is `true` the VEGAS routines use CUBA's grid refinement method.
 *
 * \see \ref vegas_cuba_refinement
 */
inline void vegas_cuba_refinement(bool enabled)
{
	vegas_cuba_refinement() = enabled;
}

/// @}

}

#endif
