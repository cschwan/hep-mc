#ifndef HEP_MC_MPI_VEGAS_CALLBACK_HPP
#define HEP_MC_MPI_VEGAS_CALLBACK_HPP

/*
 * hep-mc - A Template Library for Monte Carlo Integration
 * Copyright (C) 2014-2019  Christopher Schwan
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

#include "hep/mc/vegas_callback.hpp"
#include "hep/mc/vegas_chkpt.hpp"

#include <mpi.h>

namespace hep
{

/// \addtogroup callbacks
/// @{

/// The default callback function. This function does nothing and always returns `true`. It is the
/// MPI equivalent of \ref vegas_silent_callback.
template <typename T>
inline bool mpi_vegas_silent_callback(MPI_Comm, vegas_chkpt<T> const&)
{
    return true;
}

/// Callback function that prints a detailed summary about every iteration performed so far. This
/// function always returns `true`. It is the equivalent of \ref vegas_verbose_callback and only
/// writes an output if it was called from rank zero to avoid duplicated output.
///
/// \see vegas_callback
template <typename T>
inline bool mpi_vegas_verbose_callback(MPI_Comm communicator, vegas_chkpt<T> const& chkpt)
{
    int rank = -1;
    MPI_Comm_rank(communicator, &rank);

    if (rank == 0)
    {
        vegas_verbose_callback<T>(chkpt);
    }

    return true;
}

/// @}

}

#endif
