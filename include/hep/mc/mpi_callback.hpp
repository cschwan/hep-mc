#ifndef HEP_MC_MPI_CALLBACK_HPP
#define HEP_MC_MPI_CALLBACK_HPP

/*
 * hep-mc - A Template Library for Monte Carlo Integration
 * Copyright (C) 2019  Christopher Schwan
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

#include "hep/mc/callback.hpp"

#include <string>

#include <mpi.h>

namespace hep
{

/// \addtogroup callbacks
/// @{

/// Default callback type used in all MPI integration algorithms.
template <typename Checkpoint>
class mpi_callback
{
public:
    /// Numeric type used by the checkpoint.
    using numeric_type = typename Checkpoint::result_type::numeric_type;

    /// Constructor.
    mpi_callback(
        callback_mode mode = callback_mode::verbose,
        std::string const& filename = "",
        numeric_type target_rel_err = numeric_type()
    )
        : callback_{mode, filename, target_rel_err}
    {
    }

    /// Callback function. Calls the corresponding non-MPI function if the rank of the process is
    /// `0`.
    bool operator()(MPI_Comm communicator, Checkpoint const& chkpt)
    {
        if (callback_.mode() != callback_mode::silent)
        {
            int rank;
            MPI_Comm_rank(communicator, &rank);

            if (rank != 0)
            {
                callback_.mode(callback_mode::silent);
            }
        }

        return callback_(chkpt);
    }

private:
    callback<Checkpoint> callback_;
};

/// @}

}

#endif
