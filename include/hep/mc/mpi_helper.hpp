#ifndef HEP_MC_MPI_HELPER_HPP
#define HEP_MC_MPI_HELPER_HPP

/*
 * hep-mc - A Template Library for Monte Carlo Integration
 * Copyright (C) 2013-2014  Christopher Schwan
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

#include <mpi.h>

namespace
{

template <typename T>
MPI_Datatype mpi_datatype();

template <>
inline MPI_Datatype mpi_datatype<float>()
{
	return MPI_FLOAT;
}

template <>
inline MPI_Datatype mpi_datatype<double>()
{
	return MPI_DOUBLE;
}

template <>
inline MPI_Datatype mpi_datatype<long double>()
{
	return MPI_LONG_DOUBLE;
}

}

#endif
