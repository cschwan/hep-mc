#ifndef HEP_MC_MPI_DATAYPE_HPP
#define HEP_MC_MPI_DATAYPE_HPP

#include <mpi.h>

namespace hep
{

/**
 * Returns the correct \c MPI_Datatype for the template type \c T.
 */
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
