/*!**********************************************************************
 * \file messenger.cpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2013-08-17.
* Copyright 2013 Justin Brown. All rights reserved.
************************************************************************/

#include "messenger.hpp"
#ifdef _MPI
#include "mpi.h"
#endif // _MPI

namespace bases
{
#ifdef _MPI
	template <>
	MPI_Datatype mpi_type <double> () {
		return MPI::DOUBLE;
	}
	
	template <>
	MPI_Datatype mpi_type <float> () {
		return MPI::REAL;
	}
	
	template <>
	MPI_Datatype mpi_type <int> () {
		return MPI::INT;
	}
#endif // _MPI
} /* bases */
