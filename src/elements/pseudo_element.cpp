/*!**********************************************************************
<<<<<<< HEAD
 * \file pseudo.cpp
=======
 * \file boussinesq.cpp
>>>>>>> devel
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2014-02-03.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#include "pseudo_element.hpp"

#include "plans-solvers/solvers.hpp"

namespace pisces
{
	using namespace plans;

	template <class datatype>
	pseudo_element <datatype>::pseudo_element (grids::axis i_axis_n, grids::axis i_axis_m, int i_name, io::parameters& i_params, data::data <datatype> &i_data, mpi::messenger* i_messenger_ptr, int i_element_flags, bool load_diffusion) : 
	boussinesq_element <datatype> (i_axis_n, i_axis_m, i_name, i_params, i_data, i_messenger_ptr, i_element_flags, load_diffusion) {
		pressure.resize (m);

		for (int j = 0; j < m; ++j)
		{
			pressure [j] = 1.0 + (*grids [1]) [0] * 0.1;
		}

		// Set up the velocity constraint
		*pdiv <datatype> (equations ["pressure"], equations ["x_velocity"], equations ["z_velocity"], &pressure [0])
		==
		0.0;
		
	TRACE ("Initialized.");
	}

	template class pseudo_element <double>;
} /* pisces */
