/*!**********************************************************************
 * \file diffusion.cpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2014-02-03.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#include "diffusion.hpp"

#include <sstream>
#include <iomanip>
#include <stdlib.h>
#include <memory>

#include "io/formats/netcdf.hpp"
#include "io/formats/ascii.hpp"

#include "plans/plan.hpp"
#include "plans/advection.hpp"
#include "plans/diffusion.hpp"
#include "plans/source.hpp"
#include "plans-solvers/boundaries/implemented_boundary.hpp"
#include "plans-solvers/solvers.hpp"
#include "io/functors/div.hpp"

namespace pisces
{
	using namespace plans;
	using namespace plans::solvers;
	using namespace boundaries;

	std::string diffusion_element::class_name() {
		return "diffusion";
	}

	implemented_element::registrar <diffusion_element> diffusion_registrar ("diffusion");

	std::shared_ptr <element> diffusion_element::instance (grids::axis i_axis_n, grids::axis i_axis_m, int i_name, io::parameters& i_params, data::data &i_data, mpi::messenger* i_messenger_ptr, int i_element_flags) {
		return std::shared_ptr <element> (new diffusion_element (i_axis_n, i_axis_m, i_name, i_params, i_data, i_messenger_ptr, i_element_flags));
	}
	
	diffusion_element::diffusion_element (grids::axis i_axis_n, grids::axis i_axis_m, int i_name, io::parameters& i_params, data::data &i_data, mpi::messenger* i_messenger_ptr, int i_element_flags) : 
	implemented_element (i_axis_n, i_axis_m, i_name, i_params, i_data, i_messenger_ptr, i_element_flags) {
		TRACE ("Initializing...");
		// Set up the scalar equation
		*split_solver (equations ["scalar"], timestep, 
			dirichlet (i_params ["equations.scalar.bottom.value"].as <double> ()), 
			dirichlet (i_params ["equations.scalar.top.value"].as <double> ())) 
		== 
		params ["equations.scalar.diffusion"] * diff ();	
	TRACE ("Initialized.");
	}
	
	double diffusion_element::calculate_timestep (int i, int j, formats::virtual_file *virtual_file) {
		return 1.0 / 0.0;
	}
} /* pisces */
