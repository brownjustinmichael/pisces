/*!***********************************************************************
 * \file element_one_d.cpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-08.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include <cmath>
#include <cassert>
#include <string>
#include <memory>
#include "../config.hpp"
#include "../bases/timestep.hpp"
#include "../utils/chebyshev.hpp"
#include "element_one_d.hpp"
#include "diffusion_one_d.hpp"
#include "advection_one_d.h"
#include "solver_one_d.hpp"
#include "fftw_one_d.hpp"
#include "boundary_one_d.hpp"
	
namespace one_d
{
	namespace chebyshev
	{
		advection_diffusion_element::advection_diffusion_element (int i_n, double i_position_0, double i_position_n, std::string i_name, io::parameter_map& inputParams, int i_flags) : element (i_n, i_position_0, i_position_n, i_name, inputParams, i_flags) {
			double delta_t = inputParams["time_step_size"].asDouble; 
			double diffusion_coeff = inputParams["diffusion_coeff"].asDouble;
			double advection_coeff = inputParams["advection_coeff"].asDouble; 
			double alpha = 0.5;
			std::vector <double> init (n);
		
			assert (n > 0);
		
			TRACE (logger, "Initializing..." << logger);
		
			matrix.resize (i_n * i_n, 0.0);
			
			// Set up output
			normal_stream = std::make_shared <io::incremental_output> (io::incremental_output ("../output/test_angle_" + name + "_", ".dat", 4, new io::header, i_n, inputParams["output_every"].asInt));
			normal_stream->append (cell [0]);
			normal_stream->append ((*this) [position]);
			normal_stream->append ((*this) [velocity]);
			normal_stream->append ((*this) [rhs]);
			
			// Set up plans in order
			add_plan (std::make_shared <explicit_diffusion> (explicit_diffusion (diffusion_coeff * (1.0 - alpha), i_n, grid, velocity, position, rhs)));

			set_transform (std::make_shared <fftw_cosine> (fftw_cosine (n, velocity)));
			if (advection_coeff != 0.0) {
				add_plan (std::make_shared <advec> (advec (n, advection_coeff, velocity, rhs, grid)));
			}
			// add_plan (std::make_shared<constant_timestep> (constant_timestep (delta_t, timestep)));
			add_plan (std::make_shared <implicit_diffusion> (implicit_diffusion (- diffusion_coeff * alpha, i_n, grid, &matrix [0])));
		
			// Set up solver
			set_solver (std::make_shared <solver> (solver (n, timestep, grid->get_data (0), &matrix [0], velocity, rhs)));
			
			normal_stream->to_file ();
		
			TRACE (logger, "Initialized.");
		}
		
		double advection_diffusion_element::calculate_timestep () {
			return inputParams["time_step_size"].asDouble;
		}
	} /* chebyshev */
} /* one_d */