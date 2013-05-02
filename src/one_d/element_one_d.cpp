/*!***********************************************************************
 * \file one_d/element.cpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-08.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include <cmath>
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
		advection_diffusion_element::advection_diffusion_element (std::string i_name, int i_n, double initial_position, double *initial_velocity, int i_flags) : element (i_name, i_n, i_flags) {
			double diffusion_coeff = 1.0;
			double advection_coeff = 0.0;
			double alpha = 0.5;
		
			TRACE (logger, "Initializing..." << logger);
		
			matrix.resize (i_n * i_n, 0.0);
			
			normal_stream = std::make_shared <io::incremental_output> (io::incremental_output ("../output/test_angle_" + name + "_", ".dat", 4, new io::header, i_n));
			normal_stream->append (cell [0]);
			normal_stream->append ((*this) [position]);
			normal_stream->append ((*this) [velocity]);
			normal_stream->append ((*this) [rhs]);
									
			set_grid (std::make_shared<bases::collocation_grid> (chebyshev_grid (i_n, i_n, sqrt (2.0 / (i_n - 1.0)), logger)));
			
			set_timestep (std::make_shared<bases::calculate_timestep> (bases::calculate_timestep (0.00001, timestep)));
		
			add_implicit_plan (std::make_shared <implicit_diffusion> (implicit_diffusion (- diffusion_coeff * alpha, timestep, i_n, grid, &matrix [0])));
			
			add_explicit_grid_plan (std::make_shared <explicit_diffusion> (explicit_diffusion (diffusion_coeff * (1.0 - alpha), timestep, i_n, grid, velocity, position, rhs)));
			add_explicit_space_plan (std::make_shared <advec> (advec (n, timestep, advection_coeff, velocity, rhs, grid)));
		
			set_solver (std::make_shared <solver> (solver (n, &matrix [0], velocity, rhs)));
		
			set_transform (std::make_shared <fftw_cosine> (fftw_cosine (n, velocity)));
		
			double pioN = std::acos (-1.0) / (n - 1);
			for (int i = 0; i < i_n; ++i) {
				(*this) (position, i) = std::cos (i * pioN) + initial_position;
				(*this) (velocity, i) = initial_velocity [i];
			}

			TRACE (logger, "Initialized.");
		}
	} /* chebyshev */
} /* one_d */