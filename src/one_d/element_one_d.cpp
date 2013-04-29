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
#include "../utilities/chebyshev.hpp"
#include "element_one_d.hpp"
#include "diffusion_one_d.hpp"
#include "advection_one_d.h"
#include "solver_one_d.hpp"
#include "fftw_one_d.hpp"
#include "scale_one_d.hpp"
#include "boundary_one_d.hpp"
	
namespace one_d
{
	namespace chebyshev
	{
		advection_diffusion_element::advection_diffusion_element (std::string name, int i_n, double initial_position, double *initial_velocity, int i_flags) : element (i_n, i_flags) {
			previous_timestep = 0.0;
			double diffusion_coeff = 10.0;
			double advection_coeff = -10.0;
			double alpha = 0.5;
		
			TRACE (logger, "Initializing...");
		
			matrix.resize (i_n * i_n, 0.0);
				
			grid = std::make_shared<bases::collocation_grid> (chebyshev_grid (i_n, i_n, sqrt (2.0 / (i_n - 1.0)), logger));

			transform_stream = std::make_shared <io::incremental_output> (io::incremental_output ("../output/test_cheb" + name, ".dat", 4, new io::header, i_n, logger));
			transform_stream->append (cell [0]);
			transform_stream->append ((*this) [velocity]);
		
			normal_stream = std::make_shared <io::incremental_output> (io::incremental_output ("../output/test_angle" + name, ".dat", 4, new io::header, i_n, logger));
			normal_stream->append (cell [0]);
			normal_stream->append ((*this) [position]);
			normal_stream->append ((*this) [velocity]);
			normal_stream->append ((*this) [rhs]);
		
			failsafe_dump = std::make_shared <io::simple_output> (io::simple_output ("dump" + name + ".dat", i_n, logger));
			failsafe_dump->append ((*this) [velocity]);
		
			timestep_plan = std::make_shared<bases::calculate_timestep> (bases::calculate_timestep (0.001, timestep));
		
			add_implicit_plan (std::make_shared <implicit_diffusion> (implicit_diffusion (- diffusion_coeff * alpha, timestep, i_n, grid, &matrix [0], &flags, logger)));
		
			add_explicit_grid_plan (std::make_shared <explicit_diffusion> (explicit_diffusion (diffusion_coeff * (1.0 - alpha), timestep, i_n, grid, (*this) (velocity), (*this) (rhs), &flags, logger)));
			add_explicit_space_plan (std::make_shared <advec> (advec (n, timestep, advection_coeff, (*this) (velocity), (*this) (rhs))));
		
			matrix_solver = std::make_shared <lapack_solver> (lapack_solver (n, (*this) (velocity), (*this) (rhs), &matrix [0], (*this) (velocity), &flags, logger));
		
			transform_forward = std::make_shared <fftw_cosine> (fftw_cosine (n, (*this) (velocity), &flags, logger));
		
			double pioN = std::acos (-1.0) / (n - 1);
			for (int i = 0; i < i_n; ++i) {
				(*this) (position, i) = std::cos (i * pioN) + initial_position;
				(*this) (velocity, i) = initial_velocity [i];
			}
		
			transform_forward->execute ();

			TRACE (logger, "Initialized.");
		}
	} /* chebyshev */
} /* one_d */