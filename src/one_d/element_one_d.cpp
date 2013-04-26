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
#include "element_one_d.hpp"
#include "diffusion_one_d.hpp"
#include "../collocation/chebyshev.hpp"
#include "advection_one_d.h"
#include "solver_one_d.hpp"
#include "fftw_one_d.hpp"
#include "scale_one_d.hpp"
#include "boundary_one_d.hpp"
	
namespace one_d {
	advection_diffusion_element::advection_diffusion_element (std::string name, double i_alpha_0, double i_alpha_n, int i_n, double initial_position, double *initial_velocity, int i_flags) : element (i_n, i_flags) {
		previous_timestep = 0.0;
		double diffusion_coeff = 1.0;
		double advection_coeff = -0.0;
		double alpha = 1.0;
		alpha_0 = i_alpha_0;
		alpha_n = i_alpha_n;
		add_scalar (position);
		add_scalar (velocity);
		add_scalar (rhs);
		
		std::shared_ptr <bases::plan> temp_shared;
		
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
		
		temp_shared = std::make_shared <scale> (scale (1.0, n * n, grid->get_data (0), &matrix [0], &flags, logger));
		add_implicit_plan (temp_shared);
		
		temp_shared->execute ();
		add_implicit_plan (std::make_shared <implicit_diffusion> (implicit_diffusion (- diffusion_coeff * alpha, 0.0, 0.0, &timestep, i_n, grid, &matrix [0], &flags, logger)));
		add_explicit_grid_plan (std::make_shared <scale> (scale (0.0, n, (*this) (rhs), &flags, logger)));
		add_explicit_grid_plan (std::make_shared <explicit_diffusion> (explicit_diffusion (diffusion_coeff * (1.0 - alpha), alpha_0, alpha_n, &timestep, i_n, grid, (*this) (velocity), (*this) (rhs), &flags, logger)));
		
		/*
			TODO Because implicit diffusion isn't working at boundaries, need to compensate explicit scheme.
		*/
		
		// add_explicit_space_plan (advec::make_shared (n, &timestep, advection_coeff, (*this) (velocity), (*this) (rhs)));
		matrix_solver = std::make_shared <lapack_solver> (lapack_solver (n, (*this) (velocity), (*this) (rhs), &matrix [0], (*this) (velocity), &flags, logger));
		transform_forward = std::make_shared <fftw_cosine> (fftw_cosine (n, (*this) (velocity), &flags, logger));
		
		double pioN = std::acos (-1.0) / (n - 1);
		for (int i = 0; i < i_n + 1; ++i) {
			scalars [position] [i] = std::cos (i * pioN) + initial_position;
			scalars [velocity] [i] = initial_velocity [i];
		}
		
		transform_forward->execute ();
		// matrix_solver->solve ();

		TRACE (logger, "Initialized.");
	}
} /* one_d */