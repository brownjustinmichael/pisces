/*!***********************************************************************
 * \file one_d/element.cpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-08.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include <cmath>
#include <string>
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
	void element::calculate () {		
		TRACE (logger, "Calculating...");
				
		try {
			// Start in grid space
			
			TRACE (logger, "Updating timestep...");
			
			// Testing
			// Should be replaced by a CFL check
			timestep = 0.001;
			
			TRACE (logger, "Executing explicit grid plans...");
			
			for (int i = 0; i < n_explicit_grid_plans; ++i) {
				explicit_grid_plans [i]->execute ();
			}
			
			TRACE (logger, "Transforming to normal space...");
			
			// Switch to normal space
			if (transform_forward) {
				transform_forward->execute ();
			} else {
				WARN (logger, "Transform not defined. It is likely the element was not set up correctly.")
			}

			TRACE (logger, "Writing to file...");
			
			// Output in angle space
			if (angle_stream) {
				angle_stream->to_file ();
			}
			
			TRACE (logger, "Executing explicit space plans...");

			for (int i = 0; i < n_explicit_space_plans; ++i) {
				explicit_space_plans [i]->execute ();
			}
			
			if (timestep != previous_timestep) {
				TRACE (logger, "Executing implicit plans...");
				
				flags &= ~factorized;
				for (int i = 0; i < n_implicit_plans; ++i) {
					implicit_plans [i]->execute ();
				}
			}

			previous_timestep = timestep;
		} catch (...) {
			ERROR (logger, "Exception caught, failsafe dump to _dump.dat");
			failsafe_dump->to_file ();
			exit (EXIT_FAILURE);
		}
		
		TRACE (logger, "Calculation complete.");
	}
	
	advection_diffusion_element::advection_diffusion_element (std::string name, double i_alpha_plus, double i_alpha_minus, int i_n, double *initial_position, double *initial_velocity, int i_flags) : element (i_n, i_flags) {
		previous_timestep = 0.0;
		double diffusion_coeff = 0.1;
		double advection_coeff = -0.0;
		double alpha = 0.5;
		alpha_plus = i_alpha_plus;
		alpha_minus = i_alpha_minus;
		add_scalar (position);
		add_scalar (velocity);
		add_scalar (rhs);
		
		TRACE (logger, "Initializing...");
		
		matrix.resize (i_n * i_n, 0.0);
		
		grid = std::make_shared<bases::collocation_grid> (chebyshev_grid (i_n, i_n, sqrt (2.0 / (i_n - 1.0)), logger));

		angle_stream.reset (new io::incremental_output ("../output/test_angle" + name, ".dat", 4, new io::header, i_n, logger));
		angle_stream->append (cell [0]);
		angle_stream->append ((*this) [position]);
		angle_stream->append ((*this) [velocity]);
		
		failsafe_dump.reset (new io::simple_output ("dump" + name + ".dat", i_n, logger));
		failsafe_dump->append ((*this) [velocity]);
		
		add_implicit_plan (scale::make_unique (1.0, n * n, grid->get_data (0), &matrix [0], &flags, logger));
		add_implicit_plan (implicit_diffusion::make_unique (- diffusion_coeff * alpha, alpha_plus, alpha_minus, &timestep, i_n, grid, &matrix [0], &flags, logger));
		add_explicit_grid_plan (scale::make_unique (0.0, n, (*this) (rhs), &flags, logger));
		add_explicit_grid_plan (explicit_diffusion::make_unique (diffusion_coeff * (1.0 - alpha), &timestep, i_n, grid, (*this) (velocity), (*this) (rhs), &flags, logger));
		// add_explicit_space_plan (advec::make_unique (n, &timestep, advection_coeff, (*this) (velocity), (*this) (rhs)));
		matrix_solver = lapack_solver::make_unique (n, (*this) (velocity), (*this) (rhs), &matrix [0], (*this) (velocity), &flags, logger);
		transform_forward = fftw_cosine::make_unique (n, (*this) (velocity), &flags, logger);
		
		for (int i = 0; i < i_n + 1; ++i) {
			scalars [position] [i] = initial_position [i];
			scalars [velocity] [i] = initial_velocity [i];
		}
		
		transform_forward->execute ();

		TRACE (logger, "Initialized.");
	}
} /* one_d */