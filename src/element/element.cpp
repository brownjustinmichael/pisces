/*!***********************************************************************
 * \file element.cpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-08.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include <cmath>
#include <fftw3.h>
#include "../config.hpp"
#include "element.hpp"
#include "../io/exceptions.hpp"
#include "../diffusion/diffusion.hpp"
#include "../advection/advection.h"
#include "../boundary/boundary.hpp"
#include "../solver/solver.hpp"
#include "../fft/fft.hpp"
#include "../explicit_plans/lapack.hpp"

namespace element
{
	void element::calculate () {		
		TRACE ("Calculating...");
				
		try {
			// Start in grid space
			
			TRACE ("Updating timestep...");
			
			// Testing
			// Should be replaced by a CFL check
			timestep = 0.001;
			
			TRACE ("Executing explicit grid plans...");
			
			for (int i = 0; i < n_explicit_grid_plans; ++i) {
				explicit_grid_plans [i]->execute ();
			}
			
			TRACE ("Transforming to normal space...");
			
			// Switch to normal space
			if (transform_forward) {
				transform_forward->execute ();
			} else {
				WARN ("Transform not defined. It is likely the element was not set up correctly.")
			}
			
			TRACE ("Writing to file...");
			
			// Output in angle space
			if (angle_stream) {
				angle_stream->to_file ();
			}
			
			TRACE ("Executing explicit space plans...");

			for (int i = 0; i < n_explicit_space_plans; ++i) {
				explicit_space_plans [i]->execute ();
			}
			
			if (timestep != previous_timestep) {
				TRACE ("Executing implicit plans...");
				
				flags &= ~solver::factorized;
				for (int i = 0; i < n_implicit_plans; ++i) {
					implicit_plans [i]->execute ();
				}
			}

			previous_timestep = timestep;
		} catch (io::exceptions::file_exception &io_exception) {
			ERROR ("Exception caught, failsafe dump to _dump.dat");
			failsafe_dump->to_file ();
			exit (EXIT_FAILURE);
		}
		
		TRACE ("Calculation complete.");
	}
		
	void element::execute_boundaries () {
		TRACE ("Executing boundaries...");
		
		for (int i = 0; i < n_boundaries; ++i) {
			boundaries [i]->execute ();
		}
		
		TRACE ("Boundaries executed.");
	}
	
	void element::update () {
		TRACE ("Updating...");
		
		if (matrix_solver) {
			matrix_solver->solve ();
		} else {
			WARN ("No matrix solver defined. It is likely the element was not set up correctly.")
		}
		
		TRACE ("Update complete");
	}
	
	diffusion_element_1D::diffusion_element_1D (int i_n, int i_flags) : element_1D (i_n, i_flags) {
		previous_timestep = 0.0;
		double diffusion_coeff = 1.0;
		double advection_coeff = 1.0;
		double alpha = 0.5;
		double scale = 1.0;
		double sigma = 0.1;
		add_scalar (position);
		add_scalar (velocity);
		add_scalar (rhs);
		
		TRACE ("Initializing...");
		
		matrix.resize (i_n * i_n, 0.0);
		
		grid = std::make_shared<collocation::collocation_grid> (collocation::chebyshev_grid (i_n, i_n, sqrt (2.0 / (i_n - 1.0))));

		angle_stream.reset (new io::incremental_output ("../output/test_angle", ".dat", 4, new io::header, i_n));
		angle_stream->append (&cell [0]);
		angle_stream->append (&(scalars [position]) [0]);
		angle_stream->append (&(scalars [velocity]) [0]);
		
		failsafe_dump.reset (new io::simple_output ("_dump.dat", i_n));
		failsafe_dump->append (&(scalars [velocity]) [0]);
				
		add_boundary (one_d::boundary::make_unique (0.0, &(scalars [velocity] [0]), 0.0, &(scalars [velocity] [n - 1])));
		add_boundary (one_d::boundary::make_unique (0.0, &(scalars [rhs] [0]), 0.0, &(scalars [rhs] [n - 1])));
		add_implicit_plan (explicit_plans::copy::make_unique (n * n, grid->get_data (0), &matrix [0]));
		add_implicit_plan (one_d::implicit_plans::diffusion::make_unique (- diffusion_coeff * alpha, 0.0, 0.0, &timestep, i_n, grid, &matrix [0], &flags));
		add_explicit_grid_plan (explicit_plans::zero::make_unique (n, &(scalars [rhs]) [0]));
		add_explicit_grid_plan (one_d::explicit_plans::diffusion::make_unique (diffusion_coeff * (1.0 - alpha), &timestep, i_n, grid, &(scalars [velocity] [0]), &(scalars [rhs]) [0], &flags));
		add_explicit_space_plan (advection::advec_1D::make_unique (n, &timestep, advection_coeff, &(scalars [velocity]) [0], &(scalars [rhs]) [0]));
		matrix_solver = solver::lapack_solver::make_unique (n, &(scalars [velocity]) [0], &(scalars [rhs]) [0], &matrix [0], &(scalars [velocity]) [0], &flags);
		transform_forward = fft::fftw_cosine::make_unique (n, &(scalars [velocity]) [0], &(scalars [velocity]) [0]);
		
		double pioN = std::acos (-1.0) / i_n;
		for (int i = 0; i < i_n; ++i) {
			scalars [position] [i] = std::cos (pioN * i);
			scalars [velocity] [i] = scale * std::exp (- scalars [position] [i] * scalars [position] [i] / 2.0 / sigma / sigma) - scale * std::exp (- 1.0 / 2.0 / sigma / sigma);
		}
		
		transform_forward->execute ();
		TRACE ("Initialized.");
	}
} /* element */