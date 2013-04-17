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

namespace element
{
	void element_1D::calculate () {
		int i;
		
		TRACE ("Updating...");
				
		try {
			// Start in grid space
			// Testing
			// Should be replaced by a CFL check
			timestep = 0.01;
			
			for (i = 0; i < n_explicit_grid_plans; ++i) {
				explicit_grid_plans [i]->execute ();
			}
			
			// Switch to normal space
			if (transform_forward) {
				transform_forward->execute ();
			}
			
			// Output in angle space
			if (angle_stream) {
				angle_stream->to_file ();
			}

			for (i = 0; i < n_explicit_space_plans; ++i) {
				explicit_space_plans [i]->execute ();
			}
			
			if (timestep != previous_timestep) {
				flags &= ~solver::factorized;
				for (i = 0; i < n_implicit_plans; ++i) {
					implicit_plans [i]->execute ();
				}
			}

			previous_timestep = timestep;
		} catch (io::exceptions::file_exception &io_exception) {
			ERROR ("Exception caught, failsafe dump to _dump.dat");
			failsafe_dump->to_file ();
			exit (EXIT_FAILURE);
		}
	}
		
	void element_1D::execute_boundaries () {
		for (int i = 0; i < n_boundaries; ++i) {
			boundaries [i]->execute ();
		}
	}
	
	void element_1D::update () {
		if (matrix_solver) {
			matrix_solver->solve ();
		}
		
		TRACE ("Update complete");
	}
	
	diffusion_element_1D::diffusion_element_1D (int i_n, int i_flags) : element_1D (i_n, i_flags) {
		int i;
		previous_timestep = 0.0;
		double diffusion_coeff = 1.0;
		double advection_coeff = 1.0;
		double alpha = 0.5;
		add_scalar (position);
		add_scalar (velocity);
		add_scalar (rhs);
		
		TRACE ("Initializing...");
		
		matrix.resize (i_n * i_n, 0.0);
		
		double pioN = std::acos (-1.0) / i_n;
		for (i = 0; i < i_n; ++i) {
			scalars [position] [i] = std::cos (pioN * i);
		}
		
		scalars [velocity] [0] = 2.0;
		scalars [velocity] [2] = -1.0;
		
		grid.reset (new collocation::chebyshev_grid (i_n, i_n, sqrt (2.0 / (i_n - 1.0))));

		angle_stream.reset (new io::incremental_output ("../output/test_angle", ".dat", 4, new io::header, i_n));
		angle_stream->append (&cell [0]);
		angle_stream->append (&(scalars [position]) [0]);
		angle_stream->append (&(scalars [velocity]) [0]);
		
		failsafe_dump.reset (new io::simple_output ("_dump.dat", i_n));
		failsafe_dump->append (&(scalars [velocity]) [0]);
		
		DEBUG ("scalars [velocity] = " << &(scalars [velocity] [0]) << " [0] = " << scalars [velocity] [0])
		DEBUG ("scalars [velocity] = " << &(scalars [velocity] [n - 1]) << " [n - 1] = " << scalars [velocity] [n - 1])
		
		add_boundary (std::unique_ptr<plan> (new boundary::boundary_1D (0.0, &(scalars [velocity] [0]), 0.0, &(scalars [velocity] [n - 1]))));
		add_boundary (std::unique_ptr<plan> (new boundary::boundary_1D (0.0, &(scalars [rhs] [0]), 0.0, &(scalars [rhs] [n - 1]))));
		add_implicit_plan (std::unique_ptr<plan> (new copy (n * n, grid->get_data (0), &matrix [0])));
		add_implicit_plan (std::unique_ptr<plan> (new diffusion::implicit_methods::collocation_chebyshev_1D (- diffusion_coeff * alpha, &timestep, i_n, grid, &matrix [0], &flags)));
		add_explicit_grid_plan (std::unique_ptr<plan> (new zero (n, &(scalars [rhs]) [0])));
		add_explicit_grid_plan (std::unique_ptr<plan> (new diffusion::explicit_methods::collocation_chebyshev_1D (diffusion_coeff * (1.0 - alpha), &timestep, i_n, grid, &(scalars [velocity]) [0], &(scalars [rhs]) [0], &flags)));
		add_explicit_space_plan (std::unique_ptr<plan> (new advection::advec_1D (n, &timestep, advection_coeff, &(scalars [velocity]) [0], &(scalars [rhs]) [0])));
		matrix_solver.reset (new solver::lapack_solver (n, &(scalars [velocity]) [0], &(scalars [rhs]) [0], &matrix [0], &(scalars [velocity]) [0], &flags));
		transform_forward.reset (new fft::fftw_cosine (n, &(scalars [velocity]) [0], &(scalars [velocity]) [0]));
		TRACE ("Initialized.");
	}
} /* element */