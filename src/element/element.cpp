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
#include "../solver/solver.hpp"

namespace element
{
	diffusion_element_1D::diffusion_element_1D (int i_n, int i_flags) : element_1D (i_n) {
		int i;
		flags = i_flags;
		previous_timestep = 0.0;
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
		
		grid.reset (new collocation::chebyshev_grid (i_n, i_n, sqrt (2.0 / (i_n - 1))));

		angle_stream.reset (new io::incremental_output ("../output/test_angle", ".dat", 4, new io::header, i_n));
		angle_stream->append (&cell [0]);
		angle_stream->append (&(scalars [position]) [0]);
		angle_stream->append (&(scalars [velocity]) [0]);
		
		failsafe_dump.reset (new io::simple_output ("_dump.dat", i_n));
		failsafe_dump->append (&(scalars [velocity]) [0]);
		
		implicit_diffusion.reset (new diffusion::implicit_methods::collocation_chebyshev_1D (-5, i_n, grid, &matrix [0], i_flags));
		explicit_diffusion.reset (new diffusion::explicit_methods::collocation_chebyshev_1D (5, i_n, grid, &(scalars [velocity]) [0], &(scalars [rhs]) [0], i_flags));
		matrix_solver.reset (new solver::lapack_solver (n, &(scalars [velocity]) [0], &(scalars [rhs]) [0], &matrix [0], &(scalars [velocity]) [0]));
		fourier_plan = fftw_plan_r2r_1d (i_n, &(scalars [velocity]) [0], &(scalars [velocity]) [0], FFTW_REDFT00, FFTW_ESTIMATE);
		
		TRACE ("Initialized.");
	}
	
	void diffusion_element_1D::update () {
		int i;
		int nn = n * n;
		int ione = 1;
		std::vector<double> temp (n);
		std::vector<double> matrix_copy (n * n);
		
		TRACE ("Updating...");
				
		try {
			// Testing
			// Should be replaced by a CFL check
			double timestep = 0.01;
			
			for (i = 0; i < n; ++i) {
				(scalars [rhs]) [i] = 0.0;
			}	

			explicit_diffusion->execute (timestep, &flags);

			// Transform backward
			fftw_execute (fourier_plan);
			
			// We rescale the values to account for the factor of 2(N-1) that occurs during FFT
			for (i = 0; i < n; ++i) {
				(scalars [velocity]) [i] /= sqrt (2.0 * (n - 1));
			}
			
			// Output in angle space
			angle_stream->to_file ();

			if (timestep != previous_timestep) {
				flags &= ~solver::factorized;
				dcopy_ (&nn, grid->get_data (0), &ione, &matrix [0], &ione);
				
				// Calculate the diffusion in Chebyshev space
				implicit_diffusion->execute (timestep, &flags);
			}

			matrix_solver->solve (&flags);
			
			previous_timestep = timestep;
			
		} catch (io::exceptions::file_exception &io_exception) {
			ERROR ("Exception caught, failsafe dump to _dump.dat");
			failsafe_dump->to_file ();
			exit (EXIT_FAILURE);
		}
		
		TRACE ("Update complete");
	

	}
} /* element */