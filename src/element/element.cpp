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
#include "../fft/fft.hpp"

namespace element
{
	diffusion_element_1D::diffusion_element_1D (int i_n, int i_flags) : element_1D (i_n) {
		int i;
		flags = i_flags;
		previous_timestep = 0.0;
		double diffusion_coeff = 1.0;
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
		
		implicit_diffusion.reset (new diffusion::implicit_methods::collocation_chebyshev_1D (- diffusion_coeff * alpha, &timestep, i_n, grid, &matrix [0], &flags));
		explicit_diffusion.reset (new diffusion::explicit_methods::collocation_chebyshev_1D (diffusion_coeff * (1.0 - alpha), &timestep, i_n, grid, &(scalars [velocity]) [0], &(scalars [rhs]) [0], &flags));
		matrix_solver.reset (new solver::lapack_solver (n, &(scalars [velocity]) [0], &(scalars [rhs]) [0], &matrix [0], &(scalars [velocity]) [0], &flags));
		fourier_transform.reset (new fft::fftw_cosine (n, &(scalars [velocity]) [0], &(scalars [velocity]) [0]));
		
		TRACE ("Initialized.");
	}
	
	void diffusion_element_1D::update () {
		int nn = n * n;
		int i, ione = 1;
		
		TRACE ("Updating...");
				
		try {
			// Testing
			// Should be replaced by a CFL check
			timestep = 0.01;
			
			for (i = 0; i < n; ++i) {
				scalars [rhs] [i] = 0.0;
			}

			explicit_diffusion->execute ();
			
			fourier_transform->execute ();
			
			// Output in angle space
			angle_stream->to_file ();

			// Set up and evaluate the implicit part of the diffusion equation
			if (timestep != previous_timestep) {
				dcopy_ (&nn, grid->get_data (0), &ione, &matrix [0], &ione);
				flags &= ~solver::factorized;

				implicit_diffusion->execute ();
			}
		
			matrix_solver->solve ();

			previous_timestep = timestep;
			
		} catch (io::exceptions::file_exception &io_exception) {
			ERROR ("Exception caught, failsafe dump to _dump.dat");
			failsafe_dump->to_file ();
			exit (EXIT_FAILURE);
		}
		
		TRACE ("Update complete");
	

	}
} /* element */