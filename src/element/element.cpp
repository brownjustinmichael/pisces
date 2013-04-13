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

namespace element
{
	diffusion_element::diffusion_element (int i_n, int i_flags) : element_1D (i_n) {
		int i;
		flags = i_flags;
		add_scalar (position);
		add_scalar (velocity);
		add_scalar (rhs);
		
		TRACE ("Initializing...");
		
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
		
		implicit_diffusion.reset (new diffusion::collocation_chebyshev_implicit_1D (1., 0.5, i_n, grid, &(scalars [velocity]) [0], &(scalars [rhs]) [0], &(scalars [velocity]) [0], i_flags));
		explicit_diffusion.reset (new diffusion::collocation_chebyshev_explicit_1D (0.5, i_n, grid, &(scalars [velocity]) [0], &(scalars [rhs]) [0], i_flags));
		fourier_plan = fftw_plan_r2r_1d (i_n, &(scalars [velocity]) [0], &(scalars [velocity]) [0], FFTW_REDFT00, FFTW_ESTIMATE);
		
		TRACE ("Initialized.");
	}
	
	void diffusion_element::update () {
		int i;
		
		TRACE ("Updating...");
				
		try {
			// Testing
			// Should be replaced by a CFL check
			double timestep = 0.000001;
		
			// Calculate the diffusion in Chebyshev space
			explicit_diffusion->execute (timestep);
			
			// Transform forward
			fftw_execute (fourier_plan);

			// We rescale the values to account for the factor of 2(N-1) that occurs during FFT
			for (i = 0; i < n; ++i) {
				(scalars [velocity]) [i] /= sqrt (2.0 * (n - 1));
			}
			
			// Output in angle space
			angle_stream->to_file ();

			implicit_diffusion->execute (timestep);
			
			for (i = 0; i < n; ++i) {
				scalars [rhs] [i] = 0.0;
			}
			
		} catch (io::exceptions::file_exception &io_exception) {
			ERROR ("Exception caught, failsafe dump to _dump.dat");
			failsafe_dump->to_file ();
			exit (EXIT_FAILURE);
		}
		
		TRACE ("Update complete");
	
	}
} /* element */