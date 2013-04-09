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
	diffusion_element::diffusion_element (int i_n, int i_flags) {
		int i;
		n = i_n;
		flags = i_flags;
		cell.resize (i_n);
		position.resize (i_n);
		velocity.resize (i_n, 0.0);
		
		TRACE ("Initializing...");
		
		double pioN = std::acos (-1.0) / i_n;
		for (i = 0; i < i_n; ++i) {
			cell [i] = i;
			position [i] = std::cos (pioN * i);
		}
		
		velocity [0] = 2.0;
		velocity [2] = -1.0;

		angle_stream.reset (new io::incremental_output ("../output/test_angle", ".dat", 4, new io::header, i_n));
		angle_stream->append (&cell [0]);
		angle_stream->append (&position [0]);
		angle_stream->append (&velocity [0]);
		
		failsafe_dump.reset (new io::simple_output ("_dump.dat", i_n));
		failsafe_dump->append (&velocity [0]);
		
		diffusion_plan.reset (new diffusion::cheb_1D (2., 0.5, i_n, &velocity [0], &velocity [0], i_flags));
		fourier_plan = fftw_plan_r2r_1d (i_n, &velocity [0], &velocity [0], FFTW_REDFT00, FFTW_ESTIMATE);
		
		TRACE ("Initialized.");
	}
	
	void diffusion_element::update () {
		int i;
		
		TRACE ("Updating...");
		
		try {
			// Should be replaced by a CFL check
			double timestep = 0.1;
		
			// Calculate the diffusion in Chebyshev space
			diffusion_plan->execute (timestep);
			
			// Transform forward
			fftw_execute (fourier_plan);
		
			// Output in angle space
			angle_stream->to_file ();
					
			// Transform backward
			fftw_execute (fourier_plan);
		
			// We rescale the values to account for the factor of 2(N-1) that occurs during FFT
			for (i = 0; i < n; ++i) {
				velocity [i] /= (2 * (n - 1));
			}
			
		} catch (io::exceptions::file_exception &io_exception) {
			ERROR ("Exception caught, failsafe dump to _dump.dat");
			failsafe_dump->to_file ();
			exit (EXIT_FAILURE);
		}
		
		TRACE ("Update complete");
	
	}
} /* element */