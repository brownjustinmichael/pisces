/*!***********************************************************************
 * \file one_d/element.cpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-08.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include <cmath>
#include "element_one_d.hpp"
#include "diffusion_one_d.hpp"
#include "../collocation/collocation.hpp"
#include "advection_one_d.h"
#include "solver_one_d.hpp"
#include "fftw_one_d.hpp"
#include "scale_one_d.hpp"
#include "boundary_one_d.hpp"
	
namespace one_d {
	advection_diffusion_element::advection_diffusion_element (int i_n, double *initial_conditions, int i_flags) : element (i_n, i_flags) {
		previous_timestep = 0.0;
		double diffusion_coeff = 1.0;
		double advection_coeff = 1.0;
		double alpha = 0.5;
		add_scalar (position);
		add_scalar (velocity);
		add_scalar (rhs);
		
		TRACE ("Initializing...");
		
		matrix.resize (i_n * i_n, 0.0);
		
		grid = std::make_shared<bases::collocation_grid> (chebyshev_grid (i_n, i_n, sqrt (2.0 / (i_n - 1.0))));

		angle_stream.reset (new io::incremental_output ("../output/test_angle", ".dat", 4, new io::header, i_n));
		angle_stream->append (&cell [0]);
		angle_stream->append ((*this) [position]);
		angle_stream->append ((*this) [velocity]);
		
		failsafe_dump.reset (new io::simple_output ("_dump.dat", i_n));
		failsafe_dump->append ((*this) [velocity]);
		
		add_boundary (boundary::make_unique (0.0, (*this) (velocity), 0.0, (*this) (velocity, n - 1)));
		add_boundary (boundary::make_unique (0.0, (*this) (rhs), 0.0, (*this) (rhs, n - 1)));
		add_implicit_plan (scale::make_unique (1.0, n * n, grid->get_data (0), &matrix [0]));
		add_implicit_plan (implicit_diffusion::make_unique (- diffusion_coeff * alpha, 0.0, 0.0, &timestep, i_n, grid, &matrix [0], &flags));
		add_explicit_grid_plan (scale::make_unique (0.0, n, (*this) (rhs)));
		add_explicit_grid_plan (explicit_diffusion::make_unique (diffusion_coeff * (1.0 - alpha), &timestep, i_n, grid, (*this) (velocity), (*this) (rhs), &flags));
		add_explicit_space_plan (advec::make_unique (n, &timestep, advection_coeff, (*this) (velocity), (*this) (rhs)));
		matrix_solver = lapack_solver::make_unique (n, (*this) (velocity), (*this) (rhs), &matrix [0], (*this) (velocity), &flags);
		transform_forward = fftw_cosine::make_unique (n, (*this) (velocity));
		
		double pioN = std::acos (-1.0) / i_n;
		for (int i = 0; i < i_n; ++i) {
			scalars [position] [i] = std::cos (pioN * i);
			scalars [velocity] [i] = initial_conditions [i];
		}
		
		transform_forward->execute ();
		TRACE ("Initialized.");
	}
} /* one_d */