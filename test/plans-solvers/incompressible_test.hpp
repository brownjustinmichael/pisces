/*!**********************************************************************
 * \file incompressible_test.cpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2013-11-05.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include <cxxtest/TestSuite.h>

#include "mpi/messenger.hpp"
#include <vector>
#include <stdio.h>
#include "linalg/linalg.hpp"
#include "plans/grids/grid.hpp"
#include "plans-solvers/solvers.hpp"

#define TEST_TINY 1.e-6

class plan_solver_test_suite : public CxxTest::TestSuite
{	
public:
	void test_incompressible () {
		int n = 200, m = 300;
		int flags;

		grids::horizontal::grid grid_n (n, -1.0, 1.0);
		grids::vertical::grid grid_m (m, -1.0, 1.0);

		logger::log_config::set_severity (0);
		
		grids::variable x_velocity (grid_n, grid_m, flags), z_velocity (grid_n, grid_m, flags), data (grid_n, grid_m, flags), rhs (grid_n, grid_m, flags), rhs_compare (grid_n, grid_m, flags);

		double *z_ptr = z_velocity.ptr (real_spectral), *x_ptr = x_velocity.ptr (real_spectral);
		mpi::messenger mess;

		auto solver = plans::solvers::incompressible (&mess, NULL, NULL, data, data, x_velocity, z_velocity);

		solver.factorize ();

		for (int q = 0; q < 1; ++q)
		{
			for (int j = 0; j < m; ++j) {
				for (int i = 0; i < n; ++i) {
					x_ptr [i * m + j] = (rand () % 1000) / 1000.;
					z_ptr [i * m + j] = (rand () % 1000) / 1000.;
					data [i * m + j] = (rand () % 1000) / 1000.;
					rhs [i * m + j] = 0.0;
				}
			}
			
			double scalar = acos (-1.0) * 2.0 / 2.0;

			solver.factorize ();

			solver.execute ();

			double total = 0.;
			for (int i = 2; i < n; i += 2) {
				for (int j = 1; j < m - 1; ++j) {
					total += ((z_ptr [i * m + j + 1] - z_ptr [i * m + j - 1]) / (grid_m [j + 1] - grid_m [j - 1]) - scalar * (i / 2) * x_ptr [(i + 1) * m + j]) * ((z_ptr [i * m + j + 1] - z_ptr [i * m + j - 1]) / (grid_m [j + 1] - grid_m [j - 1]) - scalar * (i / 2) * x_ptr [(i + 1) * m + j]);
					total += ((z_ptr [(i + 1) * m + j + 1] - z_ptr [(i + 1) * m + j - 1]) / (grid_m [j + 1] - grid_m [j - 1]) + scalar * (i / 2) * x_ptr [i * m + j]) * ((z_ptr [(i + 1) * m + j + 1] - z_ptr [(i + 1) * m + j - 1]) / (grid_m [j + 1] - grid_m [j - 1]) + scalar * (i / 2) * x_ptr [i * m + j]);
				}
			}
			TSM_ASSERT_DELTA ("Incompressibility failure", 0.0, total, 1.0e-15);
		}

	}
};

